/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <models/imputation_hmm.h>
#include <models/buffer_utils.h>

inline
float horizontal_add (const __m256& a)
{
    __m128 vlow = _mm256_castps256_ps128(a);
    __m128 vhigh = _mm256_extractf128_ps(a, 1); // high 128
   vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
   __m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
   __m128 sums = _mm_add_ps(vlow, shuf);
   shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
   sums = _mm_add_ss(sums, shuf);    // (no wasted instructions, and all of them are the 4B minimum)
   return _mm_cvtss_f32(sums);
}

// Lane masks for one byte of Hvar: lane b is all ones when bit (7 - b) is set, matching the state
// order of Hvar.getByte(). Loading a mask costs less than building it from the byte with shifts.
struct hvar_lane_masks {
	alignas(32) uint32_t lanes[256][8];
	constexpr hvar_lane_masks() : lanes() {
		for (int byte = 0 ; byte < 256 ; byte ++)
			for (int b = 0 ; b < 8 ; b ++)
				lanes[byte][b] = ((byte >> (7 - b)) & 1) ? 0xFFFFFFFFu : 0u;
	}
};
static constexpr hvar_lane_masks HVAR_LANE_MASKS;

inline
__m256 hvarLaneMask(const unsigned char byte)
{
	return _mm256_castsi256_ps(_mm256_load_si256(reinterpret_cast<const __m256i *>(HVAR_LANE_MASKS.lanes[byte])));
}

imputation_hmm::imputation_hmm(conditioning_set * _C, const unsigned int _checkpoint_block) : checkpoint_block(_checkpoint_block), alpha_block(1), n_mat(0), resident_block(-1) {
	C = _C;
	modK=0;
	Emissions = aligned_vector32 < float > (2*C->n_tot_sites);
}

imputation_hmm::~imputation_hmm() {
	AlphaCheckpoints.clear();
	AlphaBlock.clear();
	AlphaSum.clear();
	Emissions.clear();
}

void imputation_hmm::resize()
{
	modK = ((C->n_states / 8) + (C->n_states % 8 ? 1 : 0))*8;
	const size_t n_sites = C->polymorphic_sites.size();
	AlphaSum.resize(n_sites, 0.0f);
	AlphaScale.resize(n_sites);
	AlphaOffset.resize(n_sites);
	AlphaAnchor.resize(n_sites);
	MatSite.resize(n_sites);
	MatFact1.resize(n_sites);
	MatFact2.resize(n_sites);
	alpha_block = std::max < size_t > (1, std::min < size_t > (checkpoint_block ? checkpoint_block : n_sites, n_sites));
	resize_discarding_contents(AlphaBlock, (size_t)alpha_block * modK);
	Beta.resize(modK);
	AB.resize(modK);
}

void imputation_hmm::init(const std::vector < float > & HL)
{
	for (int l = 0 ; l < C->n_tot_sites ; l ++)
	{
		float p0 = HL[2*l+0] * C->ee_imp + HL[2*l+1] * C->ed_imp;
		float p1 = HL[2*l+0] * C->ed_imp + HL[2*l+1] * C->ee_imp;
		Emissions[2*l+0] = p0 / (p0+p1);
		Emissions[2*l+1] = p1 / (p0+p1);
	}

}

void imputation_hmm::computePosteriors(const std::vector < float > & HL, std::vector < bool > & flat, std::vector < float > & HP) {
	resize();
	init(HL);
	forward(flat);
	backward(HL, flat, HP);
}

// Computes materialised row m into curr from row m-1 (prev) with the factors forward() folded for it,
// and returns its sum. forward() and recomputeBlock() both use it, so recomputed rows are identical.
float imputation_hmm::materialiseRow(const int m, const std::vector < bool > & flat, const float * prev, float * curr) {
	const unsigned int nstates = C->n_states;
	const unsigned int nstatesMD8 = (nstates / 8) * 8;
	const int l = MatSite[m];
	float sum = 0.0f;

	if (flat[C->polymorphic_sites[l]] || C->lq_flag[C->polymorphic_sites[l]])
	{
		if (m == 0)
		{
			std::fill(curr, curr + modK, 1.0f / nstates);
			sum = 1.0f;
		}
		else
		{
			const float fact1 = MatFact1[m];
			const float fact2 = MatFact2[m];
			const __m256 _fact1 = _mm256_set1_ps(fact1);
			const __m256 _fact2 = _mm256_set1_ps(fact2);
			__m256 _sum = _mm256_set1_ps(0.0f);
			int k = 0;
			for (; k < nstatesMD8; k += 8)
			{
				const __m256 _prob_prev = _mm256_load_ps(&prev[k]);
				const __m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
				_sum = _mm256_add_ps(_sum, _prob_curr);
				_mm256_store_ps(&curr[k], _prob_curr);
			}
			if (k) sum = horizontal_add(_sum);
			for (int offset = nstatesMD8; offset < nstates ; offset ++) {
				curr[offset] = (prev[offset] * fact2 + fact1);
				sum += curr[offset];
			}
		}
	}
	else
	{
		const std::array<float,2> emit = {Emissions[2*C->polymorphic_sites[l]+0], Emissions[2*C->polymorphic_sites[l]+1]};
		const __m256 _emit0 = _mm256_set1_ps(emit[0]);
		const __m256 _emit1 = _mm256_set1_ps(emit[1]);
		if (m == 0)
		{
			const float fact1 = 1.0f / nstates;
			const __m256 _fact1 = _mm256_set1_ps(fact1);
			__m256 _sum = _mm256_set1_ps(0.0f);
			int k = 0;
			for (; k < nstatesMD8; k += 8)
			{
				const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, hvarLaneMask(C->Hvar.getByte(l, k)));
				const __m256 _prob_curr = _mm256_mul_ps(_emiss, _fact1);
				_sum = _mm256_add_ps(_sum, _prob_curr);
				_mm256_store_ps(&curr[k], _prob_curr);
			}
			if (k) sum = horizontal_add(_sum);
			for (int offset = nstatesMD8; offset < nstates ; offset ++)
			{
				curr[offset] = emit[C->Hvar.get(l, offset)] * fact1;
				sum += curr[offset];
			}
		}
		else
		{
			const float fact1 = MatFact1[m];
			const float fact2 = MatFact2[m];
			const __m256 _fact1 = _mm256_set1_ps(fact1);
			const __m256 _fact2 = _mm256_set1_ps(fact2);
			__m256 _sum = _mm256_set1_ps(0.0f);

			int k = 0;
			for (; k < nstatesMD8; k += 8)
			{
				const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, hvarLaneMask(C->Hvar.getByte(l, k)));
				const __m256 _prob_prev = _mm256_load_ps(&prev[k]);
				const __m256 _prob_temp = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
				const __m256 _prob_curr = _mm256_mul_ps(_prob_temp, _emiss);
				_sum = _mm256_add_ps(_sum, _prob_curr);
				_mm256_store_ps(&curr[k], _prob_curr);
			}
			if (k) sum = horizontal_add(_sum);
			for (int offset = nstatesMD8; offset < nstates ; offset ++)
			{
				curr[offset] = (prev[offset]*fact2+fact1)*emit[C->Hvar.get(l, offset)];
				sum += curr[offset];
			}
		}
	}
	return sum;
}

void imputation_hmm::forward(std::vector < bool > & flat) {
	const int n_sites = C->polymorphic_sites.size();
	const unsigned int nstates = C->n_states;
	// Inside a flat run alpha_l = runA * row[n_mat-1] + runB. AlphaSum follows the closed form
	// sum_k(a_k*fact2 + fact1) = fact2*AlphaSum[l-1] + fact1*nstates. The last site is always
	// materialised so that backward() can start from its row.
	double runA = 1.0, runB = 0.0;
	n_mat = 0;

	// Checkpoints for the rows that will be materialised: sites with reads plus the first and last.
	// Underflow re-anchors can add rows; the checkpoints grow for those when they occur.
	size_t n_expected = 0;
	for (int l = 0 ; l < n_sites ; l ++)
		n_expected += (l == 0 || l == n_sites - 1 || !(flat[C->polymorphic_sites[l]] || C->lq_flag[C->polymorphic_sites[l]]));
	resize_discarding_contents(AlphaCheckpoints, ((n_expected + alpha_block - 1) / alpha_block) * modK);

	for (int l = 0 ; l < n_sites ; l ++)
	{
		const bool site_flat = flat[C->polymorphic_sites[l]] || C->lq_flag[C->polymorphic_sites[l]];
		if (l)
		{
			const float fact1 = C->t[l-1] / nstates;
			const float fact2 = C->nt[l-1] / AlphaSum[l-1];
			// When the previous site was materialised (runA=1, runB=0) these are exactly fact2 and fact1.
			const double newA = runA * fact2;
			const double newB = runB * fact2 + fact1;
			if (site_flat && l != n_sites - 1 && newA > 1e-30)
			{
				runA = newA;
				runB = newB;
				AlphaSum[l] = fact2 * AlphaSum[l-1] + fact1 * nstates;
				AlphaScale[l] = (float)runA;
				AlphaOffset[l] = (float)runB;
				AlphaAnchor[l] = n_mat - 1;
				continue;
			}
			// Materialise; the newA guard re-anchors a run whose scale would underflow float.
			MatFact1[n_mat] = (float)newB;
			MatFact2[n_mat] = (float)newA;
		}
		const int m = n_mat ++;
		MatSite[m] = l;
		const float * prev = m ? &AlphaBlock[(size_t)((m-1) % alpha_block) * modK] : nullptr;
		float * curr = &AlphaBlock[(size_t)(m % alpha_block) * modK];
		AlphaSum[l] = materialiseRow(m, flat, prev, curr);
		if (m % alpha_block == 0)
		{
			const size_t offset = (size_t)(m / alpha_block) * modK;
			if (offset + modK > AlphaCheckpoints.size()) AlphaCheckpoints.resize(offset + modK);
			std::copy(curr, curr + modK, AlphaCheckpoints.begin() + offset);
		}
		runA = 1.0;
		runB = 0.0;
		AlphaScale[l] = 1.0f;
		AlphaOffset[l] = 0.0f;
		AlphaAnchor[l] = m;
	}
	resident_block = n_mat ? (int)((n_mat - 1) / alpha_block) : -1;
}

void imputation_hmm::recomputeBlock(const int block, const std::vector < bool > & flat) {
	const int first = block * alpha_block;
	const int last = std::min(first + (int)alpha_block, n_mat) - 1;
	std::copy(AlphaCheckpoints.begin() + (size_t)block * modK, AlphaCheckpoints.begin() + (size_t)(block + 1) * modK, AlphaBlock.begin());
	for (int m = first + 1 ; m <= last ; m ++) materialiseRow(m, flat, &AlphaBlock[(size_t)(m - first - 1) * modK], &AlphaBlock[(size_t)(m - first) * modK]);
	resident_block = block;
}

// Row m of the forward pass. backward() asks for rows in non-increasing order, so each block is
// recomputed at most once.
const float * imputation_hmm::alphaRow(const int m, const std::vector < bool > & flat) {
	if ((int)(m / alpha_block) != resident_block) recomputeBlock(m / alpha_block, flat);
	return &AlphaBlock[(size_t)(m % alpha_block) * modK];
}

void imputation_hmm::backward(const std::vector < float > & HL, std::vector < bool > & flat, std::vector < float > & HP)
{
	float betaSum = 0.0f, betaSumNext = 0.0f;
	std::array<float,2> prob_hid, prob_obs;
	const unsigned int nstates = C->n_states;
	const unsigned int nstatesMD8 = (nstates / 8) * 8;
	const int n_sites = C->polymorphic_sites.size();
	fill(Beta.begin(), Beta.end(), 1.0f);
	__m256 _sum,  _prob0, _prob1;

	// Mirror of forward()'s composition: inside a flat run beta_l = runC * Beta + runD, where Beta
	// holds the row of the nearest materialised site to the right. A flat site's posterior needs
	// sum_{k:H=a} alpha_l*beta_l; expanding both affine maps, it comes from masked sums over AB
	// (anchor row * Beta), the anchor row and Beta, plus the count of ALT states. The REF side is
	// the run totals minus the ALT side. AB and the totals are rebuilt when the anchor or Beta changes.
	double runC = 1.0, runD = 0.0;
	double S1tot = 0.0, S2tot = 0.0, S3tot = 0.0;
	int p_anchor = -1;

	for (int l = n_sites-1 ; l >= 0 ; l --)
	{
		betaSum=0.0f;
		prob_hid[0]=0.0f;
		prob_hid[1]=0.0f;

		_sum = _mm256_set1_ps(0.0f);
		_prob0 = _mm256_set1_ps(0.0f);
		_prob1 = _mm256_set1_ps(0.0f);

		if (flat[C->polymorphic_sites[l]] || C->lq_flag[C->polymorphic_sites[l]])
		{
			if (l == n_sites-1)
			{
				const float * alpha_l = alphaRow(AlphaAnchor[l], flat);
				const float fact1 = 1.0f / nstates;
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256 _mask_curr = hvarLaneMask(C->Hvar.getByte(l, k));
					const __m256 _alphas = _mm256_load_ps(&alpha_l[k]);
					_prob0 = _mm256_add_ps(_prob0, _mm256_andnot_ps(_mask_curr, _alphas));
					_prob1 = _mm256_add_ps(_prob1, _mm256_and_ps(_alphas, _mask_curr));
					_mm256_store_ps(&Beta[k], _fact1);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_prob0);
					prob_hid[1] = horizontal_add(_prob1);
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Beta[offset] = fact1;
					prob_hid[C->Hvar.get(l, offset)] += alpha_l[offset];
				}
				betaSum = 1.0f;
				runC = 1.0; runD = 0.0; p_anchor = -1;
			}
			else
			{
				const float fact1 = C->t[l] / nstates;
				const float fact2 = C->nt[l] / betaSumNext;
				const double newC = runC * fact2;
				const double newD = runD * fact2 + fact1;
				if (newC <= 1e-30)
				{
					// The run's scale would underflow float: materialise Beta with the folded factors
					// and start a new run (rare; only on very long flat runs).
					const __m256 _fact1 = _mm256_set1_ps((float)newD);
					const __m256 _fact2 = _mm256_set1_ps((float)newC);
					int k = 0;
					for (; k < nstatesMD8; k += 8)
					{
						const __m256 _prob_curr = _mm256_fmadd_ps(_mm256_load_ps(&Beta[k]), _fact2, _fact1);
						_sum = _mm256_add_ps(_sum, _prob_curr);
						_mm256_store_ps(&Beta[k], _prob_curr);
					}
					if (k) betaSum = horizontal_add(_sum);
					for (int offset = nstatesMD8; offset < nstates ; offset ++)
					{
						Beta[offset] = Beta[offset] * (float)newC + (float)newD;
						betaSum += Beta[offset];
					}
					runC = 1.0; runD = 0.0; p_anchor = -1;
				}
				else
				{
					runC = newC;
					runD = newD;
					betaSum = fact2 * betaSumNext + fact1 * nstates;
				}

				if (p_anchor != AlphaAnchor[l])
				{
					p_anchor = AlphaAnchor[l];
					const float * anchor = alphaRow(p_anchor, flat);
					__m256 _s1 = _mm256_setzero_ps(), _s2 = _mm256_setzero_ps(), _s3 = _mm256_setzero_ps();
					int k = 0;
					for (; k < nstatesMD8; k += 8)
					{
						const __m256 _a = _mm256_load_ps(&anchor[k]);
						const __m256 _b = _mm256_load_ps(&Beta[k]);
						const __m256 _ab = _mm256_mul_ps(_a, _b);
						_mm256_store_ps(&AB[k], _ab);
						_s1 = _mm256_add_ps(_s1, _ab);
						_s2 = _mm256_add_ps(_s2, _a);
						_s3 = _mm256_add_ps(_s3, _b);
					}
					S1tot = k ? horizontal_add(_s1) : 0.0f;
					S2tot = k ? horizontal_add(_s2) : 0.0f;
					S3tot = k ? horizontal_add(_s3) : 0.0f;
					for (int offset = nstatesMD8; offset < nstates ; offset ++)
					{
						AB[offset] = anchor[offset] * Beta[offset];
						S1tot += AB[offset];
						S2tot += anchor[offset];
						S3tot += Beta[offset];
					}
				}

				// Each allele's share of sum_k alpha_l*beta_l, expanded over the two affine maps, is
				// AC*sum(AB) + AD*sum(anchor) + BC*sum(Beta) + BD*count over the states carrying it. Both
				// sides keep relative precision, as on the base: at low-quality sites with reads the
				// smaller side is later scaled by read likelihoods that can favour it by up to 1e10.
				const int abs_site = C->polymorphic_sites[l];
				const float * anchor = alphaRow(p_anchor, flat);
				const double A = AlphaScale[l], B = AlphaOffset[l];
				const double AC = A*runC, AD = A*runD, BC = B*runC, BD = B*runD;
				// Sums both sides directly in one masked pass over the per-state weights
				// w = AC*AB + AD*anchor + BC*Beta. The posterior depends only on the ratio of the two
				// sides, so the coefficients are scaled to a maximum of 1 to keep w well inside float range.
				auto directSides = [&]() {
					const double cmax = std::max(std::max(AC, AD), std::max(BC, BD));
					const float c1 = AC / cmax, c2 = AD / cmax, c3 = BC / cmax, c4 = BD / cmax;
					const __m256 _c1 = _mm256_set1_ps(c1), _c2 = _mm256_set1_ps(c2), _c3 = _mm256_set1_ps(c3);
					__m256 _alt = _mm256_setzero_ps(), _ref = _mm256_setzero_ps();
					int cnt1 = 0;
					int k = 0;
					for (; k < nstatesMD8; k += 8)
					{
						const unsigned int byte = C->Hvar.getByte(l, k);
						cnt1 += __builtin_popcount(byte);
						const __m256 _f = hvarLaneMask(byte);
						const __m256 _w = _mm256_fmadd_ps(_mm256_load_ps(&AB[k]), _c1, _mm256_fmadd_ps(_mm256_load_ps(&anchor[k]), _c2, _mm256_mul_ps(_mm256_load_ps(&Beta[k]), _c3)));
						_alt = _mm256_add_ps(_alt, _mm256_and_ps(_w, _f));
						_ref = _mm256_add_ps(_ref, _mm256_andnot_ps(_f, _w));
					}
					double alt = k ? horizontal_add(_alt) : 0.0f;
					double ref = k ? horizontal_add(_ref) : 0.0f;
					for (int offset = nstatesMD8; offset < nstates ; offset ++)
					{
						const float w = AB[offset] * c1 + anchor[offset] * c2 + Beta[offset] * c3;
						if (C->Hvar.get(l, offset)) { alt += w; cnt1 ++; }
						else ref += w;
					}
					prob_hid[1] = (float)(cmax * (alt + c4 * cnt1));
					prob_hid[0] = (float)(cmax * (ref + c4 * (nstates - cnt1)));
				};
				if (C->var_type[abs_site] == TYPE_RARE)
				{
					// A rare site carries its major allele in every state except the few listed in Svar.
					// Sum the minor side over those, and take the major side from the totals unless it
					// holds so little of the mass that rounding in the totals would swamp it.
					const std::vector < unsigned int > & minor = C->Svar[abs_site];
					double M1 = 0.0, M2 = 0.0, M3 = 0.0;
					for (const unsigned int k : minor)
					{
						M1 += AB[k];
						M2 += anchor[k];
						M3 += Beta[k];
					}
					const double n_minor = minor.size();
					const double minor_side = AC*M1 + AD*M2 + BC*M3 + BD*n_minor;
					const double major_side = AC*(S1tot-M1) + AD*(S2tot-M2) + BC*(S3tot-M3) + BD*(nstates-n_minor);
					if (major_side >= 1e-2 * (major_side + minor_side))
					{
						const bool major = C->major_alleles[abs_site];
						prob_hid[major] = (float)major_side;
						prob_hid[!major] = (float)minor_side;
					}
					else directSides();
				}
				else directSides();
			}

			prob_obs[0] = (prob_hid[0]*C->ee_imp + prob_hid[1]*C->ed_imp);
			prob_obs[1] = (prob_hid[0]*C->ed_imp + prob_hid[1]*C->ee_imp);
			if (!flat[C->polymorphic_sites[l]])
			{
				prob_obs[0] *= HL[2*C->polymorphic_sites[l]+0];
				prob_obs[1] *= HL[2*C->polymorphic_sites[l]+1];
			}
		}
		else
		{
			const float * alpha_l = alphaRow(AlphaAnchor[l], flat);
			std::array<float,2> emit = {Emissions[2*C->polymorphic_sites[l]+0], Emissions[2*C->polymorphic_sites[l]+1]};
			const __m256 _emit0 = _mm256_set1_ps(emit[0]);
			const __m256 _emit1 = _mm256_set1_ps(emit[1]);

			if (l == n_sites-1)
			{
				const float fact1 = 1.0f / nstates;
				const __m256 _fact1 = _mm256_set1_ps(fact1);

				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256 _mask_curr = hvarLaneMask(C->Hvar.getByte(l, k));
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mask_curr);
					const __m256 _alphas = _mm256_load_ps(&alpha_l[k]);
					const __m256 _prob_next = _mm256_mul_ps(_emiss, _fact1);
					_prob0 = _mm256_add_ps(_prob0, _mm256_andnot_ps(_mask_curr, _alphas));
					_prob1 = _mm256_add_ps(_prob1, _mm256_and_ps(_alphas, _mask_curr));
					_sum= _mm256_add_ps(_sum, _prob_next);
					_mm256_store_ps(&Beta[k], _prob_next);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_prob0);
					prob_hid[1] = horizontal_add(_prob1);
					betaSum = horizontal_add(_sum);
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					prob_hid[C->Hvar.get(l, offset)] += alpha_l[offset];
					Beta[offset] = emit[C->Hvar.get(l, offset)] * fact1;
					betaSum += Beta[offset];
				}
			}
			else
			{
				// Fold any pending flat run into this site's factors; with no run (runC=1, runD=0)
				// these are exactly fact2 and fact1.
				const float fact1 = (float)(runD * (C->nt[l] / betaSumNext) + C->t[l] / nstates);
				const float fact2 = (float)(runC * (C->nt[l] / betaSumNext));
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				const __m256 _fact2 = _mm256_set1_ps(fact2);

				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256 _mask_curr = hvarLaneMask(C->Hvar.getByte(l, k));
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mask_curr);
					const __m256 _prob_prev = _mm256_load_ps(&Beta[k]);
					const __m256 _alphas = _mm256_load_ps(&alpha_l[k]);
					const __m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
					const __m256 _dotprod = _mm256_mul_ps(_alphas, _prob_curr);
					const __m256 _prob_next = _mm256_mul_ps(_prob_curr, _emiss);
					_prob0 = _mm256_add_ps(_prob0, _mm256_andnot_ps(_mask_curr, _dotprod));
					_prob1 = _mm256_add_ps(_prob1, _mm256_and_ps(_dotprod, _mask_curr));
					_sum= _mm256_add_ps(_sum, _prob_next);
					_mm256_store_ps(&Beta[k], _prob_next);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_prob0);
					prob_hid[1] = horizontal_add(_prob1);
					betaSum = horizontal_add(_sum);
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Beta[offset] = Beta[offset] * fact2 + fact1;
					prob_hid[C->Hvar.get(l, offset)] += alpha_l[offset] * Beta[offset];
					Beta[offset] *= emit[C->Hvar.get(l, offset)];
					betaSum += Beta[offset];
				}
			}
			runC = 1.0; runD = 0.0; p_anchor = -1;
			prob_hid[0] /= emit[0];
			prob_hid[1] /= emit[1];
			prob_obs[0] = (prob_hid[0]*C->ee_imp + prob_hid[1]*C->ed_imp) * HL[2*C->polymorphic_sites[l]+0];
			prob_obs[1] = (prob_hid[0]*C->ed_imp + prob_hid[1]*C->ee_imp) * HL[2*C->polymorphic_sites[l]+1];

		}
		HP[2*C->polymorphic_sites[l]+0] = prob_obs[0] / (prob_obs[0] + prob_obs[1]);
		HP[2*C->polymorphic_sites[l]+1] = prob_obs[1] / (prob_obs[0] + prob_obs[1]);
		betaSumNext = betaSum;
	}
	// Monomorphic sites
	for (int l = 0 ; l < C->monomorphic_sites.size() ; l ++)
	{
		prob_obs[C->major_alleles[C->monomorphic_sites[l]]] = C->ee_imp;
		prob_obs[!C->major_alleles[C->monomorphic_sites[l]]] = C->ed_imp;
		if (!flat[C->monomorphic_sites[l]])
		{
			prob_obs[0] *= HL[2*C->monomorphic_sites[l]+0];
			prob_obs[1] *= HL[2*C->monomorphic_sites[l]+1];
		}
		HP[2*C->monomorphic_sites[l]+0] = prob_obs[0] / (prob_obs[0] + prob_obs[1]);
		HP[2*C->monomorphic_sites[l]+1] = prob_obs[1] / (prob_obs[0] + prob_obs[1]);
	}
}

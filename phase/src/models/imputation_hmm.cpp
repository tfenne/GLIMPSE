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

imputation_hmm::imputation_hmm(conditioning_set * _C) {
	C = _C;
	modK=0;
	Emissions = aligned_vector32 < float > (2*C->n_tot_sites);
}

imputation_hmm::~imputation_hmm() {
	Alpha.clear();
	AlphaSum.clear();
	Emissions.clear();
}

void imputation_hmm::resize()
{
	modK = ((C->n_states / 8) + (C->n_states % 8 ? 1 : 0))*8;
	AlphaSum.resize(C->polymorphic_sites.size(), 0.0f);
	Alpha.resize(C->polymorphic_sites.size() * modK, 0.0f);
	Beta.resize(modK);
	AlphaScale.resize(C->polymorphic_sites.size());
	AlphaOffset.resize(C->polymorphic_sites.size());
	AlphaAnchor.resize(C->polymorphic_sites.size());
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

void imputation_hmm::forward(std::vector < bool > & flat) {
	const __m256i _vshift_count = _mm256_set_epi32(31,30,29,28,27,26,25,24);
	const unsigned int nstates = C->n_states;
	const unsigned int nstatesMD8 = (nstates / 8) * 8;
	const int n_poly = C->polymorphic_sites.size();

	//Flat sites have no emission term, so their update is alpha_l = alpha_{l-1}*fact2 + fact1
	//with scalar factors: runs of flat sites compose into a single affine map of the last
	//materialized ("anchor") row, and their rows are never computed or stored. AlphaSum
	//follows the closed form sum_k(a_k*fact2 + fact1) = fact2*AlphaSum_prev + fact1*nstates.
	//The run scalars are tracked relative to the anchor row for backward() via
	//AlphaScale/AlphaOffset/AlphaAnchor. The last site is always materialized so that
	//backward()'s initialization can read it directly.
	int anchor = 0;
	double runA = 1.0, runB = 0.0;

	for (int l = 0 ; l < n_poly ; l ++)
	{
		AlphaSum[l] = 0.0f;
		const int abs_site = C->polymorphic_sites[l];
		const bool site_flat = flat[abs_site] || C->lq_flag[abs_site];

		if (l == 0)
		{
			if (site_flat)
			{
				fill(Alpha.begin(), Alpha.begin()+modK, 1.0f / nstates);
				AlphaSum[l] = 1.0f;
			}
			else
			{
				const std::array<float,2> emit = {Emissions[2*abs_site+0], Emissions[2*abs_site+1]};
				const __m256 _emit0 = _mm256_set1_ps(emit[0]);
				const __m256 _emit1 = _mm256_set1_ps(emit[1]);
				const float fact1 = 1.0f / nstates;
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				//Independent accumulators to break the _sum latency chain.
				__m256 _sum0 = _mm256_setzero_ps(), _sum1 = _mm256_setzero_ps();
				int k = 0;
				for (; k + 16 <= nstatesMD8; k += 16)
				{
					const __m256i _mask0 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+0)), _vshift_count);
					const __m256i _mask1 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+8)), _vshift_count);
					const __m256 _p0 = _mm256_mul_ps(_mm256_blendv_ps(_emit0, _emit1, _mm256_castsi256_ps(_mask0)), _fact1);
					const __m256 _p1 = _mm256_mul_ps(_mm256_blendv_ps(_emit0, _emit1, _mm256_castsi256_ps(_mask1)), _fact1);
					_sum0 = _mm256_add_ps(_sum0, _p0);
					_sum1 = _mm256_add_ps(_sum1, _p1);
					_mm256_store_ps(&Alpha[l*modK+k+0], _p0);
					_mm256_store_ps(&Alpha[l*modK+k+8], _p1);
				}
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _bcst = _mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k));
					const __m256i _mask = _mm256_sllv_epi32(_bcst, _vshift_count);
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask));
					const __m256 _prob_curr = _mm256_mul_ps(_emiss, _fact1);
					_sum0 = _mm256_add_ps(_sum0, _prob_curr);
					_mm256_store_ps(&Alpha[l*modK+k], _prob_curr);
				}
				if (k) AlphaSum[l] = horizontal_add(_mm256_add_ps(_sum0, _sum1));
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Alpha[l*modK+offset] = emit[C->Hvar.get(l, offset)] * fact1;
					AlphaSum[l] += Alpha[l*modK+offset];
				}
			}
			anchor = 0; runA = 1.0; runB = 0.0;
			AlphaAnchor[l] = 0; AlphaScale[l] = 1.0f; AlphaOffset[l] = 0.0f;
			continue;
		}

		const float fact1 = C->t[l-1] / nstates;
		const float fact2 = C->nt[l-1] / AlphaSum[l-1];
		//Fold the accumulated run into this site's transition; when the previous site was
		//materialized (runA=1, runB=0) these are exactly fact2 and fact1.
		const double newA = runA * fact2;
		const double newB = runB * fact2 + fact1;

		if (site_flat && l != n_poly-1 && newA > 1e-30)
		{
			runA = newA; runB = newB;
			AlphaSum[l] = fact2 * AlphaSum[l-1] + fact1 * nstates;
			AlphaScale[l] = (float)runA; AlphaOffset[l] = (float)runB;
			AlphaAnchor[l] = anchor;
			continue;
		}

		//Materialize row l from the anchor row with the folded factors. The newA guard
		//above re-materializes mid-run in the (rare) case the run scale underflows float.
		const __m256 _fact1 = _mm256_set1_ps((float)newB);
		const __m256 _fact2 = _mm256_set1_ps((float)newA);
		const float sfact1 = (float)newB;
		const float sfact2 = (float)newA;
		if (site_flat)
		{
			//Four independent accumulators: a single _sum chains one add per 8 states,
			//and its latency (plus the horizontal add and divide feeding the next site's
			//factors) dominates the loop, which otherwise runs well below port capacity.
			__m256 _sum0 = _mm256_setzero_ps(), _sum1 = _mm256_setzero_ps();
			__m256 _sum2 = _mm256_setzero_ps(), _sum3 = _mm256_setzero_ps();
			int k = 0;
			for (; k + 32 <= nstatesMD8; k += 32)
			{
				const __m256 _p0 = _mm256_fmadd_ps(_mm256_load_ps(&Alpha[anchor*modK+k+0]), _fact2, _fact1);
				const __m256 _p1 = _mm256_fmadd_ps(_mm256_load_ps(&Alpha[anchor*modK+k+8]), _fact2, _fact1);
				const __m256 _p2 = _mm256_fmadd_ps(_mm256_load_ps(&Alpha[anchor*modK+k+16]), _fact2, _fact1);
				const __m256 _p3 = _mm256_fmadd_ps(_mm256_load_ps(&Alpha[anchor*modK+k+24]), _fact2, _fact1);
				_sum0 = _mm256_add_ps(_sum0, _p0);
				_sum1 = _mm256_add_ps(_sum1, _p1);
				_sum2 = _mm256_add_ps(_sum2, _p2);
				_sum3 = _mm256_add_ps(_sum3, _p3);
				_mm256_store_ps(&Alpha[l*modK+k+0], _p0);
				_mm256_store_ps(&Alpha[l*modK+k+8], _p1);
				_mm256_store_ps(&Alpha[l*modK+k+16], _p2);
				_mm256_store_ps(&Alpha[l*modK+k+24], _p3);
			}
			for (; k < nstatesMD8; k += 8)
			{
				const __m256 _prob_prev = _mm256_load_ps(&Alpha[anchor*modK+k]);
				const __m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
				_sum0 = _mm256_add_ps(_sum0, _prob_curr);
				_mm256_store_ps(&Alpha[l*modK+k], _prob_curr);
			}
			if (k) AlphaSum[l] = horizontal_add(_mm256_add_ps(_mm256_add_ps(_sum0, _sum1), _mm256_add_ps(_sum2, _sum3)));
			for (int offset = nstatesMD8; offset < nstates ; offset ++) {
				Alpha[l*modK+offset] = (Alpha[anchor*modK+offset] * sfact2 + sfact1);
				AlphaSum[l] += Alpha[l*modK+offset];
			}
		}
		else
		{
			const std::array<float,2> emit = {Emissions[2*abs_site+0], Emissions[2*abs_site+1]};
			const __m256 _emit0 = _mm256_set1_ps(emit[0]);
			const __m256 _emit1 = _mm256_set1_ps(emit[1]);
			//Independent accumulators to break the _sum latency chain.
			__m256 _sum0 = _mm256_setzero_ps(), _sum1 = _mm256_setzero_ps();
			int k = 0;
			for (; k + 16 <= nstatesMD8; k += 16)
			{
				const __m256i _mask0 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+0)), _vshift_count);
				const __m256i _mask1 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+8)), _vshift_count);
				const __m256 _p0 = _mm256_mul_ps(_mm256_fmadd_ps(_mm256_load_ps(&Alpha[anchor*modK+k+0]), _fact2, _fact1), _mm256_blendv_ps(_emit0, _emit1, _mm256_castsi256_ps(_mask0)));
				const __m256 _p1 = _mm256_mul_ps(_mm256_fmadd_ps(_mm256_load_ps(&Alpha[anchor*modK+k+8]), _fact2, _fact1), _mm256_blendv_ps(_emit0, _emit1, _mm256_castsi256_ps(_mask1)));
				_sum0 = _mm256_add_ps(_sum0, _p0);
				_sum1 = _mm256_add_ps(_sum1, _p1);
				_mm256_store_ps(&Alpha[l*modK+k+0], _p0);
				_mm256_store_ps(&Alpha[l*modK+k+8], _p1);
			}
			for (; k < nstatesMD8; k += 8)
			{
				const __m256i _mask = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
				const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask));
				const __m256 _prob_prev = _mm256_load_ps(&Alpha[anchor*modK+k]);
				const __m256 _prob_temp = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
				const __m256 _prob_curr = _mm256_mul_ps(_prob_temp, _emiss);
				_sum0 = _mm256_add_ps(_sum0, _prob_curr);
				_mm256_store_ps(&Alpha[l*modK+k], _prob_curr);
			}
			if (k) AlphaSum[l] = horizontal_add(_mm256_add_ps(_sum0, _sum1));
			for (int offset = nstatesMD8; offset < nstates ; offset ++)
			{
				Alpha[l*modK+offset] = (Alpha[anchor*modK+offset]*sfact2+sfact1)*emit[C->Hvar.get(l, offset)];
				AlphaSum[l] += Alpha[l*modK+offset];
			}
		}
		anchor = l; runA = 1.0; runB = 0.0;
		AlphaAnchor[l] = l; AlphaScale[l] = 1.0f; AlphaOffset[l] = 0.0f;
	}
}

void imputation_hmm::backward(const std::vector < float > & HL, std::vector < bool > & flat, std::vector < float > & HP)
{
	float betaSum = 0.0f, betaSumNext = 0.0f;
	std::array<float,2> prob_hid, prob_obs;
	const __m256i _vshift_count = _mm256_set_epi32(31,30,29,28,27,26,25,24);
	const unsigned int nstates = C->n_states;
	const unsigned int nstatesMD8 = (nstates / 8) * 8;
	const __m256 _zero = _mm256_set1_ps(0.0f);
	const __m256 _one = _mm256_set1_ps(1.0f);
	const int n_poly = C->polymorphic_sites.size();
	fill(Beta.begin(), Beta.end(), 1.0f);

	//Mirror of forward()'s flat-run composition: inside a run of flat sites
	//beta_l = runC * Beta + runD, where Beta still holds the row materialized at the
	//nearest processed site to the right, so those rows are never recomputed or stored.
	//The per-site posterior needs the two masked reductions sum_{k:H=a} alpha_l*beta_l;
	//expanding both affine maps, it is recovered from masked sums over three fixed rows
	//(AB = anchorAlpha*Beta elementwise, the anchor Alpha row, and the Beta row) plus
	//the count of alt-carrying states, with the complement obtained from run totals.
	double runC = 1.0, runD = 0.0;
	double S1tot = 0.0, S2tot = 0.0, S3tot = 0.0;
	int p_anchor = -1;

	for (int l = n_poly-1 ; l >= 0 ; l --)
	{
		betaSum=0.0f;
		prob_hid[0]=0.0f;
		prob_hid[1]=0.0f;
		const int abs_site = C->polymorphic_sites[l];

		if (flat[abs_site] || C->lq_flag[abs_site])
		{
			if (l == n_poly-1)
			{
				const float fact1 = 1.0f / nstates;
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				//Two independent accumulators per reduction: a single accumulator chains one
				//add per 8 states and its latency dominates the loop (same reasoning as forward()).
				__m256 _prob0a = _mm256_setzero_ps(), _prob0b = _mm256_setzero_ps();
				__m256 _prob1a = _mm256_setzero_ps(), _prob1b = _mm256_setzero_ps();
				int k = 0;
				for (; k + 16 <= nstatesMD8; k += 16)
				{
					const __m256i _mc0 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+0)), _vshift_count);
					const __m256i _mc1 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+8)), _vshift_count);
					const __m256 _a0 = _mm256_load_ps(&Alpha[l*modK+k+0]);
					const __m256 _a1 = _mm256_load_ps(&Alpha[l*modK+k+8]);
					_prob0a = _mm256_add_ps(_prob0a, _mm256_mul_ps(_a0, _mm256_blendv_ps(_one, _zero, _mm256_castsi256_ps(_mc0))));
					_prob1a = _mm256_add_ps(_prob1a, _mm256_mul_ps(_a0, _mm256_blendv_ps(_zero, _one, _mm256_castsi256_ps(_mc0))));
					_prob0b = _mm256_add_ps(_prob0b, _mm256_mul_ps(_a1, _mm256_blendv_ps(_one, _zero, _mm256_castsi256_ps(_mc1))));
					_prob1b = _mm256_add_ps(_prob1b, _mm256_mul_ps(_a1, _mm256_blendv_ps(_zero, _one, _mm256_castsi256_ps(_mc1))));
					_mm256_store_ps(&Beta[k+0], _fact1);
					_mm256_store_ps(&Beta[k+8], _fact1);
				}
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _mask_curr = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
					const __m256 _mask0 = _mm256_blendv_ps (_one, _zero, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask1 = _mm256_blendv_ps (_zero, _one, _mm256_castsi256_ps(_mask_curr));
					const __m256 _alphas = _mm256_load_ps(&Alpha[l*modK+k]);
					_prob0a = _mm256_add_ps(_prob0a, _mm256_mul_ps(_alphas, _mask0));
					_prob1a = _mm256_add_ps(_prob1a, _mm256_mul_ps(_alphas, _mask1));
					_mm256_store_ps(&Beta[k], _fact1);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_mm256_add_ps(_prob0a, _prob0b));
					prob_hid[1] = horizontal_add(_mm256_add_ps(_prob1a, _prob1b));
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Beta[offset] = fact1;
					prob_hid[C->Hvar.get(l, offset)] += Alpha[l*modK+offset];
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
					//Run scale underflows float: re-materialize the Beta row with the folded
					//factors and re-anchor the run (rare; only on very long flat runs).
					const __m256 _fact1 = _mm256_set1_ps((float)newD);
					const __m256 _fact2 = _mm256_set1_ps((float)newC);
					__m256 _sum0 = _mm256_setzero_ps(), _sum1 = _mm256_setzero_ps();
					int k = 0;
					for (; k + 16 <= nstatesMD8; k += 16)
					{
						const __m256 _b0 = _mm256_fmadd_ps(_mm256_load_ps(&Beta[k+0]), _fact2, _fact1);
						const __m256 _b1 = _mm256_fmadd_ps(_mm256_load_ps(&Beta[k+8]), _fact2, _fact1);
						_sum0 = _mm256_add_ps(_sum0, _b0);
						_sum1 = _mm256_add_ps(_sum1, _b1);
						_mm256_store_ps(&Beta[k+0], _b0);
						_mm256_store_ps(&Beta[k+8], _b1);
					}
					for (; k < nstatesMD8; k += 8)
					{
						const __m256 _b0 = _mm256_fmadd_ps(_mm256_load_ps(&Beta[k]), _fact2, _fact1);
						_sum0 = _mm256_add_ps(_sum0, _b0);
						_mm256_store_ps(&Beta[k], _b0);
					}
					if (k) betaSum = horizontal_add(_mm256_add_ps(_sum0, _sum1));
					for (int offset = nstatesMD8; offset < nstates ; offset ++)
					{
						Beta[offset] = Beta[offset] * (float)newC + (float)newD;
						betaSum += Beta[offset];
					}
					runC = 1.0; runD = 0.0; p_anchor = -1;
				}
				else
				{
					runC = newC; runD = newD;
					betaSum = fact2 * betaSumNext + fact1 * nstates;
				}

				if (p_anchor != AlphaAnchor[l])
				{
					//New (anchor Alpha row, Beta row) pair: build AB and the run totals.
					p_anchor = AlphaAnchor[l];
					__m256 _s1a = _mm256_setzero_ps(), _s1b = _mm256_setzero_ps();
					__m256 _s2a = _mm256_setzero_ps(), _s2b = _mm256_setzero_ps();
					__m256 _s3a = _mm256_setzero_ps(), _s3b = _mm256_setzero_ps();
					int k = 0;
					for (; k + 16 <= nstatesMD8; k += 16)
					{
						const __m256 _a0 = _mm256_load_ps(&Alpha[p_anchor*modK+k+0]);
						const __m256 _a1 = _mm256_load_ps(&Alpha[p_anchor*modK+k+8]);
						const __m256 _b0 = _mm256_load_ps(&Beta[k+0]);
						const __m256 _b1 = _mm256_load_ps(&Beta[k+8]);
						const __m256 _p0 = _mm256_mul_ps(_a0, _b0);
						const __m256 _p1 = _mm256_mul_ps(_a1, _b1);
						_mm256_store_ps(&AB[k+0], _p0);
						_mm256_store_ps(&AB[k+8], _p1);
						_s1a = _mm256_add_ps(_s1a, _p0); _s1b = _mm256_add_ps(_s1b, _p1);
						_s2a = _mm256_add_ps(_s2a, _a0); _s2b = _mm256_add_ps(_s2b, _a1);
						_s3a = _mm256_add_ps(_s3a, _b0); _s3b = _mm256_add_ps(_s3b, _b1);
					}
					for (; k < nstatesMD8; k += 8)
					{
						const __m256 _a0 = _mm256_load_ps(&Alpha[p_anchor*modK+k]);
						const __m256 _b0 = _mm256_load_ps(&Beta[k]);
						const __m256 _p0 = _mm256_mul_ps(_a0, _b0);
						_mm256_store_ps(&AB[k], _p0);
						_s1a = _mm256_add_ps(_s1a, _p0);
						_s2a = _mm256_add_ps(_s2a, _a0);
						_s3a = _mm256_add_ps(_s3a, _b0);
					}
					S1tot = horizontal_add(_mm256_add_ps(_s1a, _s1b));
					S2tot = horizontal_add(_mm256_add_ps(_s2a, _s2b));
					S3tot = horizontal_add(_mm256_add_ps(_s3a, _s3b));
					for (int offset = nstatesMD8; offset < nstates ; offset ++)
					{
						const float a = Alpha[p_anchor*modK+offset];
						const float b = Beta[offset];
						AB[offset] = a * b;
						S1tot += a * b; S2tot += a; S3tot += b;
					}
				}

				//Masked reductions over the three cached rows for the alt-carrying states.
				__m256 _m1a = _mm256_setzero_ps(), _m1b = _mm256_setzero_ps();
				__m256 _m2a = _mm256_setzero_ps(), _m2b = _mm256_setzero_ps();
				__m256 _m3a = _mm256_setzero_ps(), _m3b = _mm256_setzero_ps();
				int cnt1 = 0;
				int k = 0;
				for (; k + 16 <= nstatesMD8; k += 16)
				{
					const unsigned int byte0 = C->Hvar.getByte(l, k+0);
					const unsigned int byte1 = C->Hvar.getByte(l, k+8);
					cnt1 += __builtin_popcount(byte0) + __builtin_popcount(byte1);
					const __m256 _f0 = _mm256_castsi256_ps(_mm256_srai_epi32(_mm256_sllv_epi32(_mm256_set1_epi32(byte0), _vshift_count), 31));
					const __m256 _f1 = _mm256_castsi256_ps(_mm256_srai_epi32(_mm256_sllv_epi32(_mm256_set1_epi32(byte1), _vshift_count), 31));
					_m1a = _mm256_add_ps(_m1a, _mm256_and_ps(_mm256_load_ps(&AB[k+0]), _f0));
					_m1b = _mm256_add_ps(_m1b, _mm256_and_ps(_mm256_load_ps(&AB[k+8]), _f1));
					_m2a = _mm256_add_ps(_m2a, _mm256_and_ps(_mm256_load_ps(&Alpha[p_anchor*modK+k+0]), _f0));
					_m2b = _mm256_add_ps(_m2b, _mm256_and_ps(_mm256_load_ps(&Alpha[p_anchor*modK+k+8]), _f1));
					_m3a = _mm256_add_ps(_m3a, _mm256_and_ps(_mm256_load_ps(&Beta[k+0]), _f0));
					_m3b = _mm256_add_ps(_m3b, _mm256_and_ps(_mm256_load_ps(&Beta[k+8]), _f1));
				}
				for (; k < nstatesMD8; k += 8)
				{
					const unsigned int byte0 = C->Hvar.getByte(l, k);
					cnt1 += __builtin_popcount(byte0);
					const __m256 _f0 = _mm256_castsi256_ps(_mm256_srai_epi32(_mm256_sllv_epi32(_mm256_set1_epi32(byte0), _vshift_count), 31));
					_m1a = _mm256_add_ps(_m1a, _mm256_and_ps(_mm256_load_ps(&AB[k]), _f0));
					_m2a = _mm256_add_ps(_m2a, _mm256_and_ps(_mm256_load_ps(&Alpha[p_anchor*modK+k]), _f0));
					_m3a = _mm256_add_ps(_m3a, _mm256_and_ps(_mm256_load_ps(&Beta[k]), _f0));
				}
				double S1m = horizontal_add(_mm256_add_ps(_m1a, _m1b));
				double S2m = horizontal_add(_mm256_add_ps(_m2a, _m2b));
				double S3m = horizontal_add(_mm256_add_ps(_m3a, _m3b));
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					if (C->Hvar.get(l, offset))
					{
						S1m += AB[offset];
						S2m += Alpha[p_anchor*modK+offset];
						S3m += Beta[offset];
						cnt1 ++;
					}
				}
				const double A = AlphaScale[l];
				const double B = AlphaOffset[l];
				const double AC = A*runC, AD = A*runD, BC = B*runC, BD = B*runD;
				prob_hid[1] = (float)(AC*S1m + AD*S2m + BC*S3m + BD*cnt1);
				prob_hid[0] = (float)(AC*(S1tot-S1m) + AD*(S2tot-S2m) + BC*(S3tot-S3m) + BD*(nstates-cnt1));
			}

			prob_obs[0] = (prob_hid[0]*C->ee_imp + prob_hid[1]*C->ed_imp);
			prob_obs[1] = (prob_hid[0]*C->ed_imp + prob_hid[1]*C->ee_imp);
			if (!flat[abs_site])
			{
				prob_obs[0] *= HL[2*abs_site+0];
				prob_obs[1] *= HL[2*abs_site+1];
			}
		}
		else
		{
			std::array<float,2> emit = {Emissions[2*abs_site+0], Emissions[2*abs_site+1]};
			const __m256 _emit0 = _mm256_set1_ps(emit[0]);
			const __m256 _emit1 = _mm256_set1_ps(emit[1]);

			if (l == n_poly-1)
			{
				const float fact1 = 1.0f / nstates;
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				//Independent accumulators per reduction to break the add latency chains.
				__m256 _prob0a = _mm256_setzero_ps(), _prob0b = _mm256_setzero_ps();
				__m256 _prob1a = _mm256_setzero_ps(), _prob1b = _mm256_setzero_ps();
				__m256 _suma = _mm256_setzero_ps(), _sumb = _mm256_setzero_ps();

				int k = 0;
				for (; k + 16 <= nstatesMD8; k += 16)
				{
					const __m256i _mc0 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+0)), _vshift_count);
					const __m256i _mc1 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+8)), _vshift_count);
					const __m256 _a0 = _mm256_load_ps(&Alpha[l*modK+k+0]);
					const __m256 _a1 = _mm256_load_ps(&Alpha[l*modK+k+8]);
					const __m256 _pn0 = _mm256_mul_ps(_mm256_blendv_ps(_emit0, _emit1, _mm256_castsi256_ps(_mc0)), _fact1);
					const __m256 _pn1 = _mm256_mul_ps(_mm256_blendv_ps(_emit0, _emit1, _mm256_castsi256_ps(_mc1)), _fact1);
					_prob0a = _mm256_add_ps(_prob0a, _mm256_mul_ps(_a0, _mm256_blendv_ps(_one, _zero, _mm256_castsi256_ps(_mc0))));
					_prob1a = _mm256_add_ps(_prob1a, _mm256_mul_ps(_a0, _mm256_blendv_ps(_zero, _one, _mm256_castsi256_ps(_mc0))));
					_prob0b = _mm256_add_ps(_prob0b, _mm256_mul_ps(_a1, _mm256_blendv_ps(_one, _zero, _mm256_castsi256_ps(_mc1))));
					_prob1b = _mm256_add_ps(_prob1b, _mm256_mul_ps(_a1, _mm256_blendv_ps(_zero, _one, _mm256_castsi256_ps(_mc1))));
					_suma = _mm256_add_ps(_suma, _pn0);
					_sumb = _mm256_add_ps(_sumb, _pn1);
					_mm256_store_ps(&Beta[k+0], _pn0);
					_mm256_store_ps(&Beta[k+8], _pn1);
				}
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _mask_curr = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask0 = _mm256_blendv_ps (_one, _zero, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask1 = _mm256_blendv_ps (_zero, _one, _mm256_castsi256_ps(_mask_curr));
					const __m256 _alphas = _mm256_load_ps(&Alpha[l*modK+k]);
					const __m256 _prob_next = _mm256_mul_ps(_emiss, _fact1);
					_prob0a = _mm256_add_ps(_prob0a, _mm256_mul_ps(_alphas, _mask0));
					_prob1a = _mm256_add_ps(_prob1a, _mm256_mul_ps(_alphas, _mask1));
					_suma = _mm256_add_ps(_suma, _prob_next);
					_mm256_store_ps(&Beta[k], _prob_next);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_mm256_add_ps(_prob0a, _prob0b));
					prob_hid[1] = horizontal_add(_mm256_add_ps(_prob1a, _prob1b));
					betaSum = horizontal_add(_mm256_add_ps(_suma, _sumb));
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					prob_hid[C->Hvar.get(l, offset)] += Alpha[l*modK+offset];
					Beta[offset] = emit[C->Hvar.get(l, offset)] * fact1;
					betaSum += Beta[offset];
				}
			}
			else
			{
				//Fold any accumulated flat run into this site's factors (when the previous
				//processed site was materialized these are exactly fact2 and fact1).
				const float fact1 = (float)(runD * (C->nt[l] / betaSumNext) + C->t[l] / nstates);
				const float fact2 = (float)(runC * (C->nt[l] / betaSumNext));
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				const __m256 _fact2 = _mm256_set1_ps(fact2);
				//Independent accumulators per reduction to break the add latency chains.
				__m256 _prob0a = _mm256_setzero_ps(), _prob0b = _mm256_setzero_ps();
				__m256 _prob1a = _mm256_setzero_ps(), _prob1b = _mm256_setzero_ps();
				__m256 _suma = _mm256_setzero_ps(), _sumb = _mm256_setzero_ps();

				int k = 0;
				for (; k + 16 <= nstatesMD8; k += 16)
				{
					const __m256i _mc0 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+0)), _vshift_count);
					const __m256i _mc1 = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k+8)), _vshift_count);
					const __m256 _pc0 = _mm256_fmadd_ps(_mm256_load_ps(&Beta[k+0]), _fact2, _fact1);
					const __m256 _pc1 = _mm256_fmadd_ps(_mm256_load_ps(&Beta[k+8]), _fact2, _fact1);
					const __m256 _dp0 = _mm256_mul_ps(_mm256_load_ps(&Alpha[l*modK+k+0]), _pc0);
					const __m256 _dp1 = _mm256_mul_ps(_mm256_load_ps(&Alpha[l*modK+k+8]), _pc1);
					const __m256 _pn0 = _mm256_mul_ps(_pc0, _mm256_blendv_ps(_emit0, _emit1, _mm256_castsi256_ps(_mc0)));
					const __m256 _pn1 = _mm256_mul_ps(_pc1, _mm256_blendv_ps(_emit0, _emit1, _mm256_castsi256_ps(_mc1)));
					_prob0a = _mm256_add_ps(_prob0a, _mm256_mul_ps(_dp0, _mm256_blendv_ps(_one, _zero, _mm256_castsi256_ps(_mc0))));
					_prob1a = _mm256_add_ps(_prob1a, _mm256_mul_ps(_dp0, _mm256_blendv_ps(_zero, _one, _mm256_castsi256_ps(_mc0))));
					_prob0b = _mm256_add_ps(_prob0b, _mm256_mul_ps(_dp1, _mm256_blendv_ps(_one, _zero, _mm256_castsi256_ps(_mc1))));
					_prob1b = _mm256_add_ps(_prob1b, _mm256_mul_ps(_dp1, _mm256_blendv_ps(_zero, _one, _mm256_castsi256_ps(_mc1))));
					_suma = _mm256_add_ps(_suma, _pn0);
					_sumb = _mm256_add_ps(_sumb, _pn1);
					_mm256_store_ps(&Beta[k+0], _pn0);
					_mm256_store_ps(&Beta[k+8], _pn1);
				}
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _mask_curr = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask0 = _mm256_blendv_ps (_one, _zero, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask1 = _mm256_blendv_ps (_zero, _one, _mm256_castsi256_ps(_mask_curr));
					const __m256 _prob_prev = _mm256_load_ps(&Beta[k]);
					const __m256 _alphas = _mm256_load_ps(&Alpha[l*modK+k]);
					const __m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
					const __m256 _dotprod = _mm256_mul_ps(_alphas, _prob_curr);
					const __m256 _prob_next = _mm256_mul_ps(_prob_curr, _emiss);
					_prob0a = _mm256_add_ps(_prob0a, _mm256_mul_ps(_dotprod, _mask0));
					_prob1a = _mm256_add_ps(_prob1a, _mm256_mul_ps(_dotprod, _mask1));
					_suma = _mm256_add_ps(_suma, _prob_next);
					_mm256_store_ps(&Beta[k], _prob_next);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_mm256_add_ps(_prob0a, _prob0b));
					prob_hid[1] = horizontal_add(_mm256_add_ps(_prob1a, _prob1b));
					betaSum = horizontal_add(_mm256_add_ps(_suma, _sumb));
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Beta[offset] = Beta[offset] * fact2 + fact1;
					prob_hid[C->Hvar.get(l, offset)] += Alpha[l*modK+offset] * Beta[offset];
					Beta[offset] *= emit[C->Hvar.get(l, offset)];
					betaSum += Beta[offset];
				}
			}
			runC = 1.0; runD = 0.0; p_anchor = -1;
			prob_hid[0] /= emit[0];
			prob_hid[1] /= emit[1];
			prob_obs[0] = (prob_hid[0]*C->ee_imp + prob_hid[1]*C->ed_imp) * HL[2*abs_site+0];
			prob_obs[1] = (prob_hid[0]*C->ed_imp + prob_hid[1]*C->ee_imp) * HL[2*abs_site+1];

		}
		HP[2*abs_site+0] = prob_obs[0] / (prob_obs[0] + prob_obs[1]);
		HP[2*abs_site+1] = prob_obs[1] / (prob_obs[0] + prob_obs[1]);
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

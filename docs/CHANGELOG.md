---
title: CHANGELOG
layout: default
---

# CHANGELOG

All notable changes to this project are documented in this file.

## Unreleased
	* Added `GLIMPSE2_inspect` tool to print summary statistics (region, genetic-map span, haplotype count, variant-type and allele-frequency breakdowns) for a binary reference panel (`.bin`) produced by `GLIMPSE2_split_reference`, and optionally write its reference haplotypes as VCF/BCF (`--output`).
	* GLIMPSE2_ligate: fixed all-haploid chunks (FPLOIDY=1, e.g. male chrX non-PAR imputed with --samples-file), where consecutive samples were compared as if they were one diploid sample, alleles could be swapped between two samples, and the genotype buffer was read past its end. Haploid chunks are now concatenated without phase alignment.
	* GLIMPSE2_concordance: fixed a crash with --gt-val when the validation file has no FORMAT/DP; with --min-val-dp, missing DP values are counted as missing again rather than as low DP.
	* GLIMPSE2_concordance: --out-r2-per-site now covers every --input line instead of only the last (#179). --out-rej-sites, --out-conc-sites and --out-disc-sites, which kept only the last line, now stop with an error when more than one --input line is given.
	* GLIMPSE2_concordance: fixed the second site being counted twice in runs with --ac-bins/--allele-counts or an imputed file with mixed ploidy, and a hang when the first region has no record in the frequency or imputed file.
	* GLIMPSE2_concordance: fixed every haploid genotype being rejected as NaN when a mixed-ploidy imputed file is subset to haploid samples with --samples.
	* GLIMPSE2_concordance: validation records with missing or too few FORMAT/PL values (e.g. males called haploid against diploid imputed data) are now rejected instead of crashing or being read from the wrong sample or the previous record.
	* GLIMPSE2_concordance: imputed records with missing or too few GT, GP or (with --gt-tar --min-tar-gp) PL values are now rejected instead of crashing or reusing the previous record's values.
	* GLIMPSE2_concordance: calibration bins now hold the GP range they are labelled with ([b/100, (b+1)/100)); previously most GPs were reported one bin low.
	* GLIMPSE2_phase: a --samples-file that cannot be opened, or cannot be read to the end (e.g. a truncated .gz), is now an error; all target samples were silently treated as diploid (a missing .gz path crashed). A samples file that lists no samples now gives a warning, and the number of haploid and diploid target samples it sets is reported.
	* GLIMPSE2_phase: INFO/INFO is now computed for haploid-only output (FPLOIDY=1, e.g. male chrX non-PAR), where it was 1 at every site, and for mixed-ploidy output (FPLOIDY<0), which had no INFO. Diploid output is unchanged.
	* GLIMPSE2_split_reference and GLIMPSE2_chunk: a read error, e.g. from a truncated input file, is now an error; it was taken for the end of the file, giving a partial binary panel or chunks that covered only the part read, with exit 0.
	* GLIMPSE2_phase: FORMAT/GP is now declared Number=G (it holds 2 values for haploid samples), which removes htslib's "GP should be declared as Number=G" warning on every read of the output.
	* GLIMPSE2_phase: the sparse PBWT selection now measures the gap before a common-to-common reselection, and sets the reset point, from the previous high-quality common site; it used an earlier common site once an indel or same-position record had passed. Imputation output changes slightly (within seed-to-seed variation on the tutorial and on HGDP+1KG chr22).
	* GLIMPSE2_phase: fixed a crash on a region where the genetic map is flat (0 cM, e.g. the chr9 centromere) (#304). A region with no common SNPs for the PBWT selection now stops with an error.

## v2.0.0
	* Major release. Introduced speedups and accuracy improvements. Version used for the preprint. (https://doi.org/10.1101/2022.11.28.518213)

## v1.1.1
	* Removed bug in haploid imputation
	* Removed bug reading imputing regions with spanning indels (bug reported by @MySelvan)
	* Removed bug in INFO score calculation
	* Refactory model classes (cast)
	* Added --inputGL option
	* Added --ban-repeated-sample-names option

## v1.1.0
	* Better reporting of missing genotype likelihoods and ploidy
	* Added multi-threaded input for VCF/BCF files
	* Added haploid/diploid/mixed ploidy imputation
	* Added FPLOIDY parameter in output
	* Small improvements in the PBWT selection

## v1.0.1
	* Version used for the paper (https://doi.org/10.1038/s41588-020-00756-0)
	* Change in the state  selection for phasing
	* Bugfix for compatibility with bcftools 1.10.x regarding --main

## v1.0.0
	* First release. Version used for the preprint (https://doi.org/10.1101/2020.04.14.040329)

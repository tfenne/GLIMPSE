---
title: CHANGELOG
layout: default
---

# CHANGELOG

All notable changes to this project are documented in this file.

## Unreleased
	* Added `GLIMPSE2_inspect` tool to print summary statistics (region, genetic-map span, haplotype count, variant-type and allele-frequency breakdowns) for a binary reference panel (`.bin`) produced by `GLIMPSE2_split_reference`, and optionally write its reference haplotypes as VCF/BCF (`--output`).
	* GLIMPSE2_concordance: --out-r2-per-site, --out-rej-sites, --out-conc-sites and --out-disc-sites now cover every --input line instead of only the last (#179); with several --input lines the frequency files must then have the same samples. A per-site output that cannot be opened is now an error instead of a crash. The --out-r2-per-site column `info` is renamed `af_tag_freq`: it holds the frequency used for binning (from --af-tag, or AC/AN if the tag is absent; the minor allele frequency unless --use-alt-af is set), not an INFO score.
	* GLIMPSE2_concordance: fixed a crash with --gt-val when the validation file has no FORMAT/DP; with --min-val-dp, missing DP values are counted as missing again rather than as low DP.
	* GLIMPSE2_concordance: fixed the second site being counted twice in runs with --ac-bins/--allele-counts or an imputed file with mixed ploidy, a hang when the first region has no record in the frequency or imputed file, and every haploid genotype being rejected as NaN when a mixed-ploidy imputed file is subset to haploid samples with --samples.
	* GLIMPSE2_concordance: records with missing or too few FORMAT values, and missing or out-of-range DS values, are now rejected instead of crashing or reusing the previous record's values. Calibration bins now hold the GP range they are labelled with ([b/100, (b+1)/100)); previously most GPs were reported one bin low.

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

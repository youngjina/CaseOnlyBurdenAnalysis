# Case-Only Burden Analysis

## Summary

This repository contains the scripts necessary to conduct a gene-based association test of rare variants using case-only exome data. We make use of recent estimates of mutation rates at positions throughout the genome to calculate the expected mutation rates per gene, taking into account only well covered sites where at least 90% of samples are covered at >10X, and identify genes with significantly higher burden of mutations compared to expectation. To control for the potential confounding effect due to the constraint score for a gene, the expected number of variants for a gene is adjusted for that gene’s constraint score.

## Outline

### 1. Select_covered_positions_by_percent.R
Script for selecting positions that are covered at 90% of samples
Usage: Rscript Select_covered_positions_by_percent.R

### 2. position_variantID.py
Script for converting format (Input: gene,chr,pos; Output: chr-pos"\t"gene)
Usage: python position_variantID.py 

### 3. Split_covered_positions.py
Script for splitting positions that are covered at 90% of samples in the coverage files according to the chromosomes
Usage: python Split_covered_positions.py

### 4. Generate_AllPossibleVariants_WellCoveredRegion.py
Script for extracting rare well-covered mutations
Usage: python Generate_AllPossibleVariants_WellCoveredRegion.py

### 5. Extract_MutationRate_AllPossibleVariants_WellCoveredRegion.py
Script for annotating mutation rate of variants
Usage: python Extract_MutationRate_AllPossibleVariants_WellCoveredRegion.py

# CaseOnlyBurdenAnalysis
## Adjust for Read Depth
- Bases with at least 90% of samples covered at 10X for the cohort are retained. 

- `/nfs/goldstein/software/sh/atav.sh --site-coverage-summary --sample SAMPLE --gene-boundary /nfs/goldstein/software/atav_home/data/ccds/addjusted.CCDS.genes.index.r20.hg19.r15names.txt --min-coverage 10 --out COVERAGE_OUTPUT`
- Purpose: Generate coverage per position
- Input: sample fam file
- Output: coverage per position
- Caveat: For each cohort, input file and output file names are different.

- `Rscript Select_covered_positions_by_percent.R COVERAGE_OUTPUT Coverage_Sum_Percent_"$COHORT"_min-coverage_10_site.summary.csv Covered_90_percent_Positions_"$COHORT"_min-coverage_10_site.summary.csv $NUM_SAMPLES 90 "$PHENOTYPES"`
- Purpose: Select positions that are covered at 90% of samples
- Input: coverage per position for only case
- Output 1: positions that are covered at 90% of samples
- Output 2: all positions including columns coverage.sum and coverage.percent
- Caveat: number of samples in the parameters

- `cut -d',' -f2,3 --output-delimiter='-' Covered_90_percent_Positions_"$COHORT"_min-coverage_10_site.summary.csv > Covered_90_percent_Positions_"$COHORT"_min-coverage_10_chr_pos.txt`
- Purpose: Extract chr and pos info of positions that are covered at 90% of samples
- Input: positions that are covered at 90% of samples
- Output: chr and pos info of positions that are covered at 90% of samples

- `cut -d',' -f1-3 Covered_90_percent_Positions_"$COHORT"_min-coverage_10_site.summary.csv > Covered_90_percent_Positions_Gene_Chr_Pos_"$COHORT"_min-coverage_10_site.summary.csv`
- Purpose: Extract gene, chr and pos info of positions that are covered at 90% of samples
- Input: positions that are covered at 90% of samples
- Output: gene, chr and pos info of positions that are covered at 90% of samples

- `python /media/Dazs/youngji/Scripts/position_variantID.py /media/Dazs/youngji/Coverage_Analysis/Coverage_Summary_"$COHORT"_min-coverage_10/ Covered_90_percent_Positions_Gene_Chr_Pos_"$COHORT"_min-coverage_10_site.summary.csv Covered_90_percent_Positions_Gene_variantID_"$COHORT"_min-coverage_10_site.summary.txt`
- Purpose: Make format (chr-pos"\t"gene)
- Input: gene,chr,pos
- Output: chr-pos"\t"gene

- `python /media/Dazs/youngji/Scripts/modified_split_covered_positions.py /media/Dazs/youngji/Coverage_Analysis/Coverage_Summary_"$COHORT"_min-coverage_10/ Covered_90_percent_Positions_"$COHORT"_min-coverage_10_chr_pos.txt "$COHORT"`
- Purpose: Split coverage files according to the chromosomes
- Input: In one file, positions that are covered at 90% of samples
- Output: for each chromosome file, positions that are covered at 90% of samples


## All possible variants to calculate expected variants

## Make mutation rate
- Generate_AllPossibleVariants_WellCoveredRegion.sh
- Extract_MutationRate_AllPossibleVariants_WellCoveredRegion.sh

## Run the case-only burden analysis
- CaseOnlyBurdenAnalysis_LOF.R
- CaseOnlyBurdenAnalysis_Missense.R
- CaseOnlyBurdenAnalysis_LOF_Missense.R

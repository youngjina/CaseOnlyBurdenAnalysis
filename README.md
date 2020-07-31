# CaseOnlyBurdenAnalysis
## Adjust for Read Depth
- Bases with at least 90% of samples covered at 10X for each cohort are retained.
  
- Generate coverage per position
- `/nfs/goldstein/software/sh/atav.sh --site-coverage-summary --sample SAMPLE --gene-boundary /nfs/goldstein/software/atav_home/data/ccds/addjusted.CCDS.genes.index.r20.hg19.r15names.txt --min-coverage 10 --out COVERAGE_OUTPUT`

- Select positions that are covered at 90% of samples
- `Rscript Select_covered_positions_by_percent.R COVERAGE_OUTPUT Coverage_Sum_Percent_"$COHORT"_min-coverage_10_site.summary.csv Covered_90_percent_Positions_"$COHORT"_min-coverage_10_site.summary.csv $NUM_SAMPLES 90 "$PHENOTYPES"`
  
- Extract chr and pos info of positions that are covered at 90% of samples 
- `cut -d',' -f2,3 --output-delimiter='-' Covered_90_percent_Positions_"$COHORT"_min-coverage_10_site.summary.csv > Covered_90_percent_Positions_"$COHORT"_min-coverage_10_chr_pos.txt`
    
- Extract gene, chr and pos info of positions that are covered at 90% of samples
- `cut -d',' -f1-3 Covered_90_percent_Positions_"$COHORT"_min-coverage_10_site.summary.csv > Covered_90_percent_Positions_Gene_Chr_Pos_"$COHORT"_min-coverage_10_site.summary.csv`
    
- Make format (chr-pos"\t"gene)
- `python position_variantID.py /Coverage_Summary_"$COHORT"_min-coverage_10/ Covered_90_percent_Positions_Gene_Chr_Pos_"$COHORT"_min-coverage_10_site.summary.csv Covered_90_percent_Positions_Gene_variantID_"$COHORT"_min-coverage_10_site.summary.txt`
  
- Split coverage files according to the chromosomes
- `python Split_covered_positions.py /Coverage_Summary_"$COHORT"_min-coverage_10/ Covered_90_percent_Positions_"$COHORT"_min-coverage_10_chr_pos.txt "$COHORT"`

## Extract observed variants 
### LOF
`atav.sh --list-var-geno --sample Fam_CKD_n553.txt --effect "$LOF_EFFECTS" --min-coverage 10 --exclude-artifacts --exclude-evs-qc-failed --ccds-only --include-qc-missing --filter PASS,LIKELY,INTERMEDIATE --qd 2 --qual 50 --mq 40 --gq 20 --rprs -3 --mqrs -10 --snv-fs 60 --indel-fs 200 --snv-sor 3 --indel-sor 10 --het-percent-alt-read 0.3-1 --max-qc-fail-sample 0 --gnomad-exome-rf-tp-probability-snv 0.1 --gnomad-exome-rf-tp-probability-indel 0.2 --gnomad-exome-pop global --gnomad-exome-af 0 --loo-af 5e-04 --include-rvis --out CKD_n553_LOF_gnomad_0_looaf_5e-04`

`LOF_EFFECTS="HIGH:exon_loss_variant,HIGH:frameshift_variant,HIGH:rare_amino_acid_vari
ant,HIGH:stop_gained,HIGH:stop_lost,HIGH:start_lost,HIGH:gene_fusion,HIGH:bidirection
al_gene_fusion,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant"`


### Nonsynonymous
`atav.sh --list-var-geno --sample Fam_CKD_n553.txt --effect "$NONSYNONYMOUS_EFFECTS" --min-coverage 10 --exclude-artifacts --exclude-evs-qc-failed --ccds-only --include-qc-missing --filter PASS,LIKELY,INTERMEDIATE --qd 2 --qual 50 --mq 40 --gq 20 --rprs -3 --mqrs -10 --snv-fs 60 --indel-fs 200 --snv-sor 3 --indel-sor 10 --het-percent-alt-read 0.3-1 --max-qc-fail-sample 0 --gnomad-exome-rf-tp-probability-snv 0.1 --gnomad-exome-rf-tp-probability-indel 0.2 --gnomad-exome-pop global --gnomad-exome-af 0 --loo-af 5e-04 --include-rvis --out CKD_n553_Nonsynonymous_gnomad_0_looaf_5e-04`

`NONSYNONYMOUS_EFFECTS="MODERATE:3_prime_UTR_truncation+exon_loss_variant,MODERATE:5_prime_UTR_truncation+exon_loss_variant,MODERATE:coding_sequence_variant,MODERATE:disruptive_inframe_deletion,MODERATE:disruptive_inframe_insertion,MODERATE:conservative_inframe_deletion,MODERATE:conservative_inframe_insertion,MODERATE:missense_variant+splice_region_variant,MODERATE:missense_variant,LOW:5_prime_UTR_premature_start_codon_gain_variant,LOW:initiator_codon_variant,LOW:initiator_codon_variant+non_canonical_start_codon"`


### LOF + Nonsynonymous
`atav.sh --list-var-geno --sample Fam_CKD_n553.txt --effect "$FUNCTIONAL_EFFECTS" --min-coverage 10 --exclude-artifacts --exclude-evs-qc-failed --ccds-only --include-qc-missing --filter PASS,LIKELY,INTERMEDIATE --qd 2 --qual 50 --mq 40 --gq 20 --rprs -3 --mqrs -10 --snv-fs 60 --indel-fs 200 --snv-sor 3 --indel-sor 10 --het-percent-alt-read 0.3-1 --max-qc-fail-sample 0 --gnomad-exome-rf-tp-probability-snv 0.1 --gnomad-exome-rf-tp-probability-indel 0.2 --gnomad-exome-pop global --gnomad-exome-af 0 --loo-af 5e-04 --include-rvis --out CKD_n553_LOF_Nonsynonymous_gnomad_looaf_5e-04`

`FUNCTIONAL_EFFECTS="HIGH:exon_loss_variant,HIGH:frameshift_variant,HIGH:rare_amino_ac
id_variant,HIGH:stop_gained,HIGH:start_lost,HIGH:stop_lost,HIGH:splice_acceptor_varia
nt,HIGH:splice_donor_variant,HIGH:gene_fusion,HIGH:bidirectional_gene_fusion,MODERATE
:3_prime_UTR_truncation+exon_loss_variant,MODERATE:5_prime_UTR_truncation+exon_loss_v
ariant,MODERATE:coding_sequence_variant,MODERATE:disruptive_inframe_deletion,MODERATE
:disruptive_inframe_insertion,MODERATE:conservative_inframe_deletion,MODERATE:conserv
ative_inframe_insertion,MODERATE:missense_variant+splice_region_variant,MODERATE:miss
ense_variant"`


## All possible variants to calculate expected variants
### LOF
- `python cadd_anno_LOF.py "$JOBDIR" whole_genome_SNVs_inclAnno.tsv.gz LOF_cadd_anno.txt`
- `cut -f1 LOF_cadd_anno.txt > LOF_cadd_anno_variant.txt`
- Split CADD LOF files according to the chromosomes
- `python split_cadd_variant.py /Chr_LOF_CADD/ LOF_cadd_anno_variant.txt LOF`

### Nonsynonymous
- `python cadd_anno_nonsynonymous.py "$JOBDIR" whole_genome_SNVs_inclAnno.tsv.gz Nonsynonymous_cadd_anno.txt`
- `cut -f1 Nonsynonymous_cadd_anno.txt > Nonsynonymous_cadd_anno_variant.txt`
- Split CADD Nonsynonymous files according to the chromosomes
- `python split_cadd_variant.py /Chr_Nonsynonymous_CADD/ Nonsynonymous_cadd_anno_variant.txt Nonsynonymous`





## Make mutation rate
- Generate_AllPossibleVariants_WellCoveredRegion.sh
- Extract_MutationRate_AllPossibleVariants_WellCoveredRegion.sh

## Run the case-only burden analysis
- CaseOnlyBurdenAnalysis_LOF.R
- CaseOnlyBurdenAnalysis_Missense.R
- CaseOnlyBurdenAnalysis_LOF_Missense.R

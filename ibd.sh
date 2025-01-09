#!/bin/bash
# rename columns by removing the '_' 
sed -e '/^#CHROM/ s/_//g' liftover37_merged_panGen2024_sample_filtered.norm.recode.vcf > renamed_columns.vcf
sed -e '/^#CHROM/ s/_//g' merged_panGen2019_samples_filtered.norm.vcf > renamed_columns.vcf

# ibd calculation
plink --vcf renamed_columns.vcf --genome
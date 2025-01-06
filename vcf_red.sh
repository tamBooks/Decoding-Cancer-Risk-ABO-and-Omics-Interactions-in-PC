# Convert the normalized VCF to PLINK binary format
plink2 --vcf merged_panGen2024.vcf.gz --make-bed --out merged_panGen2024 --chr 1-22

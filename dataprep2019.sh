# Data Preperation for imputation
# Run perl script that is set by the Imputatipn server
perl /local/tsaid/trabajo/dataprep/HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b /local/tsaid/trabajo/vcf/merged2024/merged_panGen2024.bim -f /local/tsaid/trabajo/vcf/merged2024/merged_panGen2024.afreq -r /local/tsaid/trabajo/dataprep/HRC-1000G-check-bim-v4.3.0/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h sh Run-plink.sh

# Edit and run the plink script made by the perl script
bash Run-plink.sh

# Compress the chromosome files for imputation
bgzip -c new_ref_merge-updated-chr1.vcf > new_ref_merge-updated-chr1.vcf.gz

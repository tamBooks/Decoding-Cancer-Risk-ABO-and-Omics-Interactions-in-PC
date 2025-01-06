

#vcftools --gzvcf merged_panGen2024.norm.vcf.gz --indv 207769820036_R01C01 --indv 208025880084_R01C01 --indv 208025880086_R02C02 --indv 208025880090_R03C01 --indv 208025880092_R04C02 --indv 208025880109_R05C01 --indv 208025880110_R06C02 --indv 208025880114_R07C01 --indv 208025880116_R08C01 --indv 208026130054_R09C02 --indv 208026130069_R10C01 --recode --recode-INFO-all --out /local/tsaid/trabajo/samples_10_only

vcftools --gzvcf merged_panGen2024.norm.vcf.gz --remove remove_samples.txt --recode --recode-INFO-all --out merged_panGen2024_sample_filtered.norm

#vcftools --gzvcf merged_panGen2024.norm.vcf.gz --keep samples.txt --recode --recode-INFO-all --out merged_panGen2024_sample_filtered.norm

plink2 --vcf merged_panGen2024_sample_filtered.norm.recode.vcf.gz --geno 0.2 --mind 0.2 --recode ped --make-bed --out merged_panGen2024_sample_filtered.norm.recode --chr 1-22
plink2 --bfile merged_panGen2024_sample_filtered.norm.recode --recode vcf --out merged_panGen2024_sample_filtered.norm.recode
plink2 --bfile merged_panGen2024_sample_filtered.norm.recode --freq --out merged_panGen2024_sample_filtered.norm.recode
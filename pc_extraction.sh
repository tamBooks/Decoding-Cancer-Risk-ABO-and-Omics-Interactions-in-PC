#!/bin/bash

for chr in {1..22}; do
    awk '{print $5 "\t" $6}' /local/tsaid/sharing/PCA_SNPlist_dosages_oncoarray_pangen.txt > /local/tsaid/sharing/pos.txt
    zfgrep -f  pos.txt /local/tsaid/trabajo/imputation/2024_Michigan_23_Oct/chr${chr}.dose.vcf.gz >> /local/tsaid/trabajo/imputation/2024_Michigan_23_Oct/GSA2024_dosagesPCA.txt
    zfgrep -f  pos.txt /local/tsaid/trabajo/imputation/2019_Michigan_28_Oct/chr${chr}.dose.vcf.gz >> /local/tsaid/trabajo/imputation/2019_Michigan_28_Oct/GSA2019_dosagesPCA.txt
done 

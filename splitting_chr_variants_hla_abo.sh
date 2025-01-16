# HLA extraction
zcat chr6.dose.vcf.gz | grep "^#"  > HLA_region.vcf
zcat chr6.dose.vcf.gz | awk '{FS="\t"} {if ($2>=27970031 && $2<=33965553) {print $0}}' >> HLA_region.vcf

zcat chr6.info.gz | grep "^#"  > HLA_region.info
zcat chr6.info.gz | awk '{FS="\t"} {if ($2>=27970031 && $2<=33965553) {print $0}}' >> HLA_region.info

# ABO extraction
zcat chr9.dose.vcf.gz | grep "^#"  > abo_region.vcf
zcat chr9.dose.vcf.gz | awk '{FS="\t"} {if ($2 >= 136115788 && $2 <= 136160617) {print $0}}' >> abo_region.vcf

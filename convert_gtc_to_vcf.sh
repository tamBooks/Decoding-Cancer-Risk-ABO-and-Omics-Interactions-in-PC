#!/bin/bash

# Set the base paths, manifest files, and reference genome
base_path_to_gtc_folder="/local/tsaid/trabajo/gtc/"
bpm_manifest_file="/local/comun/microbiomePanGenGenotypes/rawDataCEGEN2024/GSAMD-24v3-0-EA_20034606_A2_GRCh38.bpm"
csv_manifest_file="/local/comun/microbiomePanGenGenotypes/GSAMD-24v3-0-EA_20034606_A2_GRCh38.csv"
egt_cluster_file="/local/comun/microbiomePanGenGenotypes/rawDataCEGEN2024/GSA3MD_fase3_1000G_monomorficos_nonZeroed_chrY_rep.egt"
#ref="/local/comun/dbs/GRCh37/human_g1k_v37.fasta"
#output_dir="/local/tsaid/trabajo/new_ref/"
ref="/local/comun/dbs/GRCh38_NCBI_p13/GRCh38.p13_assembly.fa"
output_dir="/local/tsaid/trabajo/"


# Function to check if file exists
file_exists() {
  local file=$1
  if [ -f "$file" ]; then
    echo "File $file already exists, skipping this step."
    return 0
  else
    return 1
  fi
}

# Convert Illumina GTC files to VCF
convert_gtc_to_vcf() {
  local gtc_folder=$1
  local vcf_folder=$2
  output_file="$vcf_folder/output_file.vcf"
  if ! file_exists "$output_file"; then
    bcftools +gtc2vcf \
      --no-version -Ou \
      --do-not-check-bpm \
      --csv $csv_manifest_file \
      --egt $egt_cluster_file \
      --gtcs $gtc_folder \
      --fasta-ref $ref \
      --extra $vcf_folder/output_file.tsv | \
    bcftools sort -Ou -T ./bcftools. | \
    bcftools norm --no-version -o $output_file -Ov -c x -f $ref 
    bgzip -c $output_file > $output_file.gz 
    bcftools index $output_file.gz
    tabix -p vcf $output_file.gz
    if [ -f "$output_file" ]; then
      echo "Removing uncompressed VCF file..."
      rm $output_file
    else
      echo "Uncompressed VCF file not found, nothing to remove."
    fi
  else
    echo "VCF file already exists, skipping conversion."
  fi
}

#copy() {
#    output_file="${vcf_folder}/output_file.vcf"
 #   cp $output_file ${vcf_folder}copy.vcf 
   #bgzip -c ${vcf_folder}copy.vcf > ${vcf_folder}toindex.vcf.gz 
    #bcftools index ${vcf_folder}toindex.vcf.gz
#}

#The final VCF might contain duplicates. If this is an issue bcftools norm -d exact can be used to remove such variants. 

# Main execution
for gtc_dir in $base_path_to_gtc_folder/*/; do
  output_gtc=`echo $gtc_dir | cut -f7 -d '/'`
  vcf_folder="${output_dir}vcf_csv/${output_gtc}_vcf"  # Define output VCF folder based on GTC folder
  mkdir -p $vcf_folder  # Ensure VCF output folder exists
  echo "Converting GTC files to VCF in directory: $gtc_dir"
  convert_gtc_to_vcf $gtc_dir $vcf_folder 
done

#for vcf_folder in `ls /local/tsaid/trabajo/vcf | head -1`; do
 #   copy $vcf_folder  
#done


#!/bin/bash

# Set the base paths and manifest files
base_path_to_idat_folder="/local/comun/microbiomePanGenGenotypes/rawDataCEGEN2019/GSA_IDATs/"
bpm_manifest_file="/local/comun/microbiomePanGenGenotypes/rawDataCEGEN2019/GSAMD-24v2-0_20024620_A1.bpm"
egt_cluster_file="/local/comun/microbiomePanGenGenotypes/rawDataCEGEN2019/Cluster_GSA_v2_CEGEN_revisado_centenarios_1000GENOMAS_duplicados.egt"
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

# Convert Illumina IDAT files to GTC files
convert_idat_to_gtc() {
  local idat_folder=$1
  local gtc_folder=$2
  output_file="$gtc_folder/output_file.gtc"  # Specify the expected output file
  if ! file_exists "$output_file"; then
    bcftools +idat2gtc \
      --bpm $bpm_manifest_file \
      --egt $egt_cluster_file \
      --idats $idat_folder \
      --output $gtc_folder
  fi
}

# Main execution
for idat_dir in $base_path_to_idat_folder/*/; do
  output_idat=`echo $idat_dir | cut -f8 -d '/'`
  gtc_folder="${output_dir}gtc19/${output_idat}_gtc"  # Define GTC output folder based on IDAT folder
  mkdir -p $gtc_folder  # Ensure GTC output folder exists
  echo "Converting IDAT files to GTC in directory: $idat_dir"
  convert_idat_to_gtc $idat_dir $gtc_folder
done

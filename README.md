# Decoding-Cancer-Risk-ABO-and-Omics-Interactions-in-PC
## PanGenEU Data Processing Pipeline

This repository contains scripts and workflows to process genetic data from raw IDAT files to VCF files, including preparation for imputation using the Michigan Imputation Server. Below is a detailed description of each script and its purpose.

---

## Table of Contents
1. [Overview](#overview)  
2. [Scripts and Workflow](#scripts-and-workflow)  
    - [IDAT to GTC Conversion](#idat-to-gtc-conversion)
    - [GTC to VCF Conversion](#gtc-to-vcf-conversion)
    - [VCF Merging and Liftover](#vcf-merging-and-liftover)
    - [Sample Comparison and Filtering](#sample-comparison-and-filtering)
    - [VCF Tools and Pre-Imputation Filtering](#vcf-tools-and-pre-imputation-filtering)
    - [Imputation Data Preparation](#imputation-data-preparation)  
3. [Requirements](#requirements)  
4. [Usage](#usage)  

---

## Overview
This pipeline processes genetic data using various tools and scripts for data preparation, sample filtering, and format conversion. It supports two references:  
- **Reference 38**: PanGenEU 2024  
- **Reference 37**: PanGenEU 2019  

The scripts are organized into key steps, which are explained in detail below.

---

## Scripts and Workflow

### 1. IDAT to GTC Conversion
These scripts convert raw IDAT data files to GTC files.  
- **`convert_idat_to_gtc.sh`**: Converts IDAT files to GTC files using **Reference 38** (PanGenEU 2024).  
- **`convert_idat_to_gtc19.sh`**: Converts IDAT files to GTC files using **Reference 37** (PanGenEU 2019).  
- **`convert_idat_to_gtc19_skipsamples.sh`**: A modified version of the conversion script that skips specific samples listed in an input file.  

---

### 2. GTC to VCF Conversion
These scripts convert the GTC files into VCF files.  
- **`convert_gtc_to_vcf.sh`**: Converts GTC files to VCF using **Reference 38** (PanGenEU 2024).  
- **`convert_gtc_to_vcf19.sh`**: Converts GTC files to VCF using **Reference 37** (PanGenEU 2019).  

---

### 3. VCF Merging and Liftover
#### Merging:
- **`merge_vcfs.sh`**: Merges and compresses multiple VCF files into a single, combined VCF file.  

#### Liftover (Reference 38 → Reference 37):
- **`vcf_red.sh`**: Reduces the size of the merged VCF file to make it manageable for liftover.  
- **`vcf_liftover.py`**: Performs the liftover of the VCF file from Reference 38 to Reference 37, ensuring compatibility with downstream processes like the Michigan Imputation Server.  

---

### 4. Sample Comparison and Filtering
- **`extract_pairs_and_compare.sh`**: Identifies duplicate samples and calculates misalignment, agreement, and disagreement rates. This helps decide which samples to keep or remove for downstream analysis.  

---

### 5. VCF Tools and Pre-Imputation Filtering
- **`vcftools_script.sh`**: Filters VCF files using the following criteria:
  - Excludes samples and SNPs with high missing rates.
  - Generates new VCF files with a call rate of at least 80%.
  - Produces `.bim` and allele frequency files required for imputation.  

---

### 6. Imputation Data Preparation
Scripts for preparing the filtered VCF files for submission to the Michigan Imputation Server.  
- **`dataprep2024.sh`**: Processes data for **PanGenEU 2024** (Reference 38).  
- **`dataprep2019.sh`**: Processes data for **PanGenEU 2019** (Reference 37).  
These scripts:
  1. Use the provided Perl script from the server to prepare the data.
  2. Run the generated script for final processing.

---

## Requirements
- **Tools**:  
  - `bcftools`
  - `bcftools` + `idat2gtc`
  - `bcftools` + `gtc2vcf`
  - `vcftools`  
  - `plink`  
  - Python with required libraries for `vcf_liftover.py`  
  - Michigan Imputation Server's Perl scripts and files 

- **Input Files**:  
  - Raw IDAT files  
  - Sample lists (for skipping samples or comparisons)  

---

## Usage
1. **Convert IDAT files to GTC**:
   ```bash
   ./convert_idat_to_gtc.sh
   ./convert_idat_to_gtc19.sh
   ```

2. **Convert GTC to VCF**:
   ```bash
   ./convert_gtc_to_vcf.sh
   ./convert_gtc_to_vcf19.sh
   ```

3. **Merge and Liftover**:
   ```bash
   ./merge_vcfs.sh
   ./vcf_red.sh
   python3 vcf_liftover.py
   ```

4. **Filter and Prepare VCF for Imputation**:
   ```bash
   ./vcftools_script.sh
   ./dataprep2024.sh
   ./dataprep2019.sh
   ```

5. **Duplicate Sample Analysis**:
   ```bash
   ./extract_pairs_and_compare.sh
   ```

---

## Acknowledgments
This pipeline was developed for processing genetic data for the PanGenEU project. Special thanks to the Michigan Imputation Server for providing the tools and scripts used in data preparation.

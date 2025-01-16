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
    - [HLA and ABO Analysis](#hla-and-abo-analysis)  
    - [VCF Tools and Pre-Imputation Filtering](#vcf-tools-and-pre-imputation-filtering)  
    - [Imputation Data Preparation](#imputation-data-preparation)  
    - [Post-Imputation Analysis](#post-imputation-analysis)  
    - [Association and Meta-Analysis](#association-and-meta-analysis)  
3. [Requirements](#requirements)  
4. [Usage](#usage)
5. [References](#references)
   
---

## Overview  
This pipeline processes genetic data using various tools and scripts for data preparation, sample filtering, and format conversion. It supports two references:  
- **Reference 38**: PanGenEU 2024  
- **Reference 37**: PanGenEU 2019  

The scripts are organized into key steps, which are explained in detail below.  

## Pipeline Overview

Below is a visual representation of the pipeline from [IDAT to GTC Conversion](#idat-to-gtc-conversion) until [Imputation Data Preparation](#imputation-data-preparation):

![Pipeline Overview](![workflow](https://github.com/user-attachments/assets/f7c0992a-7ed0-4704-81e4-1b83c0f78f16))
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

### 5. HLA and ABO Analysis  
#### HLA and ABO Extraction:  
- **`splitting_chr_variants_hla_abo.sh`**: Extracts the **HLA region** from Chromosome 6 and the **ABO region** from Chromosome 9. The HLA region is reimputed using the Michigan Imputation Server.  

#### Quality Check and Filtering:  
- **`quality_check_HLA_b4_imp.py`**: Performs a quality check on the extracted HLA region, filtering out variants with an R-squared value below 0.7.  
- **`imputed_hla_filteration.py`**: After HLA imputation, filters out variants with an R-squared value below 0.3.  

#### ABO Inference:  
- **`ABO_predict_final.py`**: Infers blood group types based on the ABO region of the sample.  

---

### 6. VCF Tools and Pre-Imputation Filtering  
- **`vcftools_script.sh`**: Filters VCF files using the following criteria:  
  - Excludes samples and SNPs with high missing rates.  
  - Generates new VCF files with a call rate of at least 80%.  
  - Produces `.bim` and allele frequency files required for imputation.  

---

### 7. Imputation Data Preparation  
Scripts for preparing the filtered VCF files for submission to the Michigan Imputation Server.  
- **`dataprep2024.sh`**: Processes data for **PanGenEU 2024** (Reference 37).  
- **`dataprep2019.sh`**: Processes data for **PanGenEU 2019** (Reference 37).  
These scripts:  
  1. Use the provided Perl script from the server to prepare the data.  
  2. Run the generated script for final processing.  

---

### 8. Post-Imputation Analysis  
#### Principal Component Analysis:  
- **`pc_extraction.sh`**: Extracts dosage data for ethnicity-related SNPs from imputed VCF files.  
- **`psa_cal.R`**: Calculates the Principal Component Analysis (PCA), extracting the first five components and saving them separately.
- **`ibd.sh`**`: Identifies samples with a genetic relatedness coefficient higher than 0.9.  

#### Dataframe Merging:  
- **`abo_txt_prep.py `**: Modification and combination of the ABO blood group datasets.
- **`datatable_merge.R`**: Combines the PCA components and other variables into a single, large dataframe for downstream analysis.

---

### 9. Association and Meta-Analysis  
- **`LR_and_meta_analysis.R`**:  
  This script is used for:  
  1. Association analysis between various genetic and environmental variables.  
  2. Meta-analysis to combine results across different datasets.  
  3. Interaction analysis to identify relationships between variables like blood group, HLA alleles, and phenotype.  

---

## Requirements  
- **Tools**:  
  - `bcftools`
  - `bcftools` + `idat2gtc`
  - `bcftools` + `gtc2vcf`
  - `vcftools`  
  - `plink`
  - `plink2`  
  - Python with required libraries for `vcf_liftover.py`, `quality_check_HLA_b4_imp.py`, and `imputed_hla_filteration.py`  
  - Michigan Imputation Server's Perl scripts  
  - R with required libraries for PCA, meta-analysis, and data merging and visualization scripts  

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

4. **Extract HLA and ABO Regions**:  
   ```bash  
   ./splitting_chr_variants_hla_abo.sh  
   python3 quality_check_HLA_b4_imp.py  
   python3 imputed_hla_filteration.py  
   ```  

5. **Filter and Prepare VCF for Imputation**:  
   ```bash  
   ./vcftools_script.sh  
   ./dataprep2024.sh  
   ./dataprep2019.sh  
   ```  

6. **Run Principal Component Analysis**:  
   ```bash  
   ./pc_extraction.sh  
   Rscript psa_cal.R
   ./ibd.sh
   Rscript psa_cal.R
   ```  

7. **Combine Components and Variables**:  
   ```bash  
   Rscript datatable_merge.R
   python3 abo_txt_prep.py  
   ```  

8. **Perform Association and Meta-Analysis**:  
   ```bash  
   Rscript LR_and_meta_analysis.R  
   ```  

---

## Acknowledgments   
This pipeline was developed for processing genetic data of PanGenEU for the Master Thesis project. Special thanks to the Michigan Imputation Server for providing the tools and scripts used in data preparation.

---
## References  
Some of the codes used in this pipeline are adapted or based on tools from external repositories:  

- **IDAT to GTC Conversion**: The `convert_idat_to_gtc.sh` and `convert_idat_to_gtc19.sh` scripts are based on tools from [gtc2vcf GitHub Repository](https://github.com/freeseek/gtc2vcf.git).  
- **GTC to VCF Conversion**: The `convert_gtc_to_vcf.sh` and `convert_gtc_to_vcf19.sh` scripts are based on tools from [gtc2vcf GitHub Repository](https://github.com/freeseek/gtc2vcf.git).  


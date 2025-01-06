import pandas as pd
import gzip
import io
import re

def read_vcf(filepath):
    with gzip.open(filepath, 'rt') as file:
        lines = [line for line in file if not line.startswith('##')]
    header_line = lines[0].strip()
    data = ''.join(lines[1:])
    df_vcf = pd.read_csv(io.StringIO(data), sep='\t', comment='#', names=header_line.lstrip('#').split('\t'))
    return df_vcf

def write_vcf(df, original_vcf_path, output_vcf_path):
    # Open the original VCF file for reading the header lines
    with gzip.open(original_vcf_path, 'rt') as f:
        header_lines = [line.strip() for line in f if line.startswith('##')]
    
    # Write the header lines back to the new VCF file
    with open(output_vcf_path, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        
        # Write the column names (header line)
        f.write('#' + '\t'.join(df.columns) + '\n')
        
        # Write the data
        for index, row in df.iterrows():
            f.write('\t'.join(map(str, row.values)) + '\n')

# Path to your VCF file
#vcf_path = '/local/tsaid/sharing/HLA/HLA_4digit1_multiethnic_2024/chr6.dose.vcf.gz'
vcf_path = '/local/tsaid/sharing/HLA/HLA_4digit1_multiethnic_2019/chr6.dose.vcf.gz'
#info_path = '/local/tsaid/sharing/HLA/HLA_4digit1_multiethnic_2024/chr6.info.gz'
info_path = '/local/tsaid/sharing/HLA/HLA_4digit1_multiethnic_2019/chr6.info.gz'

# Read the VCF file
df_vcf = read_vcf(vcf_path)

# Read and filter the .info.gz file
with gzip.open(info_path, 'rt') as f:
    df_info = pd.read_csv(f, sep='\s+',skiprows=12)

# Convert the Rsq column to numeric, forcing errors to NaN
#df_info['Rsq'] = pd.to_numeric(df_info['Rsq'], errors='coerce')

# Filter the DataFrame where Rsq >= 0.3 and SNP contains "HLA"
#filtered_df = df_info[(df_info['Rsq'] >= 0.3) & (df_info['SNP'].str.contains("HLA", na=False))]

# Extract the filtered SNPs
#filtered_snps = set(filtered_df['SNP'])

# Filter the VCF DataFrame using the SNPs from the filtered info DataFrame
#filtered_vcf_df = df_vcf[df_vcf['ID'].isin(filtered_snps)]

# Function to extract R2 from the INFO column
def extract_r2(info_str):
    # Use a regular expression to search for 'R2' in the INFO field
    match = re.search(r'R2=([0-9.]+)', info_str)
    if match:
        return float(match.group(1))  # Convert the extracted value to float
    return None  # Return None if R2 is not found

# Apply the function to create a new 'R2' column
df_info['R2'] = df_info['INFO'].apply(extract_r2)

# Filter the DataFrame where R2 >= 0.3 and SNP contains "HLA"
filtered_df = df_info[(df_info['R2'] >= 0.3) & (df_info['ID'].str.contains("HLA", na=False))]

# Extract the filtered SNPs
filtered_snps = set(filtered_df['ID'])

# Filter the VCF DataFrame using the SNPs from the filtered info DataFrame
filtered_vcf_df = df_vcf[df_vcf['ID'].isin(filtered_snps)]

# Save the filtered VCF DataFrame to a new VCF file
#output_file = '/local/tsaid/sharing/HLA/HLA_4digit1_multiethnic_2024/filtered_hla.vcf'
output_file = '/local/tsaid/sharing/HLA/HLA_4digit1_multiethnic_2019/filtered_hla.vcf'
write_vcf(filtered_vcf_df, vcf_path, output_file)

'''
import pysam

# Read the filtered SNPs into a set for quick lookup
filtered_snps = set(filtered_df2['SNP'])

vcf_path= ('/local/tsaid/trabajo/imputation/hla_chrom6/chr6.dose.vcf.gz')

pysam.tabix_index(vcf_path, preset='vcf', force=True)

# Open the input VCF file
vcf_in = pysam.VariantFile(vcf_path)

# Open the output VCF file
vcf_out = pysam.VariantFile("filtered_hla.vcf.gz", "w", header=vcf_in.header)

# Iterate through the input VCF records and filter for HLA entries
for record in vcf_in:
    if record.id in filtered_snps:
        vcf_out.write(record)

# Close the files
vcf_in.close()
vcf_out.close()
'''
import pandas as pd
import gzip

def read_vcf(filepath):
    with gzip.open(filepath, 'rt') as file:
        # Skip header lines starting with '##'
        lines = [line.strip() for line in file if not line.startswith('##')]
    # The first non-header line contains the column names
    header_line = lines[0]
    # Extract column names from the header line
    columns = header_line.lstrip('#').split('\t')
    # Read the data into a DataFrame
    df = pd.read_csv(filepath, comment='#', sep='\t', names=columns, compression='gzip', low_memory=False)
    return df

def parse_info_column(info_string, desired_keys):
    # Split the INFO column into key-value pairs
    info_pairs = info_string.split(';')
    info_dict = {}
    for pair in info_pairs:
        if '=' in pair:
            key, value = pair.split('=')
            if key in desired_keys:
                info_dict[key] = float(value)
    return info_dict


def read_and_filter_info(filepath, rsq_threshold=0.7):
    # Open the gzip file and read it directly into a DataFrame
    with gzip.open(filepath, 'rt') as file:
        df = pd.read_csv(file, sep='\t', low_memory=False)
    
    # Convert the 'Rsq' column to numeric, forcing errors to NaN
    df['Rsq'] = pd.to_numeric(df['Rsq'], errors='coerce')
    
    # Drop rows where 'Rsq' is NaN
    df = df.dropna(subset=['Rsq'])
    
    # Filter the DataFrame based on the Rsq column
    filtered_df = df[df['Rsq'] >= rsq_threshold]
    
    return filtered_df

'''
def parse_info_column(info_string, key='R2'):
    """
    Extracts the value for the given key (R2) from the INFO field.
    """
    info_pairs = [pair.split('=') for pair in info_string.split(';') if '=' in pair]
    info_dict = {k: v for k, v in info_pairs}
    return float(info_dict.get(key, 'nan'))

def read_and_filter_info(filepath, rsq_threshold=0.7):
    # Open the gzip file and read it directly into a DataFrame
    with gzip.open(filepath, 'rt') as file:
        df = pd.read_csv(file, sep='\t', comment='#', header=None, low_memory=False)
    # Define column names based on the VCF format (adjust as needed)
    df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    # Extract the R2 value from the INFO column
    df['Rsq'] = df['INFO'].apply(lambda x: parse_info_column(x, key='R2'))
    # Convert 'Rsq' to numeric, forcing errors to NaN
    df['Rsq'] = pd.to_numeric(df['Rsq'], errors='coerce')
    # Drop rows where 'Rsq' is NaN
    df = df.dropna(subset=['Rsq'])
    # Filter the DataFrame based on the Rsq column
    filtered_df = df[df['Rsq'] >= rsq_threshold]
    return filtered_df

'''

vcf_path = '/local/tsaid/trabajo/imputation/2019_Michigan_28_Oct/hla_perl/HLA_region.vcf.gz'
info_path = '/local/tsaid/trabajo/imputation/2019_Michigan_28_Oct/hla_perl/HLA_region.info.gz'

vcf_df = read_vcf(vcf_path)
filtered_info_df = read_and_filter_info(info_path)

#filtered_snps = set(filtered_info_df['SNP'])

#filtered_vcf_df = vcf_df[vcf_df['SNP'].isin(filtered_snps)]

# Step 1: Extract the positions from the 'SNP' column in filtered_info_df
#filtered_positions = filtered_info_df['ID'].apply(lambda x: int(x.split(':')[1]))
# Safely extract positions from the 'ID' column
def extract_position(id_string):
    try:
        return int(id_string.split(':')[1])
    except (IndexError, ValueError):
        return None

# Apply the function to extract positions
filtered_positions = filtered_info_df['ID'].apply(extract_position)

# Remove rows with None (invalid) positions
filtered_positions = filtered_positions.dropna().astype(int)

# Step 2: Filter the vcf_df based on the positions
filtered_vcf_df = vcf_df[vcf_df['POS'].isin(filtered_positions)]

# Display filtered VCF DataFrame
print(filtered_vcf_df.head())

def write_vcf(df, input_vcf_path, output_vcf_path):
    with gzip.open(input_vcf_path, 'rt') as infile, open(output_vcf_path, 'wt') as outfile:
        for line in infile:
            if line.startswith('##'):
                outfile.write(line)
            elif line.startswith('#'):
                outfile.write(line)
                break
        df.to_csv(outfile, sep='\t', index=False, header=False)

filtered_vcf_path = '/local/tsaid/trabajo/imputation/2019_Michigan_28_Oct/hla_perl/filtered_MICHIGAN_chr6.dose.vcf'

write_vcf(filtered_vcf_df, vcf_path, filtered_vcf_path)

print("Filtered VCF file has been written to:", filtered_vcf_path)

'''
def read_and_filter_info(filepath, rsq_threshold=0.7):
    desired_keys = {'R2'}
    
    with gzip.open(filepath, 'rt') as file:
        lines = [line.strip() for line in file if not line.startswith('##')]
        
    # Extract column names from the header line
    header_line = lines[0]
    columns = header_line.lstrip('#').split('\t')
    
    # Read the data into a DataFrame
    df = pd.read_csv(filepath, comment='#', sep='\t', names=columns, compression='gzip', low_memory=False)
    
    # Parse the INFO column
    info_data = df['INFO'].apply(lambda x: parse_info_column(x, desired_keys))
    info_df = pd.DataFrame(info_data.tolist(), index=df.index)
    
    # Merge the parsed INFO columns back with the original DataFrame
    df = pd.concat([df.drop(columns='INFO'), info_df], axis=1)
    
    # Filter the DataFrame based on the R2 column
    filtered_df = df[df['R2'] >= rsq_threshold]
    
    return filtered_df
'''
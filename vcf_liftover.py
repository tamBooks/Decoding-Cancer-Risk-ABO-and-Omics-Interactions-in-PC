import pandas as pd
from pyliftover import LiftOver

# Load the chain file for GRCh37 to GRCh38 conversion
chain_file = '/local/comun/dbs/chainFilesLiftOver/GRCh37_to_GRCh38.chain.gz'
lo = LiftOver(chain_file)

# Read the VCF file into a pandas DataFrame
vcf_file = '/local/tsaid/trabajo/vcf19/merged/merged_panGen2019.vcf'

def read_vcf(filepath):
    # Open the VCF file
    with open(filepath, 'r') as file:
        # Read the header lines
        lines = [line.strip() for line in file if not line.startswith('##')]
        # The last header line contains the column names
        header_line = lines[0]
        # Extract column names from the header line
        columns = header_line.lstrip('#').split('\t')
        # Read the data into a DataFrame
        df = pd.read_csv(filepath, comment='#', sep='\s+', names=columns, low_memory=False)
        
        return df

vcf_df = read_vcf(vcf_file)

# Display the DataFrame
print(vcf_df.head())

# Function to convert positions
def convert_position(chrom, pos):
    new_pos = lo.convert_coordinate(chrom, pos)
    if new_pos:
        return new_pos[0][1]
    else:
        return None

# Apply the conversion
vcf_df['NEW_POS'] = vcf_df.apply(lambda row: convert_position(row['CHROM'], row['POS']), axis=1)

# Filter out rows where conversion was not successful
vcf_df = vcf_df.dropna(subset=['NEW_POS'])

# Update positions to new positions
vcf_df['POS'] = vcf_df['NEW_POS'].astype(int)
vcf_df = vcf_df.drop(columns=['NEW_POS'])

# Write the updated VCF DataFrame to a new file
output_vcf_file = '/local/tsaid/trabajo/vcf19/merged/merged_panGen2019_lifted.vcf'
def write_vcf(df, vcf_file, output_path):
    # Open the original VCF file for writing
    with open(vcf_file, 'r') as f:
        header_lines = [line.strip() for line in f if line.startswith('##')]
    
    # Write the header lines back to the new VCF file
    with open(output_path, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        
        # Write the column names (header line)
        f.write('#' + '\t'.join(df.columns) + '\n')
        
        # Write the data
        for index, row in df.iterrows():
            f.write('\t'.join(map(str, row.values)) + '\n')

write_vcf(vcf_df, vcf_file, output_vcf_file)
print(f"Liftover completed. Output saved to {output_vcf_file}")
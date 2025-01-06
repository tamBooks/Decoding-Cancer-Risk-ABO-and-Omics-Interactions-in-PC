import pandas as pd

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
    df = pd.read_csv(filepath, comment='#',sep='\s+', names=columns, low_memory=False) 
    
    return df

vcf_path = '/local/tsaid/trabajo/imputation/2024_Michigan_23_Oct/abo_region.vcf'
vcf_path = '/local/tsaid/trabajo/imputation/2019_Michigan_28_Oct/abo_region.vcf'
vcf_path = '/local/tsaid/sharing/PED_TO_MAP/shapeit_onco/pangen_4SNP.vcf'
vcf_path = '/local/tsaid/sharing/PED_TO_MAP/shapeit_onco/isblac_4SNP.vcf'

df_data = read_vcf(vcf_path)

# Assigning the sample columns
sample_columns = df_data.columns[9:]
for col in sample_columns:
    df_data[col] = df_data[col].str.split(':').str[0]

# Ref 38
# rs1053878 = 133256264 --> A1 --> A
# rs2519093 = 133266456 --> A2 --> T
# rs8176743 = 133256028 --> B --> T
# rs41302905 = 133255929 --> O -->T 

#pos_A1 = [133256264] # positions for detecting blood group A
#pos_A2 = [133266456]  # positions for detecting blood group A
#pos_B = [133256028]  # positions for detecting blood group B
#pos_O = [133255929]  # positions for detecting blood group O

pos_A1 = [136131651] # positions for detecting blood group A
pos_A2 = [136141870]  # positions for detecting blood group A
pos_B = [136131415]  # positions for detecting blood group B
pos_O = [136131316]  # positions for detecting blood group O


combined_positions_A = pos_A1 + pos_A2

# Filter the dataframe using the combined list of positions
df_A = df_data[df_data['POS'].isin(combined_positions_A)]
df_B = df_data[df_data['POS'].isin(pos_B)]
df_O = df_data[df_data['POS'].isin(pos_O)]

df_sub_pos = pd.concat([df_A,df_B, df_O]).drop_duplicates().reset_index(drop=True)

# Define a function to replace genotypes with REF and ALT alleles
def replace_genotypes(row):
    ref = row['REF']
    alt = row['ALT']
    for col in df_sub_pos.columns[9:]:  # Loop over genotype columns
        genotype = row[col]
        if genotype == '0|0':
            row[col] = f'{ref}|{ref}'
        elif genotype == '0|1':
            row[col] = f'{ref}|{alt}'
        elif genotype == '1|0':
            row[col] = f'{alt}|{ref}'
        elif genotype == '1|1':
            row[col] = f'{alt}|{alt}'
    return row

# Apply the replacement function to each row in the DataFrame
df_sub_pos = df_sub_pos.apply(replace_genotypes, axis=1)

# Display the updated DataFrame
print(df_sub_pos)

# Function to determine blood type from genotypes
def determine_blood_type(row):
    # Initialize A, B, O counts
    A = B = 0
    O = 1  # Set O as the default blood type
    # Retrieve the alleles at the specified positions
    A1_allele = row[pos_A1[0]]  # Allele at pos_A1
    A2_allele = row[pos_A2[0]]  # Allele at pos_A2
    B_allele = row[pos_B[0]]    # Allele at pos_B
    O_allele = row[pos_O[0]]    # Allele at pos_O
    # Split the alleles by '|' to access both alleles (e.g. A1|A2)
    A1_allele_1, A1_allele_2 = A1_allele.split('|')
    A2_allele_1, A2_allele_2 = A2_allele.split('|')
    B1, B2 = B_allele.split('|')
    O1, O2 = O_allele.split('|')
    # Check for A alleles at pos_A1 and pos_A2
    if A1_allele_1 == 'A' or A1_allele_2 == 'A':
        A = 1  # A allele detected
    if A2_allele_1 == 'T' or A2_allele_2 == 'T':
        A = 1  # A allele detected
    if ('A' in [A1_allele_1, A1_allele_2]) and ('T' in [A2_allele_1, A2_allele_2]):
        A = 2
    # Check for B alleles
    if B1 == 'T' or B2 == 'T':
        B = 1  # B allele detected
    # Check for O alleles
    if O1 == 'T' or O2 == 'T':
        O = 1  # O allele detected
    # Additional checks for specific allele combinations
    if A1_allele == 'A|A' or A2_allele == 'T|T':
        A = 2
    if B_allele == 'T|T':
        B = 2
    if O_allele == 'T|T':
        O = 2
    # Determine blood type based on A, B, and O alleles
    if A == 2:  # A1 and A2 both present
        if 'A' in [A1_allele_1, A1_allele_2] and 'T' in [A2_allele_1, A2_allele_2]:
            blood_type = 'A1A2'
        elif 'A' in [A1_allele_1, A1_allele_2]:
            blood_type = 'A1A1'
        elif 'T' in [A2_allele_1, A2_allele_2]:
            blood_type = 'A2A2'
    elif B == 2:
        blood_type = 'BB'
    elif O == 2:
        blood_type = 'OO'
    elif A == 1 and B == 1:
        if 'A' in [A1_allele_1, A1_allele_2]:
            blood_type = 'A1B'
        else:
            blood_type = 'A2B'
    elif A == 1 and O == 1:
        if 'A' in [A1_allele_1, A1_allele_2]:
            blood_type = 'A1O'
        else:
            blood_type = 'A2O'
    elif B == 1 and O == 1:
        blood_type = 'BO'
    elif A == 1: 
        if 'A' in [A1_allele_1, A1_allele_2]:
            blood_type = 'A1'  # A1 only
        else:
            blood_type = 'A2'  # A2 only
    elif B == 1:
        blood_type = 'B'
    else:
        blood_type = 'OO'  # Default to OO
    return pd.Series([A, B, O, blood_type])

# Applying the function to the transposed DataFrame
df_blood_types = df_sub_pos.T.iloc[9:]  # Transpose to have samples as rows and remove non-sample columns
df_blood_types.columns = df_sub_pos['POS'] # Set the column names as POS values
df_blood_types[['A', 'B', 'O', 'TYPE']] = df_blood_types.apply(determine_blood_type, axis=1)
df_blood_types['TYPE'].unique()

# Add a column for dosage
df_blood_types['DOSAGE'] = df_blood_types['TYPE']

# Define the mapping from the original blood types to the simplified blood types
type_mapping = {
    'BO': 'B',
    'BB': 'B',
    'AO': 'A',
    'AA': 'A',
    'AB': 'AB',
    'OO': 'O'
}

type_mapping = {
   'OO':"O",
   'A2O':"A",
   'BO':'B',
   'A2B':'AB',
   'A1O':'A',
   'A1B':'AB',
   'A2A2':'A',
   'A2A1':'A',
   'A1A1':'A',
   'A1A2':'A',
   'BB':'B'  
}
# Apply the mapping to the 'TYPE' column
df_blood_types['TYPE'] = df_blood_types['TYPE'].replace(type_mapping)

# Reset the index to have sample names as a column
df_blood_types.reset_index(inplace=True)
df_blood_types.rename(columns={'index': 'Name'}, inplace=True)

# Display the final DataFrame
print(df_blood_types[['Sample_ID', 'A', 'B', 'O', 'DOSAGE', 'TYPE']])

# save the resulted table 
#df_blood_types = df_blood_types.rename(columns={133256264:"A_POS1",133266456:"A_POS2", 133256028:"B_POS",133255929:"O_POS"})
df_blood_types = df_blood_types.rename(columns={136131651:"A_POS1",136141870:"A_POS2", 136131415:"B_POS",136131316:"O_POS"})

output_file = '/local/tsaid/trabajo/imputation/2024_Michigan_23_Oct/ABO_summary_4_SNP.txt'
output_file = '/local/tsaid/trabajo/imputation/2019_Michigan_28_Oct/ABO_summary_4_SNP.txt'
output_file = '/local/tsaid/sharing/PED_TO_MAP/shapeit_onco/ABO_PANGEN_summary_4_SNP.txt'
output_file = '/local/tsaid/sharing/PED_TO_MAP/shapeit_onco/ABO_ISBLAC_summary_4_SNP.txt'
df_blood_types.to_csv(output_file, sep=',', index=False)

#####################
# Function to calculate allele frequency
def calculate_allele_frequency(df, column_name):
    # Count occurrences of each blood type
    counts = df[column_name].value_counts()
    # Calculate total number of records
    total = len(df)
    # Calculate frequency
    frequency = (counts / total) * 100  # Percentage
    return frequency

# Calculate allele frequencies in 'blood_type' column
blood_type_frequencies = calculate_allele_frequency(df_blood_types, 'O_POS')

print("\nTYPE Allele Frequencies:")
print(blood_type_frequencies)

# Function to calculate genotype counts
def calculate_genotype_counts(df, column_name):
    # Count occurrences of each genotype in the specified column
    genotype_counts = df[column_name].value_counts()
    return genotype_counts

# Calculate genotype counts for 'A_POS1'
genotype_counts = calculate_genotype_counts(df_blood_types, 'O_POS')

# Display the genotype counts
print("Genotype Counts:")
print(genotype_counts)

# Function to calculate allele frequency for individual alleles
def calculate_individual_allele_frequency(df, column_name):
    # Create a list to hold individual alleles
    all_alleles = []
    # Split each entry in the column (e.g., G|A into G and A)
    for genotype in df[column_name]:
        alleles = genotype.split('|')  # Split by '|' to get individual alleles
        all_alleles.extend(alleles)    # Add alleles to the list
    # Count occurrences of each allele
    allele_counts = pd.Series(all_alleles).value_counts()
    # Calculate total number of alleles (each row has 2 alleles)
    total_alleles = len(all_alleles)
    # Calculate allele frequencies as percentage
    allele_frequency = (allele_counts / total_alleles) * 100
    return allele_frequency

# Calculate allele frequencies for A_POS1
allele_frequencies = calculate_individual_allele_frequency(df_blood_types, 'O_POS')

# Display the allele frequencies
print("Frequencies:")
print(allele_frequencies)

#df_blood_types.groupby(['A_POS1', 'A_POS2','B_POS','O_POS','DOSAGE', 'TYPE']).size().reset_index(name='Count').sort_values(by='DOSAGE', ascending=True)
'''
    # Determine blood type based on A, B, and O alleles
    if A == 2:
        blood_type = 'AA'
    elif B == 2:
        blood_type = 'BB'
    elif O == 2:
        blood_type = 'OO'
    elif A == 1 and B == 1:
        blood_type = 'AB'
    elif A == 1 and O == 1:
        blood_type = 'AO'
    elif B == 1 and O == 1:
        blood_type = 'BO'
    elif A == 1: 
        blood_type = 'A' # A + blood type
    elif B == 1:
        blood_type = 'B' # B + blood type
    else:
        blood_type = 'OO'
'''

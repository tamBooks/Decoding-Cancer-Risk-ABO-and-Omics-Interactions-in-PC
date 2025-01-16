import pandas as pd

# Read txt files
ABO_text2019 = pd.read_csv('/local/tsaid/trabajo/imputation/2019_Michigan_28_Oct/ABO_summary_4_SNP.txt', sep = ',')
ABO_text2024 = pd.read_csv('/local/tsaid/trabajo/imputation/2024_Michigan_23_Oct/ABO_summary_4_SNP.txt', sep = ',')

#Concat and rename sample names
df_whole = pd.concat([ABO_text2024,ABO_text2019]).drop_duplicates().reset_index(drop=True)
df_whole['Name'] = df_whole['Name'].str.replace('^0_', '', regex=True)

# Samplesheets of ids
sample_sheet = pd.read_csv('/local/tsaid/trabajo/infomation/IDs_arrays_2019and24_tomerge_Maria.txt', sep=' ')
sample_sheet['Name'] = sample_sheet['SentrixBarcode_A'].astype(str) + '_' + sample_sheet['SentrixPosition_A'].astype(str)
sample_sheet= sample_sheet.drop(['SentrixBarcode_A','SentrixPosition_A'], axis = 1)
# Replace the values in the 'Name' column according to the `replace_names` dictionary

replace_names = {"208025880116_R09C02":"208025880116_R06C02","207769820036_R07C02":"208025880116_R08C02",
      "208026130167_R12C01":"208025880116_R10C01","208025880092_R12C01":"208026130093_R05C01",
      "208026130070_R11C02":"208026130093_R07C02","208026130070_R01C01":"208026130093_R10C01",
      "208025880110_R11C02":'208026130093_R10C02', '208025880116_R02C01':'208026130141_R02C02',
      '208026130054_R08C01':'208026130141_R07C01','208025880109_R12C01':'208026130141_R09C01',
      '208026130141_R11C01':'208026130141_R09C02','208026130093_R02C02':'208026130141_R10C01',
      '208026130172_R02C01':'208026130172_R04C02',"208026130172_R09C01":"208026130172_R10C01", 
      "208026130093_R11C01":'208026130172_R12C02',"202787460056_R06C01":'202787460069_R04C02',
      "202787460056_R10C01":"202787460069_R06C02","202787460018_R09C02":"202787460069_R11C02"
      }

sample_sheet['Name'] = sample_sheet['Name'].replace(replace_names)

#Perform the merge with 'left' join to keep all rows in df_whole
#merged_data = df_whole.merge(sample_sheet[['Sample_ID', 'Name']], on='Name', how='outer', indicator=True)
merged_data = df_whole.merge(sample_sheet[['Sample_ID','Name']],on='Name',how='outer')
merged_data = merged_data.dropna()

merged_data

# Identify rows in df_whole that didn't merge by checking for NaN in the Sample_ID column
not_merged = merged_data[merged_data['Name'].isna()]

output_file = '/local/tsaid/trabajo/imputation/merged_data_abo.txt'
merged_data.to_csv(output_file, sep=',', index=False)

##########
ABO_text_isblac = pd.read_csv('/local/tsaid/sharing/PED_TO_MAP/shapeit_onco/onco_pangen_whole.txt')
ABO_text_pangen = pd.read_csv('/local/tsaid/sharing/PED_TO_MAP/shapeit_onco/ABO_ISBLAC_summary_4_SNP.txt')
#Concat and rename sample names
df_whole = pd.concat([ABO_text_pangen,ABO_text_isblac]).reset_index(drop=True)
output_file = '/local/tsaid/sharing/PED_TO_MAP/merged_data_abo.txt'
df_whole.to_csv(output_file, sep=',', index=False)




library(dplyr)
library(tidyverse)
library(data.table)
library(readxl)

load('/local/tsaid/sharing/PED_TO_MAP/ONCOARRAY/final.data.pangen1.isblac1.epicuro1.RData')
load("/local/tsaid/trabajo/principle_comp/pca_all_trans.RData")
load("/local/tsaid/trabajo/principle_comp/pca_onco_pan.RData")
load("/local/tsaid/trabajo/principle_comp/pca_gsa24.RData")
load("/local/tsaid/trabajo/principle_comp/pca_gsa19.RData")

ls()

abo_data <- read.csv("/local/tsaid/trabajo/imputation/merged_data_abo.txt")
ABO_text_pangen_isblac <- read.csv('/local/tsaid/sharing/PED_TO_MAP/merged_data_abo.txt',sep = ',')
main_data <- read_excel("/local/tsaid/trabajo/infomation/tamara_info_tfm_1.xlsx")
sample_ids_gsa <- read.table("/local/tsaid/trabajo/infomation/edited_IDs.txt")

abo_data <- abo_data %>% select(-Name)

gsa_main <- merge(main_data, sample_ids_gsa, by = "Sample_ID")
head(gsa_main)

merged_data <- merge(gsa_main, abo_data, by = "Sample_ID")
head(merged_data)

ABO_text_pangen_isblac <- ABO_text_pangen_isblac %>%
  rename(casecontrol = caco, Country = country, 
         country=Country.collapse,agec = age, sex=gender_Phenotype)

ABO_text_pangen_isblac <- ABO_text_pangen_isblac %>%
  mutate(sex = ifelse(sex == "Female", 0, ifelse(sex == "Male", 1, sex)))

ABO_text_pangen_isblac <- ABO_text_pangen_isblac %>%
  mutate(
    Name = case_when(
      study == "PANGEN" ~ gsub("^9_", "PAN_9_", Name),
      study == "ISBLAC" ~ gsub("^9_", "IS_9_", Name),
      TRUE ~ Name # Leave unchanged if neither condition is met
    )
  )

ABO_text_pangen_isblac$EuroRegion <- ABO_text_pangen_isblac$country
ABO_text_pangen_isblac$EuroRegion[ABO_text_pangen_isblac$EuroRegion == 1] = 'Britain'
ABO_text_pangen_isblac$EuroRegion[ABO_text_pangen_isblac$EuroRegion == 2] = 'South'
ABO_text_pangen_isblac$EuroRegion[ABO_text_pangen_isblac$EuroRegion == 3] = 'North'

#ABO_text_pangen_isblac= U.K. + Ireland :1, Italy + Spain: 2, Sweden+Germany: 3

# Replace numeric values with country names for GSA
merged_data$EuroRegion <- merged_data$country
merged_data$EuroRegion[merged_data$EuroRegion == 0] = 'South'
merged_data$EuroRegion[merged_data$EuroRegion == 2] = 'North'
merged_data$EuroRegion[merged_data$EuroRegion == 3] = 'Britain'

#merged_data <- merged_data %>%
#  mutate(Country = case_when(
#    Country == 0 ~ "Spain",
#    Country == 1 ~ "England",
#    Country == 2 ~ "Germany",
#    Country == 3 ~ "Ireland",
#    Country == 4 ~ "Italy",
#    Country == 5 ~ "Sweden",
#    TRUE ~ as.character(Country)  # In case of other values, keep them as is
#  ))

ABO_text_pangen_isblac_abo <- ABO_text_pangen_isblac %>% 
  select(Name, TYPE,casecontrol,EuroRegion,agec,sex)
merged_df_abo <- merged_data %>% 
  select(Name, TYPE,casecontrol,EuroRegion ,agec,sex)

ABO_text_pangen_isblac_abo$TYPE <- as.factor(ABO_text_pangen_isblac_abo$TYPE)
merged_df_abo$TYPE <- as.factor(merged_df_abo$TYPE)

ABO_text_pangen_isblac_abo$casecontrol <- as.factor(ABO_text_pangen_isblac_abo$casecontrol)
merged_df_abo$casecontrol <- as.factor(merged_df_abo$casecontrol)

combined_df_abo <- bind_rows(merged_df_abo,ABO_text_pangen_isblac_abo)
combined_df_abo$EuroRegion <- as.factor(combined_df_abo$EuroRegion)

# Create a new factor level for TYPE, replacing non-"OO" values with "NonO" and keeping "OO" as "O"
combined_df_abo$Non_O <- combined_df_abo$TYPE
levels(combined_df_abo$Non_O) <- ifelse(levels(combined_df_abo$TYPE) == "OO", "O", "NonO")

combined_df_abo$BGroup <- combined_df_abo$TYPE
combined_df_abo <- combined_df_abo %>%
  mutate(
    BGroup = fct_recode(BGroup,
                      "A" = "A1A1",
                      "A" = "A1A2",
                      "A" = "A2A2",
                      "A" = "A1O",
                      "A" = "A2O",
                      "B" = "BB",
                      "B" = "BO",
                      "AB" = "A1B",
                      "AB" = "A2B",
                      "O" = "OO")
  )

combined_df_abo$Agrouping <- combined_df_abo$TYPE
combined_df_abo <- combined_df_abo %>%
  mutate(
    Agrouping = case_when(
      TYPE %in% c("A1O", "A1A1") ~ "A1",
      TYPE %in% c("A2O", "A2A2") ~ "A2",
      TYPE == "A1A2" ~ "A1A2",
      TYPE %in% c("A1B", "A2B") ~ TYPE,  # Retain original values for A1B and A2B
      TYPE %in% c("BO", "BB") ~ "B",
      TYPE == "OO" ~ "O",
      TRUE ~ TYPE  # Retain the original value for any other unmatched types
    )
  )

combined_df_abo$Agrouping <- factor(combined_df_abo$Agrouping)

##
# Extract the first 5 PCs
# Convert the first 5 principal components to a dataframe
pca_5 <- as.data.frame(pca_all_trans$x[, 1:5])
pca_5 <- as.data.frame(pca_gsa24_trans$x[, 1:5])
pca_5 <- as.data.frame(pca_gsa19_trans$x[, 1:5])
pca_5 <- as.data.frame(pca_onco_pan_trans$x[, 1:5])

# Ensure `pca_5` has a column for rownames to facilitate merging
pca_5 <- pca_5 %>%
  mutate(Name = rownames(pca_5))

# Merge the dataframes based on the "Name" column
merged_data <- merge(pca_5, combined_df_abo, by = "Name", all.x =TRUE)

# Ensure the order of rows matches the original `pca_5`
merged_data <- merged_data[match(rownames(pca_5), merged_data$Name), ]
rownames(merged_data) <- NULL

merged_data$arraytype <- c(rep("GSA24", 313), rep("GSA19", 151), rep("ONCO",2017))
merged_data$arraytype <-as.factor(merged_data$arraytype)
whole_dataframe <- merged_data
save(whole_dataframe,file ="/local/tsaid/trabajo/infomation/whole_dataframe_pc.RData")

merged_data$arraytype <- rep("GSA24", 313)
merged_data$arraytype <-as.factor(merged_data$arraytype)
gsa2024_data <- merged_data
save(gsa2024_data,file ="/local/tsaid/trabajo/infomation/gsa2024_dataframe_pc.RData")

merged_data$arraytype <- rep("GSA19", 151)
merged_data$arraytype <-as.factor(merged_data$arraytype)
gsa2019_data <- merged_data
save(gsa2019_data,file ="/local/tsaid/trabajo/infomation/gsa2019_dataframe_pc.RData")

merged_data$arraytype <- rep("ONCO", 2017)
merged_data$arraytype <-as.factor(merged_data$arraytype)
onco_data <- merged_data
save(onco_data,file ="/local/tsaid/trabajo/infomation/onco_dataframe_pc.RData")
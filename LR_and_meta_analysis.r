# 0. Install Packagecs if necessary
#install.packagecs("dplyr")
#install.packagecs("readxl")
#install.packagecs("caret")
#install.packages("forcats")
#install.packages("webshot")

#library(devtools)
#install_version("gtsummary", version = "1.5.0")
#devtools::install_github("ddsjoberg/gtforester")

# 1. Load the necessary libraries:
library(caret)
library(dplyr)
library(readxl)
library(tidyr)
library(forcats)
library(gtsummary)
library(broom.helpers)
library(survival)
library(gt)
library(webshot2)
library(gtforester)
library(broom)
library(ggplot2)
library(metafor)
library(stringr)
library('survival')
library('compareGroups')
library(flextable)

# 2. Load the data:
load("~/Downloads/tfm/whole_dataframe_pc.RData")
load("~/Downloads/tfm/onco_dataframe_pc.RData")
load("~/Downloads/tfm/gsa2024_dataframe_pc.RData")
load("~/Downloads/tfm/gsa2019_dataframe_pc.RData")

hla_gsa = read.csv("~/Downloads/tfm/hla_2019_2024_modified_2.txt", sep='\t')
hla_onco = read.csv("~/Downloads/tfm/oncoarray_hla_2.txt", sep='\t')
hla_onco_ids = read.table("~/Downloads/tfm/ID_oncoarray.txt")

microbiome = read.csv("~/Downloads/tfm/gsa2019_Oral_Fecal.csv", sep=";")
microbiome_gut <- microbiome[grep('ST', microbiome$MMPCID),]
microbiome_gut_whole <- merge(gsa2019_data, microbiome_gut , by="Name")

# 2.5. Modify the column Name
num_rows <- nrow(hla_gsa)  # Total number of rows in the dataframe
# Add prefixes based on index
hla_gsa$Name <- ifelse(
  seq_len(num_rows) <= 151,  # Check if the row index is within the first 151
  paste0("GSA19_", hla_gsa$Name),  # Prefix with GSA19_
  paste0("GSA24_", hla_gsa$Name)   # Prefix with GSA24_ for the rest
)

names(hla_onco)[names(hla_onco) == 'Name'] <- 'ID_Pac'
hla_onco <- merge(hla_onco,hla_onco_ids, by = 'ID_Pac')

#3. HLA allele filteration and merger 
hla_gsa <- hla_gsa %>% mutate(Name = as.character(Name))
hla_onco <- hla_onco %>% mutate(Name = as.character(Name))

# Get the names of columns that are common between hla_gsa and hla_onco
common_cols <- intersect(names(hla_gsa), names(hla_onco))

# Subset hla_onco to include only the common columns
hla_onco_subset <- hla_onco %>% select(all_of(common_cols))
hla_gsa_subset <- hla_gsa %>% select(all_of(common_cols))

# Bind the rows, now both data frames have only the common columns
hla_com_whole <- bind_rows(hla_gsa_subset, hla_onco_subset)

# Select columns
hla_cols <- hla_com_whole[, 2:ncol(hla_com_whole)]  # Adjust if needed for actual columns 7 to 108

# Calculate allele frequencies
allele_frequencies <- colSums(hla_cols) / (2 * nrow(hla_cols)) * 100

# Convert to data frame and set column names properly
allele_freq_df <- data.frame(HLA_Allele = names(allele_frequencies), Frequency = allele_frequencies)

# Order by Frequency in ascending order
ordered_allele_freq_df <- allele_freq_df[order(allele_freq_df$Frequency), ]
print(ordered_allele_freq_df)

low_freq_alleles <- ordered_allele_freq_df$HLA_Allele[ordered_allele_freq_df$Frequency < 1]
high_freq_alleles <- ordered_allele_freq_df$HLA_Allele[ordered_allele_freq_df$Frequency > 1]

# Filter `hla_com_whole` to retain only high frequency alleles and the "Name" column
hla_com_filtered <- hla_com_whole[, colnames(hla_com_whole) %in% c("Name", high_freq_alleles)]

# 4. Convert categorical variables for PDAC case and control:
## The dataframe shows in the casecontrol that there is Case, Control and Pancrititis, so divide the data to those
## Filter out the case and control
filtered_data_pdac <- subset(whole_dataframe, casecontrol %in% c(0, 1))

filtered_data_pancreatitis <- subset(whole_dataframe, casecontrol %in% c(0, 2))

### Filter out rows where 'sex' is 2
filtered_data_pdac <- filtered_data_pdac %>% 
  filter(sex != 2) # %>%

# Relevel for default
filtered_data_pdac$casecontrol <- relevel(filtered_data_pdac$casecontrol, ref = '0')
filtered_data_pdac$casecontrol<- droplevels(filtered_data_pdac$casecontrol)
filtered_data_pdac$TYPE <- relevel(filtered_data_pdac$TYPE, ref = "OO")
#filtered_data_pdac$arraytype <- relevel(filtered_data_pdac$arraytype, ref = "ONCO")
filtered_data_pdac$Non_O <- relevel(filtered_data_pdac$Non_O, ref = "O")
filtered_data_pdac$BGroup <- droplevels(filtered_data_pdac$BGroup)
filtered_data_pdac$BGroup <- relevel(filtered_data_pdac$BGroup, ref = "O")
filtered_data_pdac$Agrouping <- droplevels(filtered_data_pdac$Agrouping)
filtered_data_pdac$Agrouping <- relevel(filtered_data_pdac$Agrouping, ref = "O")
filtered_data_pdac$EuroRegion <- relevel(filtered_data_pdac$EuroRegion, ref = "South")

# Perform the left join to merge with `filtered_data_pdac`, keeping "Name" and high frequency alleles
combined_df <- left_join(filtered_data_pdac, hla_com_filtered, by = "Name")

hla_micro_abo <- left_join(combined_df,microbiome_gut,by = "Name")

#write.table(combined_df, file=paste0('~/Downloads/whole_dataset_epi_abo_hla.txt'), col.names = T, quote=FALSE, row.names= F, sep='\t')

percentage_df <- filtered_data_pdac %>%
  count(TYPE) %>%
  transmute(TYPE, percentage = n / sum(n) * 100)

bloodtype_table<-filtered_data_pdac %>% select(Non_O,BGroup,TYPE) %>% tbl_summary()

bloodtype_table |>
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/blood_table.pdf")

# table to summarize everything
hla_micro_abo %>% 
  select(arraytype ,casecontrol, agec, sex, BGroup,q.0,q.1) %>% 
  tbl_summary(by = casecontrol, 
              label = agec ~ "Age (years)",
              percent = "row"
  )

# compute a time-to-overall death variable
Clinical = hla_micro_abo[,c('agec','sex','BGroup','arraytype','casecontrol')]
tab <- descrTable(arraytype ~ . , data = Clinical[,c('agec','sex','BGroup','casecontrol')], hide.no="no")
tab


tables <- compareGroups(casecontrol ~ arraytype + agec+ sex+ BGroup+q.0+q.1,data = hla_micro_abo)
all_tables <- createTable(tables,show.p.overall = FALSE)
table_19 <- createTable(update(tables, subset=arraytype=='GSA19'), show.p.overall = FALSE)
table_24 <- createTable(update(tables, subset=arraytype=='GSA24'), show.p.overall = FALSE)
table_onco <- createTable(update(tables, subset=arraytype=='ONCO'), show.p.overall = FALSE)

export2word(all_tables, file = "~/Downloads/tfm/pooled/table.docx")
export2word(table_19, file = "~/Downloads/tfm/gsa19/table.docx")
export2word(table_24, file = "~/Downloads/tfm/gsa24/table.docx")
export2word(table_onco, file = "~/Downloads/tfm/onco/table.docx")

# Extract the row descriptions (row names) from the 'descr' slot
row_names_19 <- rownames(table_19$descr)
row_names_24 <- rownames(table_24$descr)
row_names_onco <- rownames(table_onco$descr)
row_names_all <- rownames(all_tables$descr)

# Find the common row names
common_rows <- Reduce(intersect, list(row_names_19, row_names_24, row_names_onco, row_names_all))

# Subset each table to include only the common rows
table_19_aligned <- table_19$descr[common_rows, , drop = FALSE]
table_24_aligned <- table_24$descr[common_rows, , drop = FALSE]
table_onco_aligned <- table_onco$descr[common_rows, , drop = FALSE]
all_tables_aligned <- all_tables$descr[common_rows, , drop = FALSE]

# Combine the tables
combined_table <- cbind(
  "GSA19" = table_19_aligned,
  "GSA24" = table_24_aligned,
  "ONCO" = table_onco_aligned,
  "Pooled" = all_tables_aligned
)

# Convert `combined_table` to a data frame
combined_table_df <- as.data.frame(combined_table)

# Add row names as the first column
combined_table_df <- cbind(Characteristic = rownames(combined_table_df), combined_table_df)

# Fix column names for uniqueness
colnames(combined_table_df) <- c(
  "Characteristic", 
  "GSA19_Control", "GSA19_Case",
  "GSA24_Control", "GSA24_Case",
  "Oncoarray_Control", "Oncoarray_Case",
  "Pooled_Control", "Pooled_Case"
)

# Create the `gt` table
gt_table <- gt(combined_table_df) %>%
  # Add the main headers
  tab_spanner(
    label = "GSA 2019",
    columns = c("GSA19_Control", "GSA19_Case")
  ) %>%
  tab_spanner(
    label = "GSA 2024",
    columns = c("GSA24_Control", "GSA24_Case")
  ) %>%
  tab_spanner(
    label = "Oncoarray",
    columns = c("Oncoarray_Control", "Oncoarray_Case")
  ) %>%
  tab_spanner(
    label = "Pooled",
    columns = c("Pooled_Control", "Pooled_Case")
  ) %>%
  # Add sub-headers
  cols_label(
    GSA19_Control = "Control (N=58)", GSA19_Case = "Case (N=65)",
    GSA24_Control = "Control (N=218)", GSA24_Case = "Case (N=80)",
    Oncoarray_Control = "Control (N=700)", Oncoarray_Case = "Case (N=1317)",
    Pooled_Control = "Control (N=976)", Pooled_Case = "Case (N=1462)"
  ) %>%
  # Add a title
  tab_header(
    title = "Summary Descriptive Table"
  ) %>%
  # Optional: Adjust table appearance
  tab_options(
    table.font.size = 12,
    heading.align = "center",
    column_labels.border.top.width = px(2),
    column_labels.border.bottom.width = px(2),
    #table.border.width = px(1),
    data_row.padding = px(4)
  )

# Print the table
print(gt_table)

pdf('/local/tsaid/trabajo/logistic/EuroRegionPlot.pdf')
# Create a scatter plot
ggplot(filtered_data_pdac, aes(x = PC1, y = PC2, color = as.factor(EuroRegion))) +
  geom_point(size = 0.5, alpha = 0.4) +
  labs(
    title = "PC1 vs PC2 by EuroRegion",
    x = "Principal Component 1 (PC1)",
    y = "Principal Component 2 (PC2)",
    color = "EuroRegion"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

ggplot(filtered_data_pdac, aes(x = PC1, y = PC3, color = as.factor(EuroRegion))) +
  geom_point(size = 0.5, alpha = 0.4) +
  labs(
    title = "PC1 vs PC3 by EuroRegion",
    x = "Principal Component 1 (PC1)",
    y = "Principal Component 3 (PC3)",
    color = "EuroRegion"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

ggplot(filtered_data_pdac, aes(x = PC2, y = PC3, color = as.factor(EuroRegion))) +
  geom_point(size = 0.5, alpha = 0.4) +
  labs(
    title = "PC2 vs PC3 by EuroRegion",
    x = "Principal Component 2 (PC2)",
    y = "Principal Component 3 (PC3)",
    color = "EuroRegion"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

dev.off()

# 5. Univariate Logistic Regression:
## To assess the crude association between ABO blood group variants and PDAC
### Comparison with Alleles
univariate_model1 <- glm(casecontrol ~ TYPE, data = filtered_data_pdac, family=binomial(logit))
summary(univariate_model1)
 
# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
univariate_or1 <- exp(cbind(OR = coef(univariate_model1), confint(univariate_model1)))
univariate_or1

### Comparison with Non-O
univariate_model2 <- glm(casecontrol ~ Non_O, data = filtered_data_pdac, family=binomial(logit))
summary(univariate_model2)
 
# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
univariate_or2 <- exp(cbind(OR = coef(univariate_model2), confint(univariate_model2)))
univariate_or2

### Comparison with compressed Types
univariate_model3 <- glm(casecontrol ~ BGroup, data = filtered_data_pdac, family=binomial(logit))
summary(univariate_model3)
 
# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
univariate_or3 <- exp(cbind(OR = coef(univariate_model3), confint(univariate_model3)))
univariate_or3

### Comparison with compressed Types
univariate_model4 <- glm(casecontrol ~ Agrouping, data = filtered_data_pdac, family=binomial(logit))
summary(univariate_model4)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
univariate_or4 <- exp(cbind(OR = coef(univariate_model4), confint(univariate_model4)))
univariate_or4

# Table combining the univarient models togther
table1 <-  tbl_regression(univariate_model1, exponentiate = TRUE)
table2 <-  tbl_regression(univariate_model3, exponentiate = TRUE)
table3 <-  tbl_regression(univariate_model2, exponentiate = TRUE)
table4 <-  tbl_regression(univariate_model4, exponentiate = TRUE)

tbl_merge1 <-
  tbl_merge(
    tbls = list(table1, table2, table3, table4),
    tab_spanner = c("**Allele**", "**Blood Group**", "**Non-O**", "**A-Alleles Grouped**")
  )
tbl_merge1 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/univarient_tables.pdf")

forest_plot1 <- 
  table1 %>%
  modify_column_merge(
    pattern = "{estimate} (95% CI {ci}; {p.value})",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model 1**"))

forest_plot1 |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/uni_forest1.pdf")

forest_plot2 <- 
  table2 %>%
  modify_column_merge(
    pattern = "{estimate} (95% CI {ci}; {p.value})",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model 2**"))

forest_plot2 |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/uni_forest2.pdf")

forest_plot3 <- 
  table3 %>%
  modify_column_merge(
    pattern = "{estimate} (95% CI {ci}; {p.value})",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model 3**"))

forest_plot3 |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/uni_forest3.pdf")

forest_plot4 <- 
  table4 %>%
  modify_column_merge(
    pattern = "{estimate} (95% CI {ci}; {p.value})",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model 4**"))

forest_plot4 |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/uni_forest4.pdf")

### 6. Multivariate Logistic Regression:
#To adjust for potential confounders like agec, sex, and smoking status, you can run a multivariate logistic regression.
## Pooled Analysis
multivariate_model1 <- glm(casecontrol ~ TYPE + agec + sex  + PC1 + PC2 + PC3 + PC4 + PC5, 
                                     data = filtered_data_pdac, binomial(logit))
summary(multivariate_model1)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or1 <- exp(cbind(OR = coef(multivariate_model1), confint(multivariate_model1)))
multivariate_or1

multivariate_model2 <- glm(casecontrol ~ Non_O + agec + sex + PC1 + PC2 + PC3 + PC4 + PC5, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model2)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or2 <- exp(cbind(OR = coef(multivariate_model2), confint(multivariate_model2)))
multivariate_or2

multivariate_model3 <- glm(casecontrol ~ BGroup + agec + sex + PC1 + PC2 + PC3 + PC4 + PC5, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model3)
# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or3 <- exp(cbind(OR = coef(multivariate_model3), confint(multivariate_model3)))
multivariate_or3

multivariate_model4 <- glm(casecontrol ~ Agrouping + agec + sex + PC1 + PC2 + PC3 + PC4 + PC5, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model3)
# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or4 <- exp(cbind(OR = coef(multivariate_model4), confint(multivariate_model4)))
multivariate_or4

## Plot
# Define each regression table with specific data subsets and labels
table1 <- tbl_regression(multivariate_model1, exponentiate = TRUE)# %>%
#modify_header(label ~ "**Allele**")  # Label for Allele model
table2 <- tbl_regression(multivariate_model3, exponentiate = TRUE) #%>%
#modify_header(label ~ "**Non-O**")   # Label for Non-O model
table3 <- tbl_regression(multivariate_model2, exponentiate = TRUE) #%>%
#modify_header(label ~ "**Blood Group**")  # Label for Blood Group model
table4 <- tbl_regression(multivariate_model4, exponentiate = TRUE) #%>%
#modify_header(label ~ "**A-Allele Grouping**")  # Label for Blood Group model

table1 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/no_euro_multivarient_table1.pdf")

table2 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/no_euro_multivarient_table2.pdf")

table3 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/no_euro_multivarient_table3.pdf")

table4 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/no_euro_multivarient_table4.pdf")

forest_plot1 <- 
  table1 %>%
  modify_column_merge(
    pattern = "{estimate}; p-value: {p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 1**"))

forest_plot1 |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/no_euro_multi_forest1.pdf")

forest_plot2 <- 
  table2 %>%
  modify_column_merge(
    pattern = "{estimate}; p-value: {p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 2**"))

forest_plot2 |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/no_euro_multi_forest2.pdf")

forest_plot3 <- 
  table3 %>%
  modify_column_merge(
    pattern = "{estimate}; p-value: {p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 3**"))

forest_plot3 |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/no_euro_multi_forest3.pdf")

forest_plot4 <- 
  table4 %>%
  modify_column_merge(
    pattern = "{estimate}; p-value: {p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 4**"))

forest_plot4 |> 
  gt::gtsave(filename = "~/Downloads/tfm/pooled/no_euro_multi_forest4.pdf")

# Merge the tables
# To hide repeated demographic labels, we can leverage gtsummary's modify_* functions
tbl_merge2 <- tbl_merge(
  tbls = list(table1, table2, table3, table4),
  tab_spanner = c("Allele", "Non-O", "Blood Group", "A-Grouping") 
)%>%
  modify_spanning_header(everything() ~ NA_character_) %>%  # Remove existing spanning headers if needed
  modify_column_unhide(columns = c("label", "estimate_1", "ci_1", 
                                   "estimate_2", "ci_2", "estimate_3", "ci_3", 
                                   "estimate_4", "ci_4"))
# Convert to gt format if needed
gt_table <- as_gt(tbl_merge2)

# Save or display the merged table
gtsave(gt_table, filename = "~/Downloads/tfm/pooled/no_euro_multivarient_tables.html")

### GSA24
filtered_data_pdac <- subset(gsa2024_data, casecontrol %in% c(0, 1))
filtered_data_pancreatitis <- subset(gsa2024_data, casecontrol %in% c(0, 2))

### Filter out rows where 'sex' is 2
filtered_data_pdac <- filtered_data_pdac %>% 
  filter(sex != 2) # %>%

# Relevel for default
filtered_data_pdac$casecontrol <- relevel(filtered_data_pdac$casecontrol, ref = '0')
filtered_data_pdac$TYPE <- relevel(filtered_data_pdac$TYPE, ref = "OO")
#filtered_data_pdac$arraytype <- relevel(filtered_data_pdac$arraytype, ref = "ONCO")
filtered_data_pdac$Non_O <- relevel(filtered_data_pdac$Non_O, ref = "O")
filtered_data_pdac$BGroup <- relevel(filtered_data_pdac$BGroup, ref = "O")
filtered_data_pdac$Agrouping <- relevel(filtered_data_pdac$Agrouping, ref = "O")
filtered_data_pdac$EuroRegion <- relevel(filtered_data_pdac$EuroRegion, ref = "South")

bloodtype_table <-filtered_data_pdac %>% select(TYPE,BGroup,Agrouping,Non_O) %>% tbl_summary()
bloodtype_table |>
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/blood_table.pdf")


multivariate_model_gsa24_type <- glm(casecontrol ~ TYPE + agec + sex + PC1 +PC2, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model_gsa24_type)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or1_gsa24_type <- exp(cbind(OR = coef(multivariate_model_gsa24_type), confint(multivariate_model_gsa24_type)))
multivariate_or1_gsa24_type

multivariate_model2_gsa24_non <- glm(casecontrol ~ Non_O + agec + sex + PC1 +PC2, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model2_gsa24_non)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or2_gsa24_non <- exp(cbind(OR = coef(multivariate_model2_gsa24_non), 
                              confint(multivariate_model2_gsa24_non)))
multivariate_or2_gsa24_non

multivariate_model3_gsa24_group <- glm(casecontrol ~ BGroup + agec + sex + PC1 +PC2, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model3_gsa24_group)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or3_gsa24_group <- exp(cbind(OR = coef(multivariate_model3_gsa24_group), 
                              confint(multivariate_model3_gsa24_group)))
multivariate_or3_gsa24_group

multivariate_model4_gsa24_group <- glm(casecontrol ~ Agrouping + agec + sex  + PC1 +PC2, 
                                       data = filtered_data_pdac, binomial(logit))
summary(multivariate_model4_gsa24_group)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or4_gsa24_group <- exp(cbind(OR = coef(multivariate_model4_gsa24_group), 
                                          confint(multivariate_model4_gsa24_group)))
multivariate_or4_gsa24_group

## Plot
# Define each regression table with specific data subsets and labels
table1 <- tbl_regression(multivariate_model_gsa24_type, exponentiate = TRUE)# %>%
#modify_header(label ~ "**Allele**")  # Label for Allele model
table2 <- tbl_regression(multivariate_model3_gsa24_group, exponentiate = TRUE) #%>%
#modify_header(label ~ "**Non-O**")   # Label for Non-O model
table3 <- tbl_regression(multivariate_model2_gsa24_non, exponentiate = TRUE) #%>%
#modify_header(label ~ "**Blood Group**")  # Label for Blood Group model
table4 <- tbl_regression(multivariate_model4_gsa24_group, exponentiate = TRUE) #%>%
#modify_header(label ~ "**A-Allele Grouping**")  # Label for Blood Group model

table1 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/gsa24_multivarient_table1.pdf")

table2 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/gsa24_multivarient_table2.pdf")

table3 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/gsa24_multivarient_table3.pdf")

table4 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/gsa24_multivarient_table4.pdf")

forest_plot1 <- 
  table1 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 1**"))

forest_plot1 |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/gsa24_multi_forest1.pdf")

forest_plot2 <- 
  table2 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 2**"))

forest_plot2 |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/gsa24_multi_forest2.pdf")

forest_plot3 <- 
  table3 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 3**"))

forest_plot3 |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/gsa24_multi_forest3.pdf")

forest_plot4 <- 
  table4 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 4**"))

forest_plot4 |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa24/gsa24_multi_forest4.pdf")

# Merge the tables
# To hide repeated demographic labels, we can leverage gtsummary's modify_* functions
tbl_merge2 <- tbl_merge(
  tbls = list(table1, table2, table3, table4),
  tab_spanner = c("Allele", "Non-O", "Blood Group", "A-Grouping") 
)%>%
  modify_spanning_header(everything() ~ NA_character_) %>%  # Remove existing spanning headers if needed
  modify_column_unhide(columns = c("label", "estimate_1", "ci_1", 
                                   "estimate_2", "ci_2", "estimate_3", "ci_3", 
                                   "estimate_4", "ci_4"))

# Convert to gt format if needed
gt_table <- as_gt(tbl_merge2)

# Save or display the merged table
gtsave(gt_table, filename = "~/Downloads/tfm/gsa24/gsa24_multivarient_tables.html")


# GSA19 
filtered_data_pdac <- subset(gsa2019_data, casecontrol %in% c(0, 1))
filtered_data_pancreatitis <- subset(gsa2019_data, casecontrol %in% c(0, 2))

### Filter out rows where 'sex' is 2
filtered_data_pdac <- filtered_data_pdac %>% 
  filter(sex != 2) # %>%

# Relevel for default
filtered_data_pdac$casecontrol <- relevel(filtered_data_pdac$casecontrol, ref = '0')
filtered_data_pdac$TYPE <- relevel(filtered_data_pdac$TYPE, ref = "OO")
#filtered_data_pdac$arraytype <- relevel(filtered_data_pdac$arraytype, ref = "ONCO")
filtered_data_pdac$Non_O <- relevel(filtered_data_pdac$Non_O, ref = "O")
filtered_data_pdac$BGroup <- relevel(filtered_data_pdac$BGroup, ref = "O")
filtered_data_pdac$Agrouping <- relevel(filtered_data_pdac$Agrouping, ref = "O")
filtered_data_pdac$EuroRegion <- relevel(filtered_data_pdac$EuroRegion, ref = "South")

bloodtype_table <-filtered_data_pdac %>% select(TYPE,BGroup,Agrouping,Non_O) %>% tbl_summary()
bloodtype_table |>
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/blood_table.pdf")

multivariate_model_gsa19_type <- glm(casecontrol ~ TYPE + agec + sex + PC1 +PC2 +PC3 + PC5+PC5, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model_gsa19_type)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or1_gsa19_type <- exp(cbind(OR = coef(multivariate_model_gsa19_type), confint(multivariate_model_gsa19_type)))
multivariate_or1_gsa19_type

multivariate_model2_gsa19_non <- glm(casecontrol ~ Non_O + agec + sex +PC1 +PC2 +PC3 + PC5+PC5 , 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model2_gsa19_non)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or2_gsa19_non <- exp(cbind(OR = coef(multivariate_model2_gsa19_non), 
                              confint(multivariate_model2_gsa19_non)))
multivariate_or2_gsa19_non

multivariate_model3_gsa19_group <- glm(casecontrol ~ BGroup + agec + sex + PC1 +PC2 +PC3 + PC5+PC5, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model3_gsa19_group)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or3_gsa19_group <- exp(cbind(OR = coef(multivariate_model3_gsa19_group), 
                              confint(multivariate_model3_gsa19_group)))
multivariate_or3_gsa19_group

multivariate_model4_gsa19_group <- glm(casecontrol ~ Agrouping + agec + sex+ PC1 +PC2 +PC3 + PC5+PC5, 
                                       data = filtered_data_pdac, binomial(logit))
summary(multivariate_model4_gsa19_group)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or4_gsa19_group <- exp(cbind(OR = coef(multivariate_model4_gsa19_group), 
                                          confint(multivariate_model4_gsa19_group)))
multivariate_or4_gsa19_group

## Plot
# Define each regression table with specific data subsets and labels
table1 <- tbl_regression(multivariate_model_gsa19_type, exponentiate = TRUE)# %>%
#modify_header(label ~ "**Allele**")  # Label for Allele model
table2 <- tbl_regression(multivariate_model3_gsa19_group, exponentiate = TRUE) #%>%
#modify_header(label ~ "**Non-O**")   # Label for Non-O model
table3 <- tbl_regression(multivariate_model2_gsa19_non, exponentiate = TRUE) #%>%
#modify_header(label ~ "**Blood Group**")  # Label for Blood Group model
table4 <- tbl_regression(multivariate_model4_gsa19_group, exponentiate = TRUE) #%>%
#modify_header(label ~ "**A-Allele Grouping**")  # Label for Blood Group model

table1 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/gsa19_multivarient_table1.pdf")

table2 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/gsa19_multivarient_table2.pdf")

table3 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/gsa19_multivarient_table3.pdf")

table4 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/gsa19_multivarient_table4.pdf")

forest_plot1 <- 
  table1 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 1**"))

forest_plot1 |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/gsa19_multi_forest1.pdf")

forest_plot2 <- 
  table2 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 2**"))

forest_plot2 |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/gsa19_multi_forest2.pdf")

forest_plot3 <- 
  table3 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 3**"))

forest_plot3 |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/gsa19_multi_forest3.pdf")

forest_plot4 <- 
  table4 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 4**"))

forest_plot4 |> 
  gt::gtsave(filename = "~/Downloads/tfm/gsa19/gsa19_multi_forest4.pdf")

# Merge the tables
# To hide repeated demographic labels, we can leverage gtsummary's modify_* functions
tbl_merge2 <- tbl_merge(
  tbls = list(table1, table2, table3, table4),
  tab_spanner = c("Allele", "Non-O", "Blood Group", "A-Grouping") 
)%>%
  modify_spanning_header(everything() ~ NA_character_) %>%  # Remove existing spanning headers if needed
  modify_column_unhide(columns = c("label", "estimate_1", "ci_1", 
                                   "estimate_2", "ci_2", "estimate_3", "ci_3", 
                                   "estimate_4", "ci_4"))

# Convert to gt format if needed
gt_table <- as_gt(tbl_merge2)

# Save or display the merged table
gtsave(gt_table, filename = "~/Downloads/tfm/gsa19/gsa19_multivarient_tables.html")


# ONCO
filtered_data_pdac <- subset(onco_data, casecontrol %in% c(0, 1))
filtered_data_pancreatitis <- subset(onco_data, casecontrol %in% c(0, 2))

### Filter out rows where 'sex' is 2
filtered_data_pdac <- filtered_data_pdac %>% 
  filter(sex != 2) # %>%

# Relevel for default
filtered_data_pdac$casecontrol <- relevel(filtered_data_pdac$casecontrol, ref = '0')
filtered_data_pdac$TYPE <- relevel(filtered_data_pdac$TYPE, ref = "OO")
#filtered_data_pdac$arraytype <- relevel(filtered_data_pdac$arraytype, ref = "ONCO")
filtered_data_pdac$Non_O <- relevel(filtered_data_pdac$Non_O, ref = "O")
filtered_data_pdac$BGroup <- relevel(filtered_data_pdac$BGroup, ref = "O")
filtered_data_pdac$Agrouping <- relevel(filtered_data_pdac$Agrouping, ref = "O")
filtered_data_pdac$EuroRegion <- as.factor(filtered_data_pdac$EuroRegion)
filtered_data_pdac$EuroRegion <- relevel(filtered_data_pdac$EuroRegion, ref = "South")

bloodtype_table <-filtered_data_pdac %>% select(TYPE,BGroup,Agrouping,Non_O) %>% tbl_summary()
bloodtype_table |>
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/blood_table.pdf")

multivariate_model_onco_type <- glm(casecontrol ~ TYPE + agec + sex + PC1+PC2 , 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model_onco_type)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or1_onco_type <- exp(cbind(OR = coef(multivariate_model_onco_type), confint(multivariate_model_onco_type)))
multivariate_or1_onco_type

multivariate_model2_onco_non <- glm(casecontrol ~ Non_O + agec + sex+ PC1+PC2 , 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model2_onco_non)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or2_onco_non <- exp(cbind(OR = coef(multivariate_model2_onco_non), 
                              confint(multivariate_model2_onco_non)))
multivariate_or2_onco_non

multivariate_model3_onco_group <- glm(casecontrol ~ BGroup + agec + sex + PC1+PC2, 
                           data = filtered_data_pdac, binomial(logit))
summary(multivariate_model3_onco_group)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or3_onco_group <- exp(cbind(OR = coef(multivariate_model3_onco_group), 
                              confint(multivariate_model3_onco_group)))
multivariate_or3_onco_group

multivariate_model4_onco_group <- glm(casecontrol ~ Agrouping + agec + sex + PC1+PC2, 
                                      data = filtered_data_pdac, binomial(logit))
summary(multivariate_model4_onco_group)

# Calculate Odds Ratios (ORs) and 95% Confidence Intervals (CIs)
multivariate_or4_onco_group <- exp(cbind(OR = coef(multivariate_model4_onco_group), 
                                         confint(multivariate_model4_onco_group)))
multivariate_or4_onco_group

## Plot
# Define each regression table with specific data subsets and labels
table1 <- tbl_regression(multivariate_model_onco_type, exponentiate = TRUE)# %>%
#modify_header(label ~ "**Allele**")  # Label for Allele model
table2 <- tbl_regression(multivariate_model3_onco_group, exponentiate = TRUE) #%>%
#modify_header(label ~ "**Non-O**")   # Label for Non-O model
table3 <- tbl_regression(multivariate_model2_onco_non, exponentiate = TRUE) #%>%
#modify_header(label ~ "**Blood Group**")  # Label for Blood Group model
table4 <- tbl_regression(multivariate_model4_onco_group, exponentiate = TRUE) #%>%
#modify_header(label ~ "**A-Allele Grouping**")  # Label for Blood Group model

table1 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/onco_multivarient_table1.pdf")

table2 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/onco_multivarient_table2.pdf")

table3 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/onco_multivarient_table3.pdf")

table4 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/onco_multivarient_table4.pdf")

forest_plot1 <- 
  table1 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 1**"))

forest_plot1 |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/onco_multi_forest1.pdf")

forest_plot2 <- 
  table2 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 2**"))

forest_plot2 |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/onco_multi_forest2.pdf")

forest_plot3 <- 
  table3 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 3**"))

forest_plot3 |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/onco_multi_forest3.pdf")

forest_plot4 <- 
  table4 %>%
  modify_column_merge(
    pattern = "{estimate};  p-value:{p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 4**"))

forest_plot4 |> 
  gt::gtsave(filename = "~/Downloads/tfm/onco/onco_multi_forest4.pdf")

# Merge the tables
# To hide repeated demographic labels, we can leverage gtsummary's modify_* functions
tbl_merge2 <- tbl_merge(
  tbls = list(table1, table2, table3, table4),
  tab_spanner = c("Allele", "Non-O", "Blood Group", "A-Grouping") 
)%>%
  modify_spanning_header(everything() ~ NA_character_) %>%  # Remove existing spanning headers if needed
  modify_column_unhide(columns = c("label", "estimate_1", "ci_1", 
                                   "estimate_2", "ci_2", "estimate_3", "ci_3", 
                                   "estimate_4", "ci_4"))

# Convert to gt format if needed
gt_table <- as_gt(tbl_merge2)

# Save or display the merged table
gtsave(gt_table, filename = "~/Downloads/tfm/onco/onco_multivarient_tables.html")


## VIF 
#install.packages("car")
library("car")
data_vif <- filtered_data_pdac

data_vif$casecontrol<- as.numeric(filtered_data_pdac$casecontrol)
data_vif$PC1<- as.numeric(filtered_data_pdac$PC1)
data_vif$PC2<- as.numeric(filtered_data_pdac$PC2)
data_vif$PC3<- as.numeric(filtered_data_pdac$PC3)
data_vif$PC4<- as.numeric(filtered_data_pdac$PC4)
data_vif$PC5<- as.numeric(filtered_data_pdac$PC5)
data_vif$TYPE <- as.numeric(filtered_data_pdac$TYPE)
data_vif$agec <- as.numeric(filtered_data_pdac$agec)
data_vif$sex <- as.numeric(filtered_data_pdac$agec)

lm_model <- lm(casecontrol ~ TYPE + agec + sex  + PC1 + PC2 + PC3 + PC4 + PC5, 
               data = data_vif)
# Calculate VIF
vif_values <- vif(lm_model)
# Print the VIF values
print(vif_values)

# Extract logOR and SE for each population
logOR_gsa24 <- coef(multivariate_model2_gsa24_non)["Non_ONonO"]
SE_gsa24 <- sqrt(vcov(multivariate_model2_gsa24_non)["Non_ONonO", "Non_ONonO"])

logOR_gsa19 <- coef(multivariate_model2_gsa19_non)["Non_ONonO"]
SE_gsa19 <- sqrt(vcov(multivariate_model2_gsa19_non)["Non_ONonO", "Non_ONonO"])

logOR_onco <- coef(multivariate_model2_onco_non)["Non_ONonO"]
SE_onco <- sqrt(vcov(multivariate_model2_onco_non)["Non_ONonO", "Non_ONonO"])

# Prepare the data for meta-analysis
res <- data.frame(
  study = c("GSA24", "GSA19", "Onco"),
  logOR = c(logOR_gsa24, logOR_gsa19, logOR_onco),
  SE = c(SE_gsa24, SE_gsa19, SE_onco)
)

# Run the meta-analysis
meta_analysis_non_o <- rma.uni(
  yi = logOR, 
  sei = SE, 
  data = res, 
  method = "REML", 
  slab = res$study
)

# Display the meta-analysis results
summary(meta_analysis_non_o)

# Forest plot
forest(
  meta_analysis_non_o,
  xlab = "Odds Ratio",
  atransf = exp,  # Transform logOR to OR
  at = log(c(0.5, 1, 2, 5)),  # x-axis tick marks in OR scale
  xlim = c(-2, 2),            # Adjust x-axis limits as needed
  main = "Meta-Analysis: Non-O vs O"
)

# with hla
# Assuming df is your data frame, with the outcome variable named "outcome"
# The HLA allele columns start from column 7 to 107

# Initialize lists to store results
#univariate_results <- list()
#univariate_or_results <- list()
#hla_col <-names(combined_df[7:108])
# Loop over columns 7 to 107
#for (i in 7:82) {
 # hla_col <- names(combined_df)[i]  # Get the column name
  
  # Fit a logistic regression model with `TYPE` as a covariate
  #model <- glm(casecontrol ~ combined_df[[hla_col]] + TYPE, data = combined_df, family = binomial(logit))
  
  # Store the model summary or coefficients
  #univariate_results[[hla_col]] <- summary(model)
  
  # Calculate OR and confidence intervals for the model
  #or_ci <- exp(cbind(OR = coef(model), confint(model)))  # Odds Ratios with confidence intervals
  
  # Store the OR and CI in the list
  #univariate_or_results[[hla_col]] <- or_ci
#}

# Print or inspect results for each HLA column
#univariate_results
#univariate_or_results

#multivariate_results <- list()
#multivariate_or_results <- list()

#for (i in 7:108) {
#  hla_col <- names(combined_df)[i]  # Get the column name
  
  # Fit a logistic regression model with `TYPE` as a covariate
  #model <- glm(casecontrol ~ combined_df[[hla_col]] + TYPE + agec + sex + EuroRegion, data = combined_df, family = binomial(logit))
  
  # Store the model summary or coefficients
  #multivariate_results[[hla_col]] <- summary(model)
  
  # Calculate OR and confidence intervals for the model
  #or_ci <- exp(cbind(OR = coef(model), confint(model)))  # Odds Ratios with confidence intervals
  
  # Store the OR and CI in the list
  #multivariate_or_results[[hla_col]] <- or_ci
#}

# Print or inspect results for each HLA column
#multivariate_results
#multivariate_or_results

#for (i in 16:dim(combined_df)[2]) {
#  combined_df[,i] <- as.factor(combined_df[,i])
#}

# Create an empty dataframe to store results
interaction_results_df <- data.frame()

# Loop through HLA columns
hla_columns <- grep("^HLA_", names(combined_df), value = TRUE)  # Select HLA columns dynamically
for (hla_col in hla_columns) {
  
  # Skip if the column has only one level
  if (length(unique(combined_df[[hla_col]])) <= 1) {
    cat("Skipping column", hla_col, "- only one level.\n")
    next
  }
  # Build the formula dynamically
  formula <- as.formula(paste(
    "casecontrol ~ Non_O +", hla_col, "+ Non_O *", hla_col, "+ agec + sex + PC1 + PC2 + PC3 + PC4 + PC5"
  ))
  
   #Fit the logistic regression model
  model <- glm(formula, data = combined_df, family = binomial)
  
  # Get model coefficients
  coef_summary <- summary(model)$coefficients
  
  # Interaction term names
  interaction_term <- paste0('Non_ONonO:',hla_col)
  
  if (interaction_term %in% rownames(coef_summary)) {
    # Extract coefficients, calculate OR, and confidence intervals
    coef <- coef_summary[interaction_term, "Estimate"]
    se <- coef_summary[interaction_term, "Std. Error"]
    p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
    OR <- exp(coef)
    Lower_CI <- exp(coef - 1.96 * se)
    Upper_CI <- exp(coef + 1.96 * se)
    
    # Store results
    interaction_results_df <- rbind(
      interaction_results_df,
      data.frame(
        HLA_Allele = hla_col,
        Group_Level = paste0('Non-O'),
        Coef = coef,
        OR = OR,
        Lower_CI = Lower_CI,
        Upper_CI = Upper_CI,
        p.value = p_value,
        stringsAsFactors = FALSE
      )
    )
  } else {
    cat("Interaction term", interaction_term, "not found in the model.\n")
  }
}

# Display or save results
print(interaction_results_df)

hla_adjusted_pvalues <- p.adjust(interaction_results_df$p.value, method = 'BH')
interaction_results_df_adjusted <- data.frame(interaction_results_df,hla_adjusted_pvalues)
sorted_df <- interaction_results_df_adjusted[order(interaction_results_df_adjusted$p.value), ]

##############################
# Ensure results storage
interaction_results_df <- data.frame()

# Loop through HLA columns
hla_columns <- grep("^HLA_", names(combined_df), value = TRUE)

for (hla_col in hla_columns) {
  # Ensure HLA column is treated as a factor
  #combined_df[[hla_col]] <- as.factor(combined_df[[hla_col]])
  
  # Skip if the column has only one level
  if (length(unique(combined_df[[hla_col]])) <= 1) {
    cat("Skipping column", hla_col, "- only one level.\n")
    next
  }
  
  # Build the formula dynamically
  formula <- as.formula(paste(
    "casecontrol ~ BGroup +", hla_col, "+ BGroup*", hla_col, "+ agec + sex + PC1 + PC2 + PC3 + PC4 + PC5"
  ))
  
  # Fit the logistic regression model
  model <- glm(formula, data = combined_df, family = binomial)
  
  # Get model coefficients
  coef_summary <- summary(model)$coefficients
  
  # Loop through BGroup levels
  for (bgroup_level in levels(combined_df$BGroup)[-1]) {  # Exclude the reference level
    # Construct the interaction term
    interaction_term <- paste0("BGroup", bgroup_level, ":", hla_col)
    #interaction_term <- hla_col
    
    # Check if the interaction term exists in the model
    if (interaction_term %in% rownames(coef_summary)) {
      # Extract coefficients, calculate OR, and confidence intervals
      coef <- coef_summary[interaction_term, "Estimate"]
      se <- coef_summary[interaction_term, "Std. Error"]
      p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
      OR <- exp(coef)
      Lower_CI <- exp(coef - 1.96 * se)
      Upper_CI <- exp(coef + 1.96 * se)
      
      # Store results
      interaction_results_df <- rbind(
        interaction_results_df,
        data.frame(
          HLA_Allele = hla_col,
          BGroup_Level = bgroup_level,
          Coef = coef,
          OR = OR,
          Lower_CI = Lower_CI,
          Upper_CI = Upper_CI,
          p.value = p_value,
          stringsAsFactors = FALSE
        )
      )
    } else {
      cat("Interaction term", interaction_term, "not found in the model for BGroup level", bgroup_level, "and allele", hla_col, "\n")
    }
  }
}
# Display the results
print(interaction_results_df)

hla_adjusted_pvalues <- p.adjust(interaction_results_df$p.value, method = 'BH')
interaction_results_df_adjusted <- data.frame(interaction_results_df,hla_adjusted_pvalues)
sorted_df <- interaction_results_df_adjusted[order(interaction_results_df_adjusted$p.value), ]

########################################################
# Specify the range of columns to check
start_col <- 16
end_col <- 108

transformed_combined_df <- combined_df

# Replace all occurrences of 2 with 1 in the specified columns
transformed_combined_df[, start_col:end_col] <- lapply(
  transformed_combined_df[, start_col:end_col],
  function(column) ifelse(column == 2, 1, column)
)

# Ensure results storage
interaction_results_df <- data.frame()

# Loop through HLA columns
hla_columns <- grep("^HLA_", names(transformed_combined_df), value = TRUE)

for (hla_col in hla_columns) {
  # Ensure HLA column is treated as a factor
  #transformed_combined_df[[hla_col]] <- as.factor(transformed_combined_df[[hla_col]])
  
  # Skip if the column has only one level
  if (length(unique(transformed_combined_df[[hla_col]])) <= 1) {
    cat("Skipping column", hla_col, "- only one level.\n")
    next
  }
  
  # Build the formula dynamically
  formula <- as.formula(paste(
    "casecontrol ~ BGroup +", hla_col, "+ BGroup*", hla_col, "+ agec + sex + PC1 + PC2 + PC3 + PC4 + PC5"
  ))
  
  # Fit the logistic regression model
  model <- glm(formula, data = transformed_combined_df, family = binomial)
  
  # Get model coefficients
  coef_summary <- summary(model)$coefficients
  
  # Loop through BGroup levels
  for (bgroup_level in levels(transformed_combined_df$BGroup)[-1]) {  # Exclude the reference level
    # Construct the interaction term
    interaction_term <- paste0("BGroup", bgroup_level, ":", hla_col)
    #interaction_term <- hla_col
    
    # Check if the interaction term exists in the model
    if (interaction_term %in% rownames(coef_summary)) {
      # Extract coefficients, calculate OR, and confidence intervals
      coef <- coef_summary[interaction_term, "Estimate"]
      se <- coef_summary[interaction_term, "Std. Error"]
      p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
      OR <- exp(coef)
      Lower_CI <- exp(coef - 1.96 * se)
      Upper_CI <- exp(coef + 1.96 * se)
      
      # Store results
      interaction_results_df <- rbind(
        interaction_results_df,
        data.frame(
          HLA_Allele = hla_col,
          BGroup_Level = bgroup_level,
          Coef = coef,
          OR = OR,
          Lower_CI = Lower_CI,
          Upper_CI = Upper_CI,
          p.value = p_value,
          stringsAsFactors = FALSE
        )
      )
    } else {
      cat("Interaction term", interaction_term, "not found in the model for BGroup level", bgroup_level, "and allele", hla_col, "\n")
    }
  }
}
# Display the results
print(interaction_results_df)

hla_adjusted_pvalues <- p.adjust(interaction_results_df$p.value, method = 'BH')
interaction_results_df_adjusted <- data.frame(interaction_results_df,hla_adjusted_pvalues)
sorted_df <- interaction_results_df_adjusted[order(interaction_results_df_adjusted$p.value), ]
sorted_df %>%
  filter(OR > 1, p.value < 0.05)
##############################
# Create an empty dataframe to store results
interaction_results_df <- data.frame()

# Loop through HLA columns
hla_columns <- grep("^HLA_", names(transformed_combined_df), value = TRUE)  # Select HLA columns dynamically
for (hla_col in hla_columns) {
  
  # Skip if the column has only one level
  if (length(unique(transformed_combined_df[[hla_col]])) <= 1) {
    cat("Skipping column", hla_col, "- only one level.\n")
    next
  }
  # Build the formula dynamically
  formula <- as.formula(paste(
    "casecontrol ~ Non_O +", hla_col, "+ Non_O *", hla_col, "+ agec + sex + PC1 + PC2 + PC3 + PC4 + PC5"
  ))
  
  #Fit the logistic regression model
  model <- glm(formula, data = transformed_combined_df, family = binomial)
  
  # Get model coefficients
  coef_summary <- summary(model)$coefficients
  
  # Interaction term names
  interaction_term <- paste0('Non_ONonO:',hla_col)
  
  if (interaction_term %in% rownames(coef_summary)) {
    # Extract coefficients, calculate OR, and confidence intervals
    coef <- coef_summary[interaction_term, "Estimate"]
    se <- coef_summary[interaction_term, "Std. Error"]
    p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
    OR <- exp(coef)
    Lower_CI <- exp(coef - 1.96 * se)
    Upper_CI <- exp(coef + 1.96 * se)
    
    # Store results
    interaction_results_df <- rbind(
      interaction_results_df,
      data.frame(
        HLA_Allele = hla_col,
        Group_Level = paste0('Non-O'),
        Coef = coef,
        OR = OR,
        Lower_CI = Lower_CI,
        Upper_CI = Upper_CI,
        p.value = p_value,
        stringsAsFactors = FALSE
      )
    )
  } else {
    cat("Interaction term", interaction_term, "not found in the model.\n")
  }
}

# Display or save results
print(interaction_results_df)

hla_adjusted_pvalues <- p.adjust(interaction_results_df$p.value, method = 'BH')
interaction_results_df_adjusted <- data.frame(interaction_results_df,hla_adjusted_pvalues)
sorted_df <- interaction_results_df_adjusted[order(interaction_results_df_adjusted$p.value), ]
sorted_df %>%
  filter(OR > 1, p.value < 0.05)
####################################

## MICROBIOME
filtered_data_pdac <- subset(microbiome_gut_whole, casecontrol %in% c(0, 1))
filtered_data_pdac$casecontrol <- droplevels(filtered_data_pdac$casecontrol)
filtered_data_pdac <- drop_na(filtered_data_pdac)

### Filter out rows where 'sex' is 2
filtered_data_pdac <- filtered_data_pdac %>% 
  filter(sex != 2) # %>%

# Relevel for default
filtered_data_pdac$casecontrol <- relevel(filtered_data_pdac$casecontrol, ref = '0')
filtered_data_pdac$Non_O <- relevel(filtered_data_pdac$Non_O, ref = "O")
filtered_data_pdac$BGroup <- droplevels(filtered_data_pdac$BGroup)
filtered_data_pdac$BGroup <- relevel(filtered_data_pdac$BGroup, ref = "O")
filtered_data_pdac$casecontrol <- droplevels(filtered_data_pdac$casecontrol)

median(filtered_data_pdac$q.0,na.rm=TRUE) # 75.27
median(filtered_data_pdac$q.1,na.rm=TRUE) # 60.11086

filtered_data_pdac <- filtered_data_pdac %>%
  mutate(
    q.0_b = case_when(
      q.0 > 75.27 ~ "high",
      TRUE ~ "low"
    )
  )
filtered_data_pdac$q.0_b <- as.factor(filtered_data_pdac$q.0_b)
filtered_data_pdac$q.0_b <- relevel(filtered_data_pdac$q.0_b, ref = "low")

filtered_data_pdac <- filtered_data_pdac %>%
  mutate(
    q.1_b = case_when(
      q.1 > 60.11086 ~ "high",
      TRUE ~ "low"
    )
  )
filtered_data_pdac$q.1_b <- as.factor(filtered_data_pdac$q.1_b)
filtered_data_pdac$q.1_b <- relevel(filtered_data_pdac$q.1_b, ref = "low")

# Boxplot for Q0
ggplot(filtered_data_pdac, aes(x = Non_O, y = q.0, fill = Non_O)) +
  geom_boxplot() +
  labs(title = "Q0 Distribution by Non-O and O", x = "Group", y = "Q0") +
  theme_minimal()

# Boxplot for Q1
ggplot(filtered_data_pdac, aes(x = Non_O, y = q.1, fill = Non_O)) +
  geom_boxplot() +
  labs(title = "Q1 Distribution by Non-O and O", x = "Group", y = "Q1") +
  theme_minimal()

##### NonO and Microbiome Q0
# Define the formula
formula <- as.formula(paste(
  "casecontrol ~ Non_O + q.0 + (Non_O * q.0) + agec + sex + PC1 + PC2 + PC3 + PC4 + PC5"
))

# Fit the logistic regression model
model <- glm(formula, data = filtered_data_pdac, family = binomial)

# Get model coefficients
coef_summary <- summary(model)$coefficients

# Initialize the results data frame
interaction_results_df <- data.frame()

# Extract the coefficients for the interaction term
interaction_term <- "Non_ONonO:q.0"
coef <- coef_summary[interaction_term, "Estimate"]
se <- coef_summary[interaction_term, "Std. Error"]
p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
OR <- exp(coef)
Lower_CI <- exp(coef - 1.96 * se)
Upper_CI <- exp(coef + 1.96 * se)
  
  # Store results in the data frame
interaction_results_df <- rbind(
    interaction_results_df,
    data.frame(
      Group_Level = "Non-O",
      Coef = coef,
      OR = OR,
      Lower_CI = Lower_CI,
      Upper_CI = Upper_CI,
      p.value = p_value,
      stringsAsFactors = FALSE
    )
  )

# Adjust the p-value (unnecessary for a single p-value, but kept for consistency)
micro_adjusted_pvalues <- p.adjust(p_value, method = 'BH')

# Add the adjusted p-value to the results data frame
interaction_results_df$micro_adjusted_pvalues <- micro_adjusted_pvalues

# Display the results
print(interaction_results_df)

##### ABO and Microbiome Q0
# Define the formula
formula <- as.formula(paste(
  "casecontrol ~ BGroup + q.0 + (BGroup * q.0) + agec + sex + PC1 + PC2 + PC3 + PC4 + PC5"
))

# Fit the logistic regression model
model <- glm(formula, data = filtered_data_pdac, family = binomial)

# Get model coefficients
coef_summary <- summary(model)$coefficients

# Initialize the results data frame
interaction_results_df <- data.frame()

# Define BGroup levels of interest (excluding the reference level, "O")
bgroup_levels <- levels(filtered_data_pdac$BGroup)[-1]  # Exclude reference level "O"

# Loop through BGroup levels and extract interaction terms
for (level in bgroup_levels) {
  interaction_term <- paste0("BGroup", level, ":q.0")  # Construct interaction term
  
  if (interaction_term %in% rownames(coef_summary)) {
    coef <- coef_summary[interaction_term, "Estimate"]
    se <- coef_summary[interaction_term, "Std. Error"]
    p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
    OR <- exp(coef)
    Lower_CI <- exp(coef - 1.96 * se)
    Upper_CI <- exp(coef + 1.96 * se)
    
    # Store results in the data frame
    interaction_results_df <- rbind(
      interaction_results_df,
      data.frame(
        BGroup_Level = level,
        Coef = coef,
        OR = OR,
        Lower_CI = Lower_CI,
        Upper_CI = Upper_CI,
        p.value = p_value,
        stringsAsFactors = FALSE
      )
    )
  } else {
    cat("Interaction term", interaction_term, "not found in the model.\n")
  }
}

# Adjust the p-values for all interaction terms
interaction_results_df$adjusted_pvalue <- p.adjust(interaction_results_df$p.value, method = "BH")

# Display the results
print(interaction_results_df)

##### NonO and Microbiome Q1
# Define the formula
formula <- as.formula(paste(
  "casecontrol ~ Non_O + q.1 + (Non_O * q.1) + agec + sex + PC1 + PC2 + PC3 + PC4 + PC5"
))

# Fit the logistic regression model
model <- glm(formula, data = filtered_data_pdac, family = binomial)

# Get model coefficients
coef_summary <- summary(model)$coefficients

# Initialize the results data frame
interaction_results_df <- data.frame()

# Extract the coefficients for the interaction term
interaction_term <- "Non_ONonO:q.1"
coef <- coef_summary[interaction_term, "Estimate"]
se <- coef_summary[interaction_term, "Std. Error"]
p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
OR <- exp(coef)
Lower_CI <- exp(coef - 1.96 * se)
Upper_CI <- exp(coef + 1.96 * se)

# Store results in the data frame
interaction_results_df <- rbind(
  interaction_results_df,
  data.frame(
    Group_Level = "Non-O",
    Coef = coef,
    OR = OR,
    Lower_CI = Lower_CI,
    Upper_CI = Upper_CI,
    p.value = p_value,
    stringsAsFactors = FALSE
  )
)

# Adjust the p-value (unnecessary for a single p-value, but kept for consistency)
micro_adjusted_pvalues <- p.adjust(p_value, method = 'BH')

# Add the adjusted p-value to the results data frame
interaction_results_df$micro_adjusted_pvalues <- micro_adjusted_pvalues

# Display the results
print(interaction_results_df)

##### ABO and Microbiome Q1
# Define the formula
formula <- as.formula(paste(
  "casecontrol ~ BGroup + q.1 + (BGroup * q.1) + agec + sex + PC1 + PC2 + PC3 + PC4 + PC5"
))

# Fit the logistic regression model
model <- glm(formula, data = filtered_data_pdac, family = binomial)

# Get model coefficients
coef_summary <- summary(model)$coefficients

# Initialize the results data frame
interaction_results_df <- data.frame()

# Define BGroup levels of interest (excluding the reference level, "O")
bgroup_levels <- levels(filtered_data_pdac$BGroup)[-1]  # Exclude reference level "O"

# Loop through BGroup levels and extract interaction terms
for (level in bgroup_levels) {
  interaction_term <- paste0("BGroup", level, ":q.1")  # Construct interaction term
  
  if (interaction_term %in% rownames(coef_summary)) {
    coef <- coef_summary[interaction_term, "Estimate"]
    se <- coef_summary[interaction_term, "Std. Error"]
    p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
    OR <- exp(coef)
    Lower_CI <- exp(coef - 1.96 * se)
    Upper_CI <- exp(coef + 1.96 * se)
    
    # Store results in the data frame
    interaction_results_df <- rbind(
      interaction_results_df,
      data.frame(
        BGroup_Level = level,
        Coef = coef,
        OR = OR,
        Lower_CI = Lower_CI,
        Upper_CI = Upper_CI,
        p.value = p_value,
        stringsAsFactors = FALSE
      )
    )
  } else {
    cat("Interaction term", interaction_term, "not found in the model.\n")
  }
}

# Adjust the p-values for all interaction terms
interaction_results_df$adjusted_pvalue <- p.adjust(interaction_results_df$p.value, method = "BH")

# Display the results
print(interaction_results_df)

## Stratified analysis with continuous values
# Stratify the data into Non-O and O groups
#non_o_data <- filtered_data_pdac %>% filter(Non_O == "NonO")
#o_data <- filtered_data_pdac %>% filter(Non_O == "O")

# Fit models for Non-O and O
model_q0 <- glm(log(q.0) ~ casecontrol+ Non_O + agec + sex, 
                   data = filtered_data_pdac)
model_q1 <- glm(log(q.1) ~ casecontrol+ Non_O + agec + sex, 
                   data = filtered_data_pdac)

# Summarize results
summary(model_q0)
summary(model_q1)

## Plot
# Define each regression table with specific data subsets and labels
table1 <- tbl_regression(model_q0, exponentiate = FALSE)
table2 <- tbl_regression(model_q1, exponentiate = FALSE) 

table1 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/micro/model_0.pdf")

table2 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/micro/model_1.pdf")

forest_plot1 <- 
  table1 %>%
  modify_column_merge(
    pattern = "{estimate}; p-value: {p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 1**"))

forest_plot1 |> 
  gt::gtsave(filename = "~/Downloads/tfm/micro/model_0_forest.pdf")

forest_plot2 <- 
  table2 %>%
  modify_column_merge(
    pattern = "{estimate}; p-value: {p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 2**"))

forest_plot2 |> 
  gt::gtsave(filename = "~/Downloads/tfm/micro/model_1_forest.pdf")


# Merge the tables
# To hide repeated demographic labels, we can leverage gtsummary's modify_* functions
tbl_merge2 <- tbl_merge(
  tbls = list(table1, table2),
  tab_spanner = c("Q0", "Q1") 
)%>%
  modify_spanning_header(everything() ~ NA_character_) %>%  # Remove existing spanning headers if needed
  modify_column_unhide(columns = c("label", "estimate_1", "ci_1", 
                                   "estimate_2", "ci_2"))
# Convert to gt format if needed
gt_table <- as_gt(tbl_merge2)
gtsave(gt_table, filename = "~/Downloads/tfm/micro/microq_tables.html")

# Perform the t-tests
q0_test <- t.test(filtered_data_pdac$q.0, filtered_data_pdac$q.0, na.rm = TRUE, var.equal = FALSE)  # Q0 comparison
q1_test <- t.test(filtered_data_pdac$q.1, filtered_data_pdac$q.1, na.rm = TRUE, var.equal = FALSE)  # Q1 comparison

# Extract results into a data frame
t_test_results <- data.frame(
  Test = c("Q0 Comparison", "Q1 Comparison"),
  Mean_Group1 = c(q0_test$estimate[1], q1_test$estimate[1]),
  Mean_Group2 = c(q0_test$estimate[2], q1_test$estimate[2]),
  t_value = c(q0_test$statistic, q1_test$statistic),
  p_value = c(q0_test$p.value, q1_test$p.value),
  Confidence_Interval_Lower = c(q0_test$conf.int[1], q1_test$conf.int[1]),
  Confidence_Interval_Upper = c(q0_test$conf.int[2], q1_test$conf.int[2])
)

# Display the results table
print(t_test_results)

# Create a well-designed table with gt
t_test_results_gt <- t_test_results %>%
  gt() %>%
  tab_header(
    title = "T-Test Results for Q0 and Q1",
    subtitle = "Comparison between Non-O and O groups"
  ) %>%
  fmt_number(
    columns = c(Mean_Group1, Mean_Group2, t_value, p_value, Confidence_Interval_Lower, Confidence_Interval_Upper),
    decimals = 3
  ) %>%
  cols_label(
    Test = "Test",
    Mean_Group1 = "Mean (Group 1)",
    Mean_Group2 = "Mean (Group 2)",
    t_value = "T-Value",
    p_value = "P-Value",
    Confidence_Interval_Lower = "Lower CI",
    Confidence_Interval_Upper = "Upper CI"
  )

# Save the table as an HTML or PDF file
gtsave(t_test_results_gt, filename = "~/Downloads/tfm/micro/t_test_results.html")  # Save as HTML
gtsave(t_test_results_gt, filename = "~/Downloads/tfm/micro/t_test_results.pdf")   # Save as PDF

# Filter the data to include only controls (casecontrol == 0)
control_data <- filtered_data_pdac %>% filter(casecontrol == 0)

# Stratify the controls into Non-O and O groups
#non_o_controls <- control_data %>% filter(Non_O == "NonO")
#o_controls <- control_data %>% filter(Non_O == "O")

# Fit models for Non-O and O controls
control_model_0 <- glm(log(q.0) ~ Non_O + agec + sex , 
                           data = control_data)
control_model_1 <- glm(log(q.1) ~ Non_O + agec + sex, 
                           data = control_data)

# Summarize results
summary(control_model_0)
summary(control_model_1)

## Plot
# Define each regression table with specific data subsets and labels
table1 <- tbl_regression(control_model_0, exponentiate = FALSE)
table2 <- tbl_regression(control_model_1, exponentiate = FALSE) 

table1 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/micro/control_model_0.pdf")

table2 |> 
  as_gt() |> 
  gt::gtsave(filename = "~/Downloads/tfm/micro/control_model_1.pdf")

forest_plot1 <- 
  table1 %>%
  modify_column_merge(
    pattern = "{estimate}; p-value: {p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 1**"))

forest_plot1 |> 
  gt::gtsave(filename = "~/Downloads/tfm/micro/control_model_0_forest.pdf")

forest_plot2 <- 
  table2 %>%
  modify_column_merge(
    pattern = "{estimate}; p-value: {p.value}",
    rows = !is.na(estimate)
  ) %>%
  add_forest() %>%
  gt::tab_header(title = gt::md("**Forest Plot of Model Multi 2**"))

forest_plot2 |> 
  gt::gtsave(filename = "~/Downloads/tfm/micro/control_model_1_forest.pdf")

# Merge the tables
# To hide repeated demographic labels, we can leverage gtsummary's modify_* functions
tbl_merge2 <- tbl_merge(
  tbls = list(table1, table2),
  tab_spanner = c("Q0", "Q1") 
)%>%
  modify_spanning_header(everything() ~ NA_character_) %>%  # Remove existing spanning headers if needed
  modify_column_unhide(columns = c("label", "estimate_1", "ci_1", 
                                   "estimate_2", "ci_2"))
# Convert to gt format if needed
gt_table <- as_gt(tbl_merge2)
gtsave(gt_table, filename = "~/Downloads/tfm/micro/controls_micro_tables.html")

# Compare Q0/Q1 between Non-O and O controls
# Perform the t-tests
q0_control_test <- t.test(control_data$q.0, control_data$q.0, na.rm = TRUE, var.equal = FALSE)  # Q0 comparison
q1_control_test <- t.test(control_data$q.1, control_data$q.1, na.rm = TRUE, var.equal = FALSE)  # Q1 comparison

t_test_results <- data.frame(
  Test = c("Q0 Comparison", "Q1 Comparison"),
  Mean_Group1 = c(q0_control_test$estimate[1], q1_control_test$estimate[1]),
  Mean_Group2 = c(q0_control_test$estimate[2], q1_control_test$estimate[2]),
  t_value = c(q0_control_test$statistic, q1_control_test$statistic),
  p_value = c(q0_control_test$p.value, q1_control_test$p.value),
  Confidence_Interval_Lower = c(q0_control_test$conf.int[1], q1_control_test$conf.int[1]),
  Confidence_Interval_Upper = c(q0_control_test$conf.int[2], q1_control_test$conf.int[2])
)

# Display the results table
print(t_test_results)

# Create a well-designed table with gt
t_test_results_gt <- t_test_results %>%
  gt() %>%
  tab_header(
    title = "T-Test Results for Q0 and Q1",
    subtitle = "Comparison between Non-O and O groups"
  ) %>%
  fmt_number(
    columns = c(Mean_Group1, Mean_Group2, t_value, p_value, Confidence_Interval_Lower, Confidence_Interval_Upper),
    decimals = 3
  ) %>%
  cols_label(
    Test = "Test",
    Mean_Group1 = "Mean (Group 1)",
    Mean_Group2 = "Mean (Group 2)",
    t_value = "T-Value",
    p_value = "P-Value",
    Confidence_Interval_Lower = "Lower CI",
    Confidence_Interval_Upper = "Upper CI"
  )

# Save the table as an HTML or PDF file
gtsave(t_test_results_gt, filename = "~/Downloads/tfm/micro/t_test_controls_results.html")  # Save as HTML
gtsave(t_test_results_gt, filename = "~/Downloads/tfm/micro/t_test_controls_results.pdf")   # Save as PDF

# Create a well-designed table with flextable
t_test_results_flextable <- t_test_results %>%
  flextable() %>%
  set_header_labels(
    Test = "Test",
    Mean_Group1 = "Mean (Group 1)",
    Mean_Group2 = "Mean (Group 2)",
    t_value = "T-Value",
    p_value = "P-Value",
    Confidence_Interval_Lower = "Lower CI",
    Confidence_Interval_Upper = "Upper CI"
  ) %>%
  add_header_row(
    values = c("T-Test Results for Q0 and Q1"), colwidths = 7
  ) %>%
  theme_vanilla() %>%
  autofit()

# Save the table as a Word or HTML file
save_as_docx(t_test_results_flextable, path = "~/Downloads/tfm/micro/flex_t_test_results.docx")  # Save as Word
save_as_html(t_test_results_flextable, path = "~/Downloads/tfm/micro/flex_t_test_results.html")  # Save as HTML

all_data_test_q0<-t.test(filtered_data_pdac$q.0, filtered_data_pdac$q.0, na.rm = TRUE) # All data
controls_test_q0<-t.test(control_data$q.0, control_data$q.0, na.rm = TRUE) # Controls only

all_data_test_q1<-t.test(filtered_data_pdac$q.1, filtered_data_pdac$q.1, na.rm = TRUE) # All data
controls_test_q1<-t.test(control_data$q.1, control_data$q.1, na.rm = TRUE) # Controls only

non_o_data <- filtered_data_pdac %>% filter(Non_O == "NonO")
o_data <- filtered_data_pdac %>% filter(Non_O == "O")

# Combine data into one data frame
all_data <- rbind(
  data.frame(Group = "Non-O", Q = non_o_data$q.0, Variable = "Q0"),
  data.frame(Group = "O", Q = o_data$q.0, Variable = "Q0"),
  data.frame(Group = "Non-O", Q = non_o_data$q.1, Variable = "Q1"),
  data.frame(Group = "O", Q = o_data$q.1, Variable = "Q1")
)

ggplot(all_data, aes(x = Group, y = Q, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(
    title = "Comparison of Q0 and Q1 between Non-O and O Groups",
    x = "Group",
    y = "Q Values"
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Non-O" = "#1f77b4", "O" = "#ff7f0e"))
ggsave("~/Downloads/tfm/micro/q_comparison_barplot.png", width = 8, height = 5)


summary_data <- all_data %>%
  group_by(Variable, Group) %>%
  summarise(
    Mean = mean(Q, na.rm = TRUE),
    SD = sd(Q, na.rm = TRUE),
    n = sum(!is.na(Q)),
    SE = SD / sqrt(n)
  ) %>%
  mutate(
    CI_Lower = Mean - 1.96 * SE,
    CI_Upper = Mean + 1.96 * SE
  )

# Data preparation
bar_data <- data.frame(
  Group = rep(c("Non-O", "O"), 4),
  Test = rep(c("q.0 (All)", "q.0 (Controls)", "q.1 (All)", "q.1 (Controls)"), each = 2),
  Mean = c(
    mean(non_o_data$q.0, na.rm = TRUE), mean(o_data$q.0, na.rm = TRUE),
    mean(non_o_controls$q.0, na.rm = TRUE), mean(o_controls$q.0, na.rm = TRUE),
    mean(non_o_data$q.1, na.rm = TRUE), mean(o_data$q.1, na.rm = TRUE),
    mean(non_o_controls$q.1, na.rm = TRUE), mean(o_controls$q.1, na.rm = TRUE)
  ),
  P_Value = c(
    all_data_test_q0$p.value, all_data_test_q0$p.value,
    controls_test_q0$p.value, controls_test_q0$p.value,
    all_data_test_q1$p.value, all_data_test_q1$p.value,
    controls_test_q1$p.value, controls_test_q1$p.value
  )
)

library(ggpubr)

# Barplot
# Load ggpubr library
library(ggpubr)

# Create the bar plot
p <- ggbarplot(
  bar_data,
  x = "Test",  # x-axis
  y = "Mean",  # y-axis
  fill = "Group",  # Grouping variable
  position = position_dodge(0.7),  # Dodge position for bars
  width = 0.7,  # Bar width
  palette = c("skyblue", "orange"),  # Custom colors
  add = "mean",  # Add the mean values
  ylab = "Mean Value",
  xlab = "Welch's Test",
  title = "Comparison of Means for q.0 and q.1 by Non-O vs O Group"
) +
  stat_compare_means(
    aes(group = Group),  # Grouping for comparisons
    label = "p.signif",  # Use significant symbols (e.g., *, **, etc.)
    method = "t.test",  # Specify the test method
    label.y = max(bar_data$Mean) + 5,  # Position of p-value labels above bars
    bracket.size = 0.5  # Size of comparison brackets
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

# Display the plot
print(p)

ggsave("~/Downloads/tfm/micro/q_comparison_barplot_test.png", width = 8, height = 5)

## ALBERTOS CODE:
for (i in 7:82) {
  hla_col <- names(combined_df)[i]  # Get the column name
  
  formula <- as.formula(paste(
      "Class1_LOSS ~", allele_valid, "* G12V + Subtype...Moffitt + TP53 + Lymphnode_stage_collapsed + Radiation.Therapy"
    )) 
    ASE_interaction <- glm(PDACcaco ~ ABO + HLA[i] + ABO*HLA[i] + covars , data = PDAC_data_final, family=binomial(logit))
  summary(ASE_interaction)

}

# Ensure G12V is a factor with the correct levels
PDAC_data_final$G12V <- factor(
  PDAC_data_final$G12V,
  levels = c("Other_mutation", "WT", "G12V")
)
# Initialize results data frame
interaction_results_df <- data.frame()
# Loop over each allele
for (i in seq_along(alleles)) {
  allele <- alleles[i]
  allele_valid <- alleles_valid[i]
  # Proceed only if allele variable has more than one level
  if (length(unique(PDAC_data_final[[allele_valid]])) > 1) {
    # Build the formula
    formula <- as.formula(paste(
      "Class1_LOSS ~", allele_valid, "* G12V + allele_valed+ G12V+ Subtype...Moffitt + TP53 + Lymphnode_stage_collapsed + Radiation.Therapy"
    ))
    # Fit the model
    model <- glm(formula, data = PDAC_data_final, family = binomial)
    # Extract coefficients
    coef_summary <- summary(model)$coefficients
    # Print coefficient names for debugging
    print(paste("Coefficients for allele:", allele))
    print(rownames(coef_summary))
    # G12V levels (excluding reference level)
    G12V_levels <- levels(PDAC_data_final$G12V)[-1]
    # Loop over G12V levels
    for (G12V_level in G12V_levels) {
      # Construct the interaction term name
      interaction_term <- paste0(allele_valid, "1:G12V", G12V_level)
      if (interaction_term %in% rownames(coef_summary)) {
        coef <- coef_summary[interaction_term, "Estimate"]
        se <- coef_summary[interaction_term, "Std. Error"]
        p_value <- coef_summary[interaction_term, "Pr(>|z|)"]
        # Calculate Odds Ratio and Confidence Intervals
        OR <- exp(coef)
        Lower_CI <- exp(coef - 1.96 * se)
        Upper_CI <- exp(coef + 1.96 * se)
        # Store the results
        interaction_results_df <- rbind(
          interaction_results_df,
          data.frame(
            Allele = allele,
            Allele_Level = "1",
            G12V_Level = G12V_level,
            Interaction_Term = interaction_term,
            Coef = coef,
            OR = OR,
            Lower_CI = Lower_CI,
            Upper_CI = Upper_CI,
            p.value = p_value,
            stringsAsFactors = FALSE
          )
        )
      } else {
        cat("Interaction term", interaction_term, "not found in the model for allele", allele, "\n")
      }
    }
  } else {
    cat("Allele", allele, "has only one level and will be skipped.\n")
  }
}
###################

#for (hla_name in colnames(hla_com_filtered[2:dim(hla_com_filtered)[2]])) {
 # print(hla_name)
  #print(hla_com_filtered[,hla_name])
  
#}
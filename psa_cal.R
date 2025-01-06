# Load required libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(readxl)

# Load the data
gsa_2019 <- read.csv("/local/tsaid/sharing/GSA2019_dosagesPCA.txt", sep='\t')
gsa_2024 <- read.csv("/local/tsaid/sharing/GSA2024_dosagesPCA.txt", sep='\t')

onco_is_pan <- read.csv("/local/tsaid/sharing/PED_TO_MAP/oncoarray_combined.txt",sep=',')

# Create a function to extract the genotype and convert it to numeric code
convert_genotype <- function(value) {
  # Split the value on the ":" and keep only the genotype part
  genotype <- strsplit(value, ":")[[1]][1]
  
  # Convert the genotype to numeric code
  if (genotype == "0|0") {
    return(0)
  } else if (genotype %in% c("1|0", "0|1")) {
    return(1)
  } else if (genotype == "1|1") {
    return(2)
  } else {
    return(NA)  # Return NA if genotype is missing or unknown
  }
}

# Apply the function to the entire data frame
# Assuming columns from the 10th column onward contain the genotype data
genotype_data19 <- gsa_2019[, 10:ncol(gsa_2019)]
genotype_data24 <- gsa_2024[, 10:ncol(gsa_2024)]
genotype_data19 <- apply(genotype_data19, c(1, 2), convert_genotype)
genotype_data24 <- apply(genotype_data24, c(1, 2), convert_genotype)

gsa_2019_converted <- cbind(gsa_2019[, 1:2, drop = FALSE], genotype_data19)
colnames(gsa_2019_converted) <- gsub("^X0_", "GSA19_", colnames(gsa_2019_converted))
gsa_2024_converted <- cbind(gsa_2024[, 1:2, drop = FALSE], genotype_data24)
colnames(gsa_2024_converted) <- gsub("^X0_", "GSA24_", colnames(gsa_2024_converted))

# List to remove pairs before pca 
#gsa24_list <- read.csv("/local/tsaid/trabajo/principle_comp/filter_ibd_gsa.txt", sep='\t')
gsa24_list <- c("GSA24_208025880084_R10C02", 'GSA24_208026130070_R06C01',
                'GSA24_208025880084_R10C01', 'GSA24_208025880090_R09C01',
                'GSA24_208025880084_R09C01','GSA24_208025880086_R10C02',
                'GSA24_208025880116_R06C01','GSA24_208026130093_R06C02',
                'GSA24_207769820036_R08C02','GSA24_208025880114_R02C02',
                'GSA24_208025880090_R02C02','GSA24_208026130069_R12C02',
                'GSA24_208025880110_R02C01','GSA24_208026130093_R09C01',
                'GSA24_207769820036_R08C01','GSA24_208025880090_R07C02',
                'GSA24_208026130070_R04C01','GSA24_208026130070_R09C01',
                'GSA24_208025880090_R05C01','GSA24_208025880116_R04C02')

gsa_2024_converted <- gsa_2024_converted %>% select(-all_of(gsa24_list))

# Merge the dataframes sequentially based on CHROM and POS
merged_df <- merge(gsa_2024_converted, gsa_2019_converted, by = c("X.CHROM", "POS"))
names(merged_df)[names(merged_df) == "X.CHROM"] <- "CHROM"

## Calculate the row sums and filter out rows with a sum of 0 and merge
merge_all <- merge(merged_df, onco_is_pan, by = c("CHROM", "POS"))

# Drop unnecessary columns, need only the genotypes
onco_is_pan_geno <- onco_is_pan[,3:ncol(onco_is_pan)] %>% mutate_all(as.integer)
onco_is_pan_geno <- onco_is_pan_geno[rowSums(onco_is_pan) != 0, ]

merged_df_geno19 <- gsa_2019_converted[,3:ncol(gsa_2019_converted)]
merged_df_geno19 <- merged_df_geno19[rowSums(merged_df_geno19) != 0, ]

merged_df_geno24 <- gsa_2024_converted[,3:ncol(gsa_2024_converted)]
merged_df_geno24 <- merged_df_geno24[rowSums(merged_df_geno24) != 0, ]

merge_all_geno <- merge_all[,3:ncol(merge_all)]

snp_data_filtered3 <- merge_all_geno[rowSums(merge_all_geno) != 0, ]

# Calculate PCA for each data set
pca_gsa19_trans <- prcomp(t(merged_df_geno19), scale = FALSE)
pca_gsa24_trans <- prcomp(t(merged_df_geno24), scale = FALSE)

pca_onco_pan_trans <- prcomp(t(onco_is_pan_geno), scale = TRUE)
pca_all_trans <- prcomp(t(snp_data_filtered3), scale = TRUE)

# Save the PCA result to an .RData file
save(pca_gsa19_trans, file = "/local/tsaid/trabajo/principle_comp/pca_gsa19.RData")
save(pca_gsa24_trans, file = "/local/tsaid/trabajo/principle_comp/pca_gsa24.RData")
save(pca_onco_pan_trans, file = "/local/tsaid/trabajo/principle_comp/pca_onco_pan.RData")
save(pca_all_trans, file = "/local/tsaid/trabajo/principle_comp/pca_all_trans.RData")

pca_gsa19_trans$x[, 1][order(pca_gsa19_trans$x[, 1])]
pca_gsa19_trans$x[, 2][order(pca_gsa19_trans$x[, 2])]

# View summary of PCA results
#summary(pca_gsa)
#summary(pca_onco_pan)
summary(pca_all_trans)

# Reverse the signs of PC
#pca_all_rot<- pca_all
##pca_all_rot$x <- -1*pca_all_rot$x

# Extract the principal component scores
#pca_gsa_scores <- as.data.frame(pca_gsa_rot$x)
#pca_onco_pan_scores <- as.data.frame(pca_onco_pan_rot$x)
pca_all_scores <- as.data.frame(pca_all_trans$x)

# Display the first few rows of the PCA scores
head(pca_all_scores)

# Plot 
varPC1 = round(100*pca_gsa24_trans$sdev[1]^2/sum(pca_gsa24_trans$sdev^2),digits=2)
varPC2 = round(100*pca_gsa24_trans$sdev[2]^2/sum(pca_gsa24_trans$sdev^2),digits=2)
varPC3 = round(100*pca_gsa24_trans$sdev[3]^2/sum(pca_gsa24_trans$sdev^2),digits=2)


pdf('/local/tsaid/trabajo/principle_comp/pca_plot_bad2.pdf')
# you can also get the proportion of variance of the 3 first principal components like this
100*(summary(pca_gsa24_trans)$importance["Proportion of Variance",1:2])
barplot(100*(summary(pca_gsa24_trans)$importance["Proportion of Variance",1:10]),ylab="Proportion of Variance", 
main="% Var of the 10 first principal components",)

# two first principal components
#colors_plot_array_pop <- c(rep("red",333), rep("green",151), rep("blue",153),rep("yellow",1864))# 333+151+153+1864=2499

#colors_plot_array <- c(rep("red",333), rep("green",151), rep("blue",2017))

colors_plot_array <- rep("blue",315)

plot(pca_gsa24_trans$x[,1], pca_gsa24_trans$x[,2], col=colors_plot_array, 
  main="PCA of 2024 Populations",
  xlab=paste0("P1 (",varPC1," %)"), ylab=paste0("PC2 (",varPC2," %)"),
  pch=20, col.axis="grey", fg="grey")
  #text(pca_gsa24_trans$x[,1], pca_gsa24_trans$x[,2], labels=rownames(pca_gsa24_trans$x),offset=1,cex=0.5 )
legend("top", legend="GSA24",
  text.col="blue", bg="transparent", cex=0.8)
dev.off()

plot(pca_gsa24_trans$x[,1], pca_gsa24_trans$x[,2], col=colors_plot_array, 
  main="PCA of 3 Arrays",
  xlab=paste0("PC1 (",varPC1," %)"), ylab=paste0("PC2 (",varPC2," %)"),
  pch=20, col.axis="grey", fg="grey")
  
legend("top", legend="GSA19",
  text.col="green", bg="transparent", cex=0.8)

plot(pca_gsa24_trans$x[,1], pca_gsa24_trans$x[,3], col=colors_plot_array, 
  main="PCA of 3 Arrays",
  xlab=paste0("PC1 (",varPC1," %)"), ylab=paste0("PC3 (",varPC3," %)"),
  pch=20, col.axis="grey", fg="grey")
  
legend("top", legend="GSA19",
  text.col="green", bg="transparent", cex=0.8)

plot(pca_gsa24_trans$x[,2], pca_gsa24_trans$x[,3], col=colors_plot_array, 
  main="PCA of 3 Arrays",
  xlab=paste0("PC2 (",varPC2," %)"), ylab=paste0("PC3 (",varPC3," %)"),
  pch=20, col.axis="grey", fg="grey")
  
legend("top", legend="GSA19",
  text.col="green", bg="transparent", cex=0.8)

dev.off()

#biplot(pca_all_trans, main="", xlab="", ylab=paste0("PC2 (",varPC2," %)"), cex=0.9, arrow.len=0.05, xlabs=rep(".",33723))
biplot(pca_gsa19_trans, main="", xlab="", ylab=paste0("PC2 (",varPC2," %)"), cex=0.9, arrow.len=0.05, 
  xlabs=rep(".", dim(t(merged_df_geno19))[1]), expand=1, scale=0.5, col.axis="grey", fg="grey")
legend("topleft", legend=c("Samples","SNPs"), text.col=c("black","red"), bg="transparent", cex=0.8)



sample_ids_gsa <- read.csv('/local/tsaid/trabajo/infomation/edited_IDs_arrays_2019and24_tomerge_Maria.txt', sep=',')

# Extract the prefix (e.g., GSA24) and the part after the first underscore
prefix <- gsub("_.*", "", colnames(merged_df))  # Extract the part before the first underscore
suffix <- gsub("^[^_]*_", "", colnames(merged_df))  # Extract the part after the last underscore

# Combine prefix and suffix to recreate the original column names
colname_map <- setNames(colnames(merged_df), suffix)

# Update the Name column in sample_ids_gsa
sample_ids_gsa <- sample_ids_gsa %>%
  mutate(
    Name = ifelse(
      Name %in% names(colname_map), # Check if Name exists in the suffix values
      colname_map[Name],           # Replace with the full column name (e.g., GSA24_208026130141_R03C01)
      Name                         # Keep the original value if no match is found
    )
  )
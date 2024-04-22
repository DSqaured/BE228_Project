#load libraries
library(dplyr)
library(tidyverse)
library (GEOquery)
library(DESeq2)
library(ggplot2)

#read the data
data <- read.csv("data/GSE263588_Expression_Profile.GRCh38.gene.txt", sep ="\t")

#get metadata using GEOquery
gse <- getGEO(GEO = 'GSE263588', GSEMatrix = TRUE)
#retrieving metadata
metadata <- pData(phenoData(gse[[1]]))
#grouping only counts data
data_counts <- data %>%
  select(3,11:59)             #only gene symbol and TPM values are taken from original dataset

#retreiving column information
columnData_counts <- metadata[c(1:2, 13:14)]
C_names1 <- rownames(columnData_counts)
colnames(data_counts)[2:50]<- C_names1              #identifying the corresponding gene accession IDs of samples


# Replace "Patient_derived GBM cell line" with "GBM"
columnData_counts$title <- gsub("Patient-derived GBM cell line-\\d+", "GBM", columnData_counts$title)

# Replace "Patient_derived normal tissue" with "Normal"
columnData_counts$title <- gsub("Patient-derived normal tissue-\\d+", "Normal", columnData_counts$title)

colData <- columnData_counts[1]
countData <- data_counts[2:50]
#Creating a DESeq2 object 
dds1 <- DESeqDataSetFromMatrix(countData = countData,
                               colData = colData,
                               design = ~ title)
#Performing differential gene expression analysis
keep1 <- rowSums (counts(dds1)) >= 10 #keeping rows that have at least 10 reads total
dds1 <- dds1[keep1,]

#Setting the factor level
dds1$title <- relevel(dds1$title, ref = "Normal")
dds1$title
dds1 <- DESeq(dds1)

#vsdata (how good the data is)
#vsdata1 <- vst(dds1, blind = FALSE)
#plotPCA(vsdata1, intgroup = "title")
#plotDispEsts(dds1)

#Getting differential expression results
res1 <- results(dds1)
res_df1 <- as.data.frame(res1)
#res

#Filtering results for significantly DEGs
res_sig1 <- res1[which(res1$padj <= 0.05 & abs (res1$log2FoldChange)>=1), ]
res_sig_df1 <- as.data.frame(res_sig1)
#r_t_g_1 <- rownames(res_sig_df1)
#res_sig

#viewing top DEG with the type of regulation
top_genes_1 <- head(res_sig_df1[order(res_sig_df1$log2FoldChange,
                                      decreasing = TRUE), ], 1000) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up","down"))

r_t_g1<- rownames(top_genes_1)#extracting decreasing order of top genes

# Creating a new vector to store gene symbols corresponding to gene IDs
gene_symbols <- character(length(r_t_g1))

# Looping over each gene ID in r_t_g1
for (i in seq_along(r_t_g1)) {
  # Finding the index of the gene ID in column 1 of the "data" DataFrame
  index <- which(data[, 1] == r_t_g1[i])
  
  # Checking if the gene ID was found in the "data" DataFrame
  if (length(index) > 0) {
    # If found, storing the corresponding gene symbol in the gene_symbols vector
    gene_symbols[i] <- data[index, 3]
  } else {
    # If not found, leaving the gene symbol as NA
    gene_symbols[i] <- NA
  }
}

# Replacing NA values with original gene IDs
gene_symbols[is.na(gene_symbols)] <- r_t_g1[is.na(gene_symbols)]

# Now gene_symbols contains the gene symbols corresponding to gene IDs in r_t_g1  ---- for pathway analysis

#MA plot showing statistically significant genes
plotMA(res1, main = "MA Plot to obtain DEGs")
legend("topright", legend = c("Not Significant", "Significant"), col = c("grey", "blue"), pch = c(1, 1))

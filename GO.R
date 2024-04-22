#load libraries
library(BiocManager)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

#taking significant genes
#genes_to_test <- rownames(res_sig[res_sig$log2FoldChange>0.5,])
genes_to_test <- gene_symbols
Go_results11 <- enrichGO(gene = genes_to_test, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")#biological processes
Go_results22 <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",  ont = "MF")#molecular functions
#c11<- as.data.frame(Go_results11)
#c12<- as.data.frame(Go_results22)

#plots
barplot(Go_results11, showCategory = 20, title = "Biological processes")
barplot(Go_results22, showCategory = 20, title = "Molecular Functions")

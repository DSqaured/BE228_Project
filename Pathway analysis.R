# Load required libraries if not already loaded
# install.packages("org.Hs.eg.db")
# install.packages("clusterProfiler")
library(org.Hs.eg.db)
library(clusterProfiler)

# Perform pathway enrichment analysis using enrichKEGG function
enrich_result1 <- enrichKEGG(gene = r_t_g1,
                             organism = 'hsa',  # Set to 'hsa' for Homo sapiens
                             pvalueCutoff = 0.05,  # Set significance cutoff
                             pAdjustMethod = "BH", # Benjamini-Hochberg method
                             qvalueCutoff = 0.05)  # Set adjusted p-value cutoff

# Write the enrichment results to a CSV file
write.csv(enrich_result1, file = "enrich_result_hsa.csv", row.names = TRUE)


# Visualize the enriched pathways
barplot(enrich_result1, showCategory = 20, title = "Top Enriched pathways", font.size = 10)  # Show the top 20 enriched pathways

#cnetplot(enrich_result1)
#dotplot(enrich_result1, showCategory = c("Signal transduction"))  # Show the top 20 enriched pathways

# Filter the enrich_result to contain only signal transduction pathways
#signal_transduction_enrichment <- enrich_result1[grepl("Signal Transduction", er1$subcategory)]

# Visualize the enriched signal transduction pathways
#dotplot(signal_transduction_enrichment, title = "Enriched Signal Transduction Pathways", font.size = 10)

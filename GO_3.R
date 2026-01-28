library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)


load("downstream_inputs/DE_DKO_adj_inputs.RData")

# GO enrichment - fyi it does NOT account for gene counts, only labels -----------
# GO enrichment - all DE genes
ego_DKO_vs_WT <- enrichGO(
  gene = rownames(DE_DKO_adj), # input vector of genes that appeared in our DE DEult
  universe = counts_file$Geneid, # vector of all genes that appeared in our experiment
  OrgDb = org.Mm.eg.db, # mouse ENSEMBL and GO term db
  keyType = "ENSEMBL", # our identifier we are providing
  ont = "BP", # "biological pathway" ontology branch
  pAdjustMethod = "BH", # Benjamini–Hochberg FDR
  pvalueCutoff = 0.05,
  readable = TRUE # converts ENSEMBL IDs to readable gene symbols
)

# GO enrichment - upregulated genes
ego_up <- enrichGO(
  gene = rownames(DE_DKO_adj_upregulated),
  universe = counts_file$Geneid,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# GO enrichment - downregulated genes
ego_down <- enrichGO(
  gene = rownames(DE_DKO_adj_downregulated),
  universe = counts_file$Geneid,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# GO dot plots ----------------------------------------------------------
# Plot all DE genes
all_DE_DKO_adj_dotplot <- dotplot(ego_DKO_vs_WT, showCategory = 15) +
  scale_y_discrete(labels = function(x) ifelse(nchar(x) > 40, 
                                               paste0(substr(x, 1, 37), "..."), 
                                               x)) +
  ggtitle("DKO Case vs Control") +
  labs(subtitle = "Adjusted with respect to WT case vs control")

# Plot upregulated genes
up_DE_DKO_adj_dotplot <- dotplot(ego_up, showCategory = 15) +
  scale_y_discrete(labels = function(x) ifelse(nchar(x) > 40, 
                                               paste0(substr(x, 1, 37), "..."), 
                                               x)) +
  ggtitle("Upregulated in DKO Case vs Control") +
  labs(subtitle = "Adjusted with respect to WT case vs control")

# Plot downregulated genes
down_DE_DKO_adj_dotplot <- dotplot(ego_down, showCategory = 15) +
  scale_y_discrete(labels = function(x) ifelse(nchar(x) > 40, 
                                               paste0(substr(x, 1, 37), "..."), 
                                               x)) +
  ggtitle("Downregulated in DKO Case vs Control") +
  labs(subtitle = "Adjusted with respect to WT case vs control")

# save outputs ------------------------------------------------------------
# save plot to pdf file
pdf(file = "plots/GO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches

print(all_DE_DKO_adj_dotplot)
print(up_DE_DKO_adj_dotplot)
print(down_DE_DKO_adj_dotplot)

dev.off()

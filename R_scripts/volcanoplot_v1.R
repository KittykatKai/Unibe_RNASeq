library(ggplot2)
library(enrichplot)
library(ggrepel)

load("DE_DKO_adj_inputs.RData")

# data framing ------------------------------------------------------------
# adding columns to original df for easier coloring
res_DKO_adj$direction <- "Not significant"
res_DKO_adj$direction[rownames(res_DKO_adj) %in% rownames(DE_DKO_adj_upregulated)] <- "Upregulated"
res_DKO_adj$direction[rownames(res_DKO_adj) %in% rownames(DE_DKO_adj_downregulated)] <- "Downregulated"

# creating new df out of DE_DKO_adj_renamed (gene IDs instead of ENSEMBL) for sig genes labeling
DE_DKO_adj_df <- as.data.frame(DE_DKO_adj_renamed)
volcano_labels <- DE_DKO_adj_df[DE_DKO_adj_df$padj < 0.0001 & abs(DE_DKO_adj_df$log2FoldChange) > 2,]

# volcano plot ------------------------------------------------------------
DE_DKO_adj_volcano <- ggplot(data = res_DKO_adj, aes(x = log2FoldChange, y = -log10(padj), col = direction)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  geom_text_repel(data = volcano_labels, 
                  aes(label=rownames(volcano_labels)), size = 3, 
                  show.legend = FALSE) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  coord_cartesian(ylim = c(0, 80), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Direction', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"),
       title = 'DKO Case vs Control (adjusted with respect to WT case vs control)', subtitle = 'padj < 0.0001 and absolute log2FoldChange > 2') +
  scale_x_continuous(breaks = seq(-10, 10, 2)) # to customise the breaks in the x axis

# save outputs ------------------------------------------------------------
# save plot to pdf file
pdf(file = "plots/volcanoplot.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches

print(DE_DKO_adj_volcano)

dev.off()

# nuclear option
save.image("volcanoplot.RData")
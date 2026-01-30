library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Mm.eg.db)

# Define directories and load data ----------------------------------------

input_dir <- "/Users/kaifrances/Unibe_RNASeqR/R_scripts/downstream_inputs"
plot_dir <- "/Users/kaifrances/Unibe_RNASeqR/R_scripts/plots"

load(file.path(input_dir, "DE_WT_inputs.RData"))

# Define genes of interest and map to Ensembl IDs -------------------------

genes_of_interest <- c("Ifit1", "Gbp2", "Gbp5", "Cxcl9")

gene_mapping <- AnnotationDbi::select(org.Mm.eg.db,
                                      keys = genes_of_interest,
                                      columns = c("ENSEMBL", "SYMBOL"),
                                      keytype = "SYMBOL")

# Define sample metadata --------------------------------------------------

sample_metadata <- data.frame(
  sample = c("SRR7821949", "SRR7821950", "SRR7821951", "SRR7821952", "SRR7821953",
             "SRR7821968", "SRR7821969", "SRR7821970",
             "SRR7821954", "SRR7821955", "SRR7821956", "SRR7821957",
             "SRR7821971", "SRR7821972", "SRR7821973"),
  genotype = c(rep("WT", 8), rep("DKO", 7)),
  condition = c(rep("Infected", 5), rep("Control", 3),
                rep("Infected", 4), rep("Control", 3))
)

# Extract counts for genes of interest ------------------------------------

sample_names <- gsub(".*/|\\.sorted\\.bam", "", colnames(counts_file)[7:21])

gene_counts_df <- do.call(rbind, lapply(1:nrow(gene_mapping), function(i) {
  idx <- which(counts_file$Geneid == gene_mapping$ENSEMBL[i])
  if (length(idx) > 0) {
    data.frame(
      gene = gene_mapping$SYMBOL[i],
      sample = sample_names,
      counts = as.numeric(counts_file[idx[1], 7:21])
    )
  }
}))

gene_counts_df <- merge(gene_counts_df, sample_metadata, by = "sample")
gene_counts_df$group <- factor(paste(gene_counts_df$genotype, gene_counts_df$condition, sep = "\n"),
                               levels = c("WT\nControl", "WT\nInfected", "DKO\nControl", "DKO\nInfected"))

# Calculate summary statistics --------------------------------------------

summary_stats <- gene_counts_df %>%
  group_by(gene, genotype, condition) %>%
  summarise(mean_counts = round(mean(counts), 1),
            sd_counts = round(sd(counts), 1),
            n = n(), .groups = "drop")

report_table <- summary_stats %>%
  mutate(Group = paste(genotype, condition),
         Value = paste0(mean_counts, " ± ", sd_counts)) %>%
  dplyr::select(gene, Group, Value) %>%
  pivot_wider(names_from = Group, values_from = Value)

# Create plot -------------------------------------------------------------

p <- ggplot(gene_counts_df, aes(x = group, y = counts, fill = interaction(genotype, condition))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 12, face = "italic"),
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  labs(title = "Expression of Interferon-Response Genes",
       subtitle = "Raw counts from RNA-seq",
       x = "", y = "Raw Counts") +
  scale_fill_manual(values = c("WT.Control" = "#4DAF4A", "WT.Infected" = "#E41A1C",
                               "DKO.Control" = "#377EB8", "DKO.Infected" = "#984EA3"))

# Save outputs ------------------------------------------------------------

write.csv(summary_stats, file.path(input_dir, "genes_of_interest_summary.csv"), row.names = FALSE)
write.csv(report_table, file.path(input_dir, "genes_of_interest_report_table.csv"), row.names = FALSE)
ggsave(file.path(plot_dir, "genes_of_interest_expression.pdf"), p, width = 8, height = 8)

print(p)
print(report_table)
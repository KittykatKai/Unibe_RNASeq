library("DESeq2")
library(biomaRt)

# loading in featureCounts data -------------------------------------------

# load counts.txt as a data frame
counts_file <- read.table("counts.txt", header=TRUE, comment.char="#", check.names=FALSE)
# create df of raw counts, no annotation info
raw_counts_df <- counts_file[7:ncol(counts_file)]
# replace shitty full path colnames using pre-made samplenames.txt
samplenames <- read.delim('samplenames.txt')
colnames(raw_counts_df) <- samplenames$Sample
# set row names = the gene ID
rownames(raw_counts_df) <- counts_file$Geneid


# create colData ----------------------------------------------------------

# trim the name to isolate the type (WT or KO) and condition (case or control)
type_trim <- sub("_.*", "", sub("Blood_", "", samplenames$Group))
condition_trim <- sub(".*_", "", samplenames$Group)

# construct coldata df, strings as factors
coldata <- data.frame(
  sample = samplenames$Sample,
  type = type_trim,
  condition = condition_trim,
  stringsAsFactors = TRUE
)

# also need to set the "reference level" for DESeq to know what our control is. 
# R automatically picks a reference level via alphabetical order, so use relevel to explicitly set ref level
coldata$type <- relevel(coldata$type, ref = "WT")
coldata$condition <- relevel(coldata$condition, ref = "Control")


# DESeqDataSetFromMatrix, design term -------------------------------------

# make dds = "DEseq data set" object
dds <- DESeqDataSetFromMatrix(countData = raw_counts_df,
                              colData = coldata,
                              design = ~ type + condition + type:condition)
# NOTE: the design is two factored, controls for type
# "type" = difference between WT and DKO at the reference condition, control
# "condition" = difference between case and control, allows model to adjust it by the baseline difference that we already observed between WT and DKO
# "type:condition" means the "condition" output will specifically show "changes based on condition, for each type, compared to the reference type (which is WT)". it also adds a new resultName which lets us specifically view "changes based on condition, for DKO, beyond the reference type"

dds <- DESeq(dds)

# DESeq initial visualizations: DispEsts, pca -------------------

# open pdf device for saving plots
pdf(file = "plots/deseq.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches

# variance stabilizing transform, plot disp ests for quality control / fit
ntd <- vst(dds, blind=TRUE) # ntd = "normalized transformed data"
plotDispEsts(dds, main = "dispersion estimates following VST")

# pca
plot_pca <- plotPCA(ntd, intgroup=c("type","condition")) + # intgroup just chooses how to color the dots, no impact on PCA
  ggtitle("PCA of normalized transformed DESeq2 data (VST)") 
print(plot_pca)

# DESeq results, cross-comparisons ----------------------------------------

resultsNames(dds) # check this output. it will tell you names of each valid effect / comparison that DESeq found

# checking res: "effect of case vs control, for condition (WT)". due to the interaction it only extracts the main effect by default, which is only in WT
res_WT <- results(dds, name = "condition_Case_vs_Control")
# filtering by padj < 0.05 and log2foldchange greater than 1
DE_WT <- res_WT[!is.na(res_WT$padj) & res_WT$padj < 0.05 & abs(res_WT$log2FoldChange) > 1, ]

# checking res: "effect of case vs control, for DKO". this is the BULK difference and is not adjusted by the observations in WT
res_DKO <- results(dds, contrast = list(
  c("condition_Case_vs_Control", "typeDKO.conditionCase")))
# filtering by padj < 0.05 and log2foldchange greater than 1
DE_DKO <- res_DKO[!is.na(res_DKO$padj) & res_DKO$padj < 0.05 & abs(res_DKO$log2FoldChange) > 1, ]


# checking res for interaction term: "effect of case vs control, for DKO, ADJUSTED by effect in WT". this is the ADDITIONAL difference, beyond what changes in WT
res_DKO_adj <- results(dds, name = "typeDKO.conditionCase")

# filtering by padj < 0.05 and log2foldchange greater than 1
DE_DKO_adj <- res_DKO_adj[!is.na(res_DKO_adj$padj) & res_DKO_adj$padj < 0.05 & abs(res_DKO_adj$log2FoldChange) > 1, ]

# split df into up and down regulated
DE_DKO_adj_upregulated <- DE_DKO_adj[DE_DKO_adj$log2FoldChange > 1, ]
DE_DKO_adj_downregulated <- DE_DKO_adj[DE_DKO_adj$log2FoldChange < -1, ]

# adding columns to original df for easier coloring
DE_DKO_adj$direction[rownames(DE_DKO_adj) %in% rownames(DE_DKO_adj_upregulated)] <- "Upregulated"
DE_DKO_adj$direction[rownames(DE_DKO_adj) %in% rownames(DE_DKO_adj_downregulated)] <- "Downregulated"

res_DKO_adj$direction <- "Not significant"
res_DKO_adj$direction[rownames(res_DKO_adj) %in% rownames(DE_DKO_adj_upregulated)] <- "Upregulated"
res_DKO_adj$direction[rownames(res_DKO_adj) %in% rownames(DE_DKO_adj_downregulated)] <- "Downregulated"

# Check how many up vs down
paste0("Up-regulated in DKO: ", dim(DE_DKO_adj_upregulated)[1])
paste0("Down-regulated in DKO: ", dim(DE_DKO_adj_downregulated)[1])

# new results objects: replacing ENSMUS whatever BS ensembl IDs with real names ----------------
# connect to Ensembl
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

gene_ids <- rownames(dds)

# map Ensembl IDs to gene symbols
annot <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),
               filters='ensembl_gene_id',
               values=gene_ids,
               mart=ensembl)

# make a named vector
gene_symbols <- annot$mgi_symbol
names(gene_symbols) <- annot$ensembl_gene_id

# Replace rownames
DE_DKO_adj_mapped_symbols <- gene_symbols[rownames(DE_DKO_adj)]
DE_DKO_adj_renamed <- DE_DKO_adj
rownames(DE_DKO_adj_renamed) <- ifelse(is.na(DE_DKO_adj_mapped_symbols), rownames(DE_DKO_adj), DE_DKO_adj_mapped_symbols)

# save all workspace objects ----------------------------------------------
# save pdf device
dev.off()

# minimal inputs needed for analysis of DE_DKO_adj in GO.R, volcanoplot.R
save(res_DKO_adj,
     DE_DKO_adj,
     DE_DKO_adj_upregulated,
     DE_DKO_adj_downregulated,
     DE_DKO_adj_renamed,
     counts_file,
     file = "DE_DKO_adj_inputs.RData")

# nuclear option
save.image("deseq.RData")

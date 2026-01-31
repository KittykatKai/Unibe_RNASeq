library("DESeq2")
library(biomaRt)
library(ggplot2)

# loading in featureCounts data -------------------------------------------

# load counts.txt as a data frame
counts_file <- read.table("raw_inputfiles/counts.txt", header=TRUE, comment.char="#", check.names=FALSE)
# create df of raw counts, no annotation info
raw_counts_df <- counts_file[7:ncol(counts_file)]
# replace full path colnames using pre-made samplenames.txt
samplenames <- read.delim('raw_inputfiles/samplenames.txt')
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

# set up ENSEMBL -> gene symbol lookup for making human-readable labels --------
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

# DESeq results, cross-comparisons ----------------------------------------

resultsNames(dds) # check this output. it will tell you names of each valid effect / comparison that DESeq found

### define extractres() function to extract nice results for each DESeq2 contrast name ###
extractres <- function(contrastTerm) {
  # inputs: contrastTerm (str or vector): build it using resultsNames outputs only
  # outputs a list:
  #   res: DESeqResultObject
  #   DE: DESeqResultObject filtered by p < 0.05 and abs(log2fc) > 1, "differentially expressed" genes
  #   DE_upregulated: DE filtered by log2fc > 1, "upregulated" genes
  #   DE_downregulated: DE filtered by log2fc < -1, "downregulated" genes
  #   DE_renamed: copy of DE df with ENSEMBL IDs replaced by human-readable gene symbols, for easier labeling
  
  # filtering by padj < 0.05 and log2foldchange greater than 1
  res <- results(dds, contrast = list(contrastTerm))
  DE <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  
  # split df into up and down regulated
  DE_upregulated <- DE[DE$log2FoldChange > 1, ]
  DE_downregulated <- DE[DE$log2FoldChange < -1, ]
  
  # adding directional columns to original df for easier coloring
  DE$direction[rownames(DE) %in% rownames(DE_upregulated)] <- "Upregulated"
  DE$direction[rownames(DE) %in% rownames(DE_downregulated)] <- "Downregulated"
  
  res$direction <- "Not significant"
  res$direction[rownames(res) %in% rownames(DE_upregulated)] <- "Upregulated"
  res$direction[rownames(res) %in% rownames(DE_downregulated)] <- "Downregulated"
  
  # create human-readable gene symbols df for easier labeling
  DE_mapped_symbols <- gene_symbols[rownames(DE)]
  DE_renamed <- DE
  rownames(DE_renamed) <- ifelse(is.na(DE_mapped_symbols), rownames(DE), DE_mapped_symbols)
  
  return(list(res = res, 
              DE = DE, 
              DE_upregulated = DE_upregulated, 
              DE_downregulated = DE_downregulated,
              DE_renamed = DE_renamed))
}

# checking res: "effect of DKO vs WT, for healthy case (default condition)". 
res_healthy_list <- extractres(c("type_DKO_vs_WT"))

# checking res: "effect of DKO vs WT, for diseased case". 
res_diseased_list <- extractres(c("type_DKO_vs_WT", "typeDKO.conditionCase"))

# checking res: "effect of case vs control, for condition (WT)". due to the interaction it only extracts the main effect by default, which is only in WT
res_WT_list <- extractres("condition_Case_vs_Control")

# checking res: "effect of case vs control, for DKO". this is the BULK difference and is not adjusted by the observations in WT
res_DKO_list <- extractres(c("condition_Case_vs_Control", "typeDKO.conditionCase"))

# checking res for interaction term: "effect of case vs healthy control, for DKO, ADJUSTED by effect in WT". this is the ADDITIONAL difference, beyond what changes in WT
res_DKO_adj_list <- extractres("typeDKO.conditionCase")

# save all workspace objects ----------------------------------------------
# save pdf device for plots
dev.off()

# minimal inputs needed for analysis of DE_DKO_adj in GO.R, volcanoplot.R
res_DKO_adj <- res_DKO_adj_list$res
DE_DKO_adj <- res_DKO_adj_list$DE
DE_DKO_adj_upregulated <- res_DKO_adj_list$DE_upregulated
DE_DKO_adj_downregulated <- res_DKO_adj_list$DE_downregulated
DE_DKO_adj_renamed <- res_DKO_adj_list$DE_renamed

save(res_DKO_adj,DE_DKO_adj, DE_DKO_adj_upregulated, DE_DKO_adj_downregulated, DE_DKO_adj_renamed, counts_file,
     file = "downstream_inputs/DE_DKO_adj_inputs.RData")

# save res_healthy
res_healthy <- res_healthy_list$res
DE_healthy <- res_healthy_list$DE
DE_healthy_upregulated <- res_healthy_list$DE_upregulated
DE_healthy_downregulated <- res_healthy_list$DE_downregulated
DE_healthy_renamed <- res_healthy_list$DE_renamed

save(res_healthy,DE_healthy, DE_healthy_upregulated, DE_healthy_downregulated, DE_healthy_renamed, counts_file,
     file = "downstream_inputs/DE_healthy_inputs.RData")

# save res_diseased
res_diseased <- res_diseased_list$res
DE_diseased <- res_diseased_list$DE
DE_diseased_upregulated <- res_diseased_list$DE_upregulated
DE_diseased_downregulated <- res_diseased_list$DE_downregulated
DE_diseased_renamed <- res_diseased_list$DE_renamed

save(res_diseased,DE_diseased, DE_diseased_upregulated, DE_diseased_downregulated, DE_diseased_renamed, counts_file,
     file = "downstream_inputs/DE_diseased_inputs.RData")

# save res_WT
res_WT <- res_WT_list$res
DE_WT <- res_WT_list$DE
DE_WT_upregulated <- res_WT_list$DE_upregulated
DE_WT_downregulated <- res_WT_list$DE_downregulated
DE_WT_renamed <- res_WT_list$DE_renamed

save(res_WT,DE_WT, DE_WT_upregulated, DE_WT_downregulated, DE_WT_renamed, counts_file,
     file = "downstream_inputs/DE_WT_inputs.RData")

# save res_DKO
res_DKO <- res_DKO_list$res
DE_DKO <- res_DKO_list$DE
DE_DKO_upregulated <- res_DKO_list$DE_upregulated
DE_DKO_downregulated <- res_DKO_list$DE_downregulated
DE_DKO_renamed <- res_DKO_list$DE_renamed

save(res_DKO,DE_DKO, DE_DKO_upregulated, DE_DKO_downregulated, DE_DKO_renamed, counts_file,
     file = "downstream_inputs/DE_DKO_inputs.RData")

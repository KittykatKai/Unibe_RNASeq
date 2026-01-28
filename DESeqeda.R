library("DESeq2")
library("pheatmap")
library("biomaRt")
library("org.Mm.eg.db")
library("clusterProfiler")

BiocManager::install("org.Mm.eg.db")

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

# replacing ENSMUS whatever BS ensembl IDs with real names ----------------

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

# replace rownames with gene symbols
rownames(dds) <- gene_symbols[rownames(dds)]


# DESeq initial visualizations: DispEsts, pheatmap, pca -------------------

# variance stabilizing transform, plot disp ests for quality control / fit
ntd <- vst(dds, blind=TRUE) # ntd = "normalized transformed data"
plotDispEsts(dds, main = "dispersion estimates following VST")

# heatmap of topgenes
# topgenes top 100 expressed genes
topgenes_exp <- order(rowMeans(counts(dds, normalized=TRUE)),
                      decreasing=TRUE)[1:100]

# Get the variance of each gene across samples
gene_var <- apply(assay(ntd), 1, var) # use assay() to extract matrix from ntd... cuz ntd is not a data frame...

# topgenes top 20 VARIANT genes
topgenes_var <- order(gene_var,
                      decreasing=TRUE)[1:20]

annotation_df <- as.data.frame(colData(dds)[,c("type","condition")])
# make pheatmap
pheatmap(assay(ntd)[topgenes_var,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=annotation_df, main='top 20 variant genes')

# pca
plotPCA(ntd, intgroup=c("type","condition")) # intgroup just chooses how to color the dots, no impact on PCA

# DESeq results, cross-comparisons ----------------------------------------

resultsNames(dds) # check this output. it will tell you names of each valid effect / comparison that DESeq found

# checking res: "effect of case vs control, for condition (WT)". due to the interaction it only extracts the main effect by default, which is only in WT
res_WT <- results(dds, name = "condition_Case_vs_Control")

DE_WT <- res_WT[!is.na(res_WT$padj) & res_WT$padj < 0.05, ]
DEsum_WT <- sum(res_WT$padj < 0.05, na.rm=TRUE)
up_WT <- sum(DE_WT$log2FoldChange > 0, na.rm=TRUE)
down_WT <- sum(DE_WT$log2FoldChange < 0, na.rm=TRUE)

# checking res: "effect of case vs control, for DKO". this is the BULK difference and is not adjusted by the observations in WT
res_DKO <- results(dds, contrast = list(
  c("condition_Case_vs_Control", "typeDKO.conditionCase")))

DE_DKO <- res_DKO[!is.na(res_DKO$padj) & res_DKO$padj < 0.05, ]
DEsum_DKO <- sum(res_DKO$padj < 0.05, na.rm=TRUE)
up_DKO <- sum(DE_DKO$log2FoldChange > 0, na.rm=TRUE)
down_DKO <- sum(DE_DKO$log2FoldChange < 0, na.rm=TRUE)


# checking res for interaction term: "effect of case vs control, for DKO, ADJUSTED by effect in WT". this is the ADDITIONAL difference, beyond what changes in WT
res_DKO_adj <- results(dds, name = "typeDKO.conditionCase")

DE_DKO_adj <- res_DKO_adj[!is.na(res_DKO_adj$padj) & res_DKO_adj$padj < 0.05, ]
DEsum_DKO_adj <- sum(res_DKO_adj$padj < 0.05, na.rm=TRUE)
up_DKO_adj <- sum(DE_DKO_adj$log2FoldChange > 0, na.rm=TRUE)
down_DKO_adj <- sum(DE_DKO_adj$log2FoldChange < 0, na.rm=TRUE)

res_summary <- rbind(
  c(up_WT, down_WT, DEsum_WT),
  c(up_DKO, down_DKO, DEsum_DKO),
  c(up_DKO_adj, down_DKO_adj, DEsum_DKO_adj)
)

rownames(res_summary) <- c("WT", "DKO", "DKO - WT")
colnames(res_summary) <- c("upregulated", "downregulated", "total count")

print(res_summary)

paste0('DE genes caused by disease in WT: ', DEsum_WT, '. upregulated: ', up_WT, ', downregulated: ', down_WT)
paste0('DE genes caused by disease in DKO: ', DEsum_DKO, '. upregulated: ', up_DKO, ', downregulated: ', down_DKO)
paste0('DE genes caused by disease in DKO, novel as compared to WT: ', DEsum_DKO_adj, '. upregulated: ', up_DKO_adj, ', downregulated: ', down_DKO_adj)

save.image("deseq.RData")


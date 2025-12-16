This repository contains code for RNA Sequencing
All analyses are carried out using a 2-factor design

### /R_scripts
Raw data input files for deseq_v2.R:
1. **counts.txt**: contains read counts, output from FeatureCounts
2. **samplenames.txt**: copy-pasted tab-delimited table of sample IDs and experimental group (ex. "Blood_WT_Case")

Scripts:
1. **deseq_v2.R**: DESeq2 2-factor design (WT/DKO, Case/Control): type + condition + type:condition
    a. outputs **DE_DKO_adj_inputs.RData** which contains the minimally required objects for all downstream scripts. this only contains data for DE_DKO_adj, which is differential expression only within the interaction term of DKO case vs control (subtracted from WT case vs control).
2. **GO_v3.R**: GO enrichment analysis and dotplots for all DE genes. repeats the analysis after subsetting the data into upregulated and downregulated. works on DE_DKO_adj
3. **volcanoplot_v1.R**: volcano plot for all DE genes in DE_DKO_adj

### /plots
deseq:
1. DispEsts plot to assess good fit after variance stabilizing transform (VST)
2. PCA

GO:
1. GO enrichment analysis dotplot of all DE genes in DE_DKO_adj
2. GO enrichment analysis dotplot of only upregulated genes
3. GO enrichment analysis dotplot of only downregulated genes

volcanoplot:
1. volcano plot of DE_DKO_adj
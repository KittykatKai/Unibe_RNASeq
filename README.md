This repository contains code for RNA Sequencing


## Part 1: BASH Scripts

Pipeline:
1. Trim (Calliefastp)
2. Quality Control (Calliefastqc & Calliemultiqc)
3. Map reads to reference genome (Calliehisat)
4. Quality Control (Calliefastqc & Calliemultiqc)
5. Count number of reads per gene via FeatureCounts (Calliesubread)


## Part 2: R Scripts 

All analyses are carried out using a 2-factor design


**/raw_inputfiles**

Raw data input files are used by deseq_v3.R and GO_v3.R:

- counts.txt: contains read counts, output from FeatureCounts
- counts.txt.summary: Contains a summarised version of counts.txt for ease of use
- samplenames.txt: copy-pasted tab-delimited table of sample IDs and experimental group (ex. "Blood_WT_Case")

Pipline:
1. Deseq_v3.R: DESeq2 2-factor design (WT/DKO, Case/Control): type + condition + type:condition. Saves DESeqResultsObjects for all contrasts into /downstream_inputs that can be used for downstream scripts
2. GO_v3.R: GO enrichment analysis and dotplots for DE_DKO_adj. Creates three plots: one for all genes, one within subset of downregulated genes, one within subset of upregulated genes
3. GeneExtractor.R: Creates counts, summary statistics, and boxplots for pre selected genes of interest.
4. volcanoplot_v1.R: creates volcano plots for DE_DKO_adj and DE_diseased. Also creates an alternate version which overlays the significant type I IFN genes on the plot.

**/downstream_inputs**

- DE_healthy_inputs.RData: "effect of DKO vs WT, for healthy case". type_DKO_vs_WT
- DE_diseased_inputs.RData: "effect of DKO vs WT, for diseased case". c("type_DKO_vs_WT", "typeDKO.conditionCase")
- DE_WT_inputs.RData: "effect of case vs control, for WT". condition_Case_vs_Control
- DE_DKO_inputs.RData: "effect of case vs control, for DKO". this is the BULK difference and is not adjusted by the observations in WT. c("condition_Case_vs_Control", "typeDKO.conditionCase")
- DE_DKO_adj_inputs.RData: "effect of case vs healthy control, for DKO, ADJUSTED by effect in WT". this is the ADDITIONAL difference, beyond what changes in WT. typeDKO.conditionCase

**/plots**

Location for plots created by R Scripts

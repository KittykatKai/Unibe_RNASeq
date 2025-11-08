#!/bin/bash
#SBATCH --job-name=hisat2_mapping
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=pibu_el8
#SBATCH --output=hisat2_%j.out
#SBATCH --error=hisat2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Container path
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Define paths
GENOME_FA="/data/users/ksales/Unibe_RNASeq/Mus_musculus.GRCm39.dna.primary_assembly.fa"
ANNOTATION_GTF="/data/users/ksales/Unibe_RNASeq/Mus_musculus.GRCm39.115.gtf"
INDEX_DIR="/data/users/ksales/Unibe_RNASeq/hisat2index"
INDEX_PREFIX="${INDEX_DIR}/genome"
READS_DIR="/data/users/ksales/Unibe_RNASeq/trimmedreads_Blood"
OUTPUT_DIR="/data/users/ksales/Unibe_RNASeq/hisat2_output"

# Build HISAT2 index
apptainer exec --bind /data ${CONTAINER} hisat2-build -p ${SLURM_CPUS_PER_TASK} ${GENOME_FA} ${INDEX_PREFIX}

# Loop through paired-end samples
for R1 in ${READS_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename ${R1} _1.fastq.gz)
    R2="${READS_DIR}/${SAMPLE}_2.fastq.gz"

    # Run HISAT2 with RF strandedness
    apptainer exec --bind /data ${CONTAINER} hisat2 -p ${SLURM_CPUS_PER_TASK} \
        --rna-strandness RF \
        -x ${INDEX_PREFIX} \
        -1 ${R1} \
        -2 ${R2} \
        -S ${OUTPUT_DIR}/${SAMPLE}.sam \
        --summary-file ${OUTPUT_DIR}/${SAMPLE}_summary.txt

    # Convert SAM to BAM
    apptainer exec --bind /data ${CONTAINER} samtools view -@ ${SLURM_CPUS_PER_TASK} -bS ${OUTPUT_DIR}/${SAMPLE}.sam -o ${OUTPUT_DIR}/${SAMPLE}.bam

    # Sort BAM
    apptainer exec --bind /data ${CONTAINER} samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${OUTPUT_DIR}/${SAMPLE}.sorted.bam ${OUTPUT_DIR}/${SAMPLE}.bam

    # Index BAM file
    apptainer exec --bind /data ${CONTAINER} samtools index ${OUTPUT_DIR}/${SAMPLE}.sorted.bam

    # Remove intermediate files to save space
    rm ${OUTPUT_DIR}/${SAMPLE}.sam
    rm ${OUTPUT_DIR}/${SAMPLE}.bam

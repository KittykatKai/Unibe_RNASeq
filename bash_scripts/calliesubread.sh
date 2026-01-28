#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=pibu_el8
#SBATCH --output=featureCounts_%j.out
#SBATCH --error=featureCounts_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Container path
CONTAINER="/containers/apptainer/subread_2.0.6.sif"

# Define paths
ANNOTATION_GTF="/data/users/ksales/Unibe_RNASeq/Mus_musculus.GRCm39.115.gtf"
BAM_DIR="/data/users/ksales/Unibe_RNASeq/hisat2/hisat2_output/attempt3"
OUTPUT_DIR="/data/users/ksales/Unibe_RNASeq/subread_output"
OUTPUT_FILE="${OUTPUT_DIR}/counts.txt"

# Run featureCounts with explicit file list
apptainer exec --bind /data ${CONTAINER} featureCounts \
    -p \ # Paired-end mode
    -s 2 \ # Reverse stranded RF library
    -T ${SLURM_CPUS_PER_TASK} \ # Number of threads 
    -a ${ANNOTATION_GTF} \ # Annotation file
    -o ${OUTPUT_FILE} \ # Output file 
    ${BAM_DIR}/SRR7821949.sorted.bam \
    ${BAM_DIR}/SRR7821950.sorted.bam \
    ${BAM_DIR}/SRR7821951.sorted.bam \
    ${BAM_DIR}/SRR7821952.sorted.bam \
    ${BAM_DIR}/SRR7821953.sorted.bam \
    ${BAM_DIR}/SRR7821968.sorted.bam \
    ${BAM_DIR}/SRR7821969.sorted.bam \
    ${BAM_DIR}/SRR7821970.sorted.bam \
    ${BAM_DIR}/SRR7821954.sorted.bam \
    ${BAM_DIR}/SRR7821955.sorted.bam \
    ${BAM_DIR}/SRR7821956.sorted.bam \
    ${BAM_DIR}/SRR7821957.sorted.bam \
    ${BAM_DIR}/SRR7821971.sorted.bam \
    ${BAM_DIR}/SRR7821972.sorted.bam \
    ${BAM_DIR}/SRR7821973.sorted.bam
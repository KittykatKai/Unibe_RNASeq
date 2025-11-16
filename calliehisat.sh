#!/bin/bash
#SBATCH --job-name=hisat2_full
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=pibu_el8
#SBATCH --output=hisat2_%A_%a.out
#SBATCH --error=hisat2_%A_%a.err
#SBATCH --array=0-14
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Container path
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Define paths
GENOME_FA="/data/users/ksales/Unibe_RNASeq/Mus_musculus.GRCm39.dna.primary_assembly.fa"
INDEX_DIR="/data/users/ksales/Unibe_RNASeq/hisat2index/attempt3"
INDEX_PREFIX="${INDEX_DIR}/genome"
READS_DIR="/data/users/ksales/Unibe_RNASeq/trimmedreads_Blood"
OUTPUT_DIR="/data/users/ksales/Unibe_RNASeq/hisat2_output/attempt3"

# Build index only once (only array task 0 does this)
INDEX_LOCK="${INDEX_DIR}/.index_complete"
if [[ ${SLURM_ARRAY_TASK_ID} -eq 0 ]]; then
    # Check if index already exists
    if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
        apptainer exec --bind /data ${CONTAINER} hisat2-build \
            -p ${SLURM_CPUS_PER_TASK} \
            ${GENOME_FA} \
            ${INDEX_PREFIX}
        
        # Create lock file to signal completion
        touch ${INDEX_LOCK}
    else
        touch ${INDEX_LOCK}
    fi
else
    # Other tasks wait for task 0 to finish building the index
    while [[ ! -f ${INDEX_LOCK} ]]; do
        sleep 30
    done
    # Extra wait to ensure all index files are fully written
    sleep 10
fi

# Array of sample IDs
samples=(
    SRR7821949
    SRR7821950
    SRR7821951
    SRR7821952
    SRR7821953
    SRR7821968
    SRR7821969
    SRR7821970
    SRR7821954
    SRR7821955
    SRR7821956
    SRR7821957
    SRR7821971
    SRR7821972
    SRR7821973
)

# Get sample for THIS array task
SAMPLE=${samples[$SLURM_ARRAY_TASK_ID]}

# Define read files for this sample
R1="${READS_DIR}/${SAMPLE}_1_trimmed.fastq.gz"
R2="${READS_DIR}/${SAMPLE}_2_trimmed.fastq.gz"

# Check if files exist
if [[ ! -f ${R1} ]] || [[ ! -f ${R2} ]]; then
    exit 1
fi

# Run HISAT2 with RF strandedness
apptainer exec --bind /data ${CONTAINER} hisat2 \
    -p ${SLURM_CPUS_PER_TASK} \
    --rna-strandness RF \
    -x ${INDEX_PREFIX} \
    -1 ${R1} \
    -2 ${R2} \
    -S ${OUTPUT_DIR}/${SAMPLE}.sam \
    --summary-file ${OUTPUT_DIR}/${SAMPLE}_summary.txt

# Convert SAM to BAM
apptainer exec --bind /data ${CONTAINER} samtools view \
    -@ ${SLURM_CPUS_PER_TASK} -hbS \
    ${OUTPUT_DIR}/${SAMPLE}.sam > ${OUTPUT_DIR}/${SAMPLE}.bam

# Sort BAM
apptainer exec --bind /data ${CONTAINER} samtools sort \
    -m 4G -@ ${SLURM_CPUS_PER_TASK} \
    -o ${OUTPUT_DIR}/${SAMPLE}.sorted.bam \
    -T ${OUTPUT_DIR}/temp_${SAMPLE} \
    ${OUTPUT_DIR}/${SAMPLE}.bam

# Index BAM file
apptainer exec --bind /data ${CONTAINER} samtools index \
    ${OUTPUT_DIR}/${SAMPLE}.sorted.bam

# Remove intermediate files to save space
rm ${OUTPUT_DIR}/${SAMPLE}.sam
rm ${OUTPUT_DIR}/${SAMPLE}.bam
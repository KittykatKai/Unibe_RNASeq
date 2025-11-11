#!/bin/bash
#SBATCH --job-name=sample1blood
#SBATCH --output=/data/users/ksales/Unibe_RNASeq/output_dir/attempt2/fastqc_%A_%a.out
#SBATCH --error=/data/users/ksales/Unibe_RNASeq/output_dir/attempt2/fastqc_%A_%a.err
#SBATCH --time=02:10:00
#SBATCH --mem=1G
#SBATCH --array=0-14
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Set paths
CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"
INPUT_DIR="/data/users/ksales/Unibe_RNASeq/trimmedreads_Blood"
OUTPUT_DIR="/data/users/ksales/Unibe_RNASeq/output_dir/attempt2"

# Array of sample IDs (using array in this situation is more efficent)
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

#  sample for array
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Run FastQC using Apptainer container
apptainer exec "$CONTAINER" fastqc \
    "$INPUT_DIR"/${sample}_1_trimmed.fastq.qz \
    "$INPUT_DIR"/${sample}_2_trimmed.fastq.qz \
    -o "$OUTPUT_DIR" \
    -t $SLURM_CPUS_PER_TASK
    
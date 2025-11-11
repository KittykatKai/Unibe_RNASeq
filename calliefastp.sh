#!/bin/bash
#SBATCH --job-name=fastp_bloods
#SBATCH --output=/data/users/ksales/output_dir/fastp_%A_%a.out
#SBATCH --error=/data/users/ksales/output_dir/fastp_%A_%a.err
#SBATCH --time=02:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --array=0-14
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kai.sales@students.unibe.ch

# Set Paths
CONTAINER="/containers/apptainer/fastp_0.24.1.sif"
INPUT_DIR="/data/users/ksales/reads_Blood"
OUTPUT_DIR="/data/users/ksales/Unibe_RNASeq/trimmedreads_Blood"

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

# Get sample for array
sample=${samples[$SLURM_ARRAY_TASK_ID]}

apptainer exec --bind /data:/data "$CONTAINER" fastp \
    -i "$INPUT_DIR"/${sample}_1.fastq.gz \
    -I "$INPUT_DIR"/${sample}_2.fastq.gz \
    -o "$OUTPUT_DIR"/${sample}_1_trimmed.fastq.gz \
    -O "$OUTPUT_DIR"/${sample}_2_trimmed.fastq.gz \
    -h "$OUTPUT_DIR"/${sample}_fastp.html \
    -j "$OUTPUT_DIR"/${sample}_fastp.json \
    --thread 4 \
    --detect_adapter_for_pe \
    --length_required 20 \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20
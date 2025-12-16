#!/bin/bash
#SBATCH --job-name=sample1blood
#SBATCH --output=/data/users/ksales/output_dir/multiqc_%j.out
#SBATCH --error=/data/users/ksales/output_dir/multiqc_%j.err
#SBATCH --time=02:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8

# Set paths
CONTAINER="/containers/apptainer/multiqc-1.19.sif"
INPUT_DIR="/data/users/ksales/output_dir"
OUTPUT_DIR="/data/users/ksales/output_dir"

# Run FastQC using Apptainer container
apptainer exec "$CONTAINER" multiqc \
    "$INPUT_DIR"/*fastqc.* 

echo "MultiQC  analysis complete!"



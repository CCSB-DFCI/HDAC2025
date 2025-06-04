#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=salmon_quant
#SBATCH --output=/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4/salmon_quant.log
#SBATCH --error=/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4/salmon_quant.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00

module load salmon/1.8.0

INDEX_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_reference_genome/salmon_index_grch38_r97"
FASTQ_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/fastq_bulk"
OUTPUT_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4"

mkdir -p "$OUTPUT_DIR"

cd "$FASTQ_DIR"

for R1 in *_R1_001.fastq.gz; do
    # Find matching R2
    R2="${R1/_R1_/_R2_}"

    # Extract sample name: e.g., DMSO-1_5days -> DMSO1_5days
    sample=$(echo "$R1" | sed -E 's/^([A-Za-z]+)-([0-9]+)_([0-9]+days).*/\1\2_\3/')

    echo "Processing: $sample | R1: $R1 | R2: $R2"

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 "$FASTQ_DIR/$R1" \
        -2 "$FASTQ_DIR/$R2" \
        -p 8 \
        --validateMappings \
        -o "$OUTPUT_DIR/$sample"
done

echo "Batch 4 processing complete."

#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=salmon_quant
#SBATCH --output=/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3/salmon_quant.log
#SBATCH --error=/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3/salmon_quant.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00

# Load Salmon
module load salmon/1.8.0

# Paths
INDEX_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_reference_genome/salmon_index_grch38_r97"
FASTQ_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/fastq_bulk"
OUTPUT_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3"

mkdir -p "$OUTPUT_DIR"

# Loop over R1 files
for R1 in ${FASTQ_DIR}/*_R1_001.fastq.gz; do
    # Find R2 mate
    R2="${R1/_R1_/_R2_}"

    # Clean sample name extraction (e.g., DMSO1_16hrs)
    sample=$(basename "$R1" | sed 's/-.*//')  # Keep everything before first hyphen
    sample="${sample}_16hrs"  # Append timepoint

    echo "Processing: $sample | R1: $(basename $R1) | R2: $(basename $R2)"

    # Salmon quantification
    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        -p 8 \
        --validateMappings \
        -o "${OUTPUT_DIR}/${sample}"

done

echo "Batch 3 processing complete."

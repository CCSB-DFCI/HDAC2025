#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=salmon_quant
#SBATCH --output=/n/scratch/users/o/old656/RNA_seq/yeast/salmon_quant.log
#SBATCH --error=/n/scratch/users/o/old656/RNA_seq/yeast/salmon_quant.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00

#Load Salmon module
module load salmon/1.8.0

#Define paths
FASTQ_DIR="/n/scratch/users/o/old656/RNA_seq/yeast/fastq_bulk"
INDEX_DIR="/n/scratch/users/o/old656/RNA_seq/yeast/reference_genome/yeast_salmon_index"
OUTPUT_DIR="/n/scratch/users/o/old656/RNA_seq/yeast/salmon_output"

mkdir -p ${OUTPUT_DIR}

#Loop through the paired-end FASTQ files
for sample in M208-1 M208-2 M208-3 ume6-1 ume6-2 ume6-3; do
    echo "Processing sample: ${sample}"

    #Set prefix exactly like M208 for ume6 samples
    prefix="JO-${sample}"

    #Define R1 and R2 file paths correctly
    R1_FILE=$(find ${FASTQ_DIR} -type f -name "${prefix}-NGS24-*_*_R1_001.fastq.gz" | head -n 1)
    R2_FILE=$(find ${FASTQ_DIR} -type f -name "${prefix}-NGS24-*_*_R2_001.fastq.gz" | head -n 1)

    echo "Looking for: ${FASTQ_DIR}/${prefix}-NGS24-*_*_R1_001.fastq.gz"
    echo "R1 found: ${R1_FILE}"
    echo "R2 found: ${R2_FILE}"

    #Check if both R1 and R2 exist
    if [[ -n "$R1_FILE" && -n "$R2_FILE" ]]; then
        salmon quant -i ${INDEX_DIR} \
            -l A \
            -1 ${R1_FILE} \
            -2 ${R2_FILE} \
            --validateMappings \
            --gcBias \
            -p 8 \
            -o ${OUTPUT_DIR}/${sample}_quant
    else
        echo "Missing FASTQ files for ${sample}. Skipping."
    fi
done

echo "Salmon quantification complete!"

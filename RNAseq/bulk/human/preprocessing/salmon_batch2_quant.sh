#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=salmon_quant_batch2
#SBATCH --output=/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2/salmon_quant.log
#SBATCH --error=/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2/salmon_quant.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00

module load salmon/1.8.0

INDEX_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_reference_genome/salmon_index_grch38_r97"
FASTQ_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/fastq_bulk"
OUTPUT_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2"

mkdir -p "$OUTPUT_DIR"
cd "$FASTQ_DIR"

for R1_full in *_R1_001.fastq.gz; do
    sample_base=$(echo "$R1_full" | sed -E 's/_L00[12]_R1_001.fastq.gz//')

    # Parse condition and replicate
    if [[ "$sample_base" =~ SKNBE2C-DMSO-([0-9]+)_ ]]; then
        rep="${BASH_REMATCH[1]}"
        sample="DMSO-${rep}"
    elif [[ "$sample_base" =~ SKNBE2C-E6R-([0-9]+)_ ]]; then
        rep="${BASH_REMATCH[1]}"
        sample="E6R-${rep}"
    elif [[ "$sample_base" =~ SKNBE2C-TSA-001-([0-9]+)_ ]]; then
        rep="${BASH_REMATCH[1]}"
        sample="TSA001_b2_${rep}"
    elif [[ "$sample_base" =~ SKNBE2C-TSA-1-([0-9]+)_ ]]; then
        rep="${BASH_REMATCH[1]}"
        sample="TSA1_b1_${rep}"
    else
        echo "Skipping unrecognized file: $R1_full"
        continue
    fi

    # Collect all lanes for R1 and R2
    R1_files=$(ls "${sample_base}"_L00[12]_R1_001.fastq.gz 2>/dev/null | sort)
    R2_files=$(ls "${sample_base}"_L00[12]_R2_001.fastq.gz 2>/dev/null | sort)

    if [ -z "$R1_files" ] || [ -z "$R2_files" ]; then
        echo "Error: Missing R1 or R2 files for $sample"
        continue
    fi

    echo "Processing: $sample"
    echo "R1 files: $R1_files"
    echo "R2 files: $R2_files"

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 <(zcat $R1_files) \
        -2 <(zcat $R2_files) \
        -p 8 \
        --validateMappings \
        -o "$OUTPUT_DIR/$sample"
done

echo "Batch 2 processing complete."

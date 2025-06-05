#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=salmon_quant_batch1
#SBATCH --output=/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1/salmon_quant.log
#SBATCH --error=/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1/salmon_quant.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00

module load salmon/1.8.0

INDEX_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_reference_genome/salmon_index_grch38_r97"
FASTQ_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/fastq_bulk"
OUTPUT_DIR="/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1"

mkdir -p "$OUTPUT_DIR"

cd "$FASTQ_DIR"

# Get unique sample names by stripping lane info (L001/L002) and R1 suffix
for R1_full in *_R1_001.fastq.gz; do
    # Extract base sample name (e.g., Batch1-Neurob-25-DMSO-1_NGS19-K520_BHJKHVDRXX_S25)
    sample_base=$(echo "$R1_full" | sed -E 's/_L00[12]_R1_001.fastq.gz//')
    
    # Extract clean sample ID (e.g., DMSO1_batch1)
    if [[ "$sample_base" =~ Batch1-Neurob-25-(DMSO|E6R|TSA)-([0-9]+)_ ]]; then
        condition="${BASH_REMATCH[1]}"
        rep="${BASH_REMATCH[2]}"
        sample="${condition}${rep}_batch1"
    elif [[ "$sample_base" =~ Neurob-1-(TSA)-([0-9]+)_ ]]; then
        condition="${BASH_REMATCH[1]}"
        rep="${BASH_REMATCH[2]}"
        sample="${condition}${rep}_batch1"
    else
        echo "Skipping unrecognized file: $R1_full"
        continue
    fi

    # Collect all lanes for this sample
    R1_files=$(ls "${sample_base}"_L00[12]_R1_001.fastq.gz 2>/dev/null)
    
    if [ -z "$R1_files" ]; then
        echo "Error: No R1 files found for $sample"
        continue
    fi

    echo "Processing: $sample"
    echo "Merging lanes: $R1_files"

    # Salmon quantification (single-end)
    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -r <(zcat $R1_files) \
        -p 8 \
        --validateMappings \
        -o "$OUTPUT_DIR/$sample"
done

echo "Batch 1 processing complete."

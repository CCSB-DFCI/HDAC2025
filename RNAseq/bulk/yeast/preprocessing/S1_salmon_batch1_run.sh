#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=salmon_quant
#SBATCH --output=/n/scratch/users/o/old656/yeast_31032025_Tina/salmon_output/salmon_quant.log
#SBATCH --error=/n/scratch/users/o/old656/yeast_31032025_Tina/salmon_output/salmon_quant.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00

# Load Salmon module
module load salmon/1.8.0

# Define paths
FASTQ_DIR="/n/scratch/users/o/old656/yeast_31032025_Tina/fastq_bulk/"
INDEX_DIR="/n/scratch/users/o/old656/RNA_seq/yeast/reference_genome/yeast_salmon_index"
OUTPUT_DIR="/n/scratch/users/o/old656/yeast_31032025_Tina/salmon_output"

mkdir -p ${OUTPUT_DIR}

# Conditions and specific replicate numbers
declare -A conditions_replicates=(
    ["DMSO"]="3 4 5"
    ["rpd3D"]="3 4 5"
    ["sin3D"]="2 3 4"
    ["TSA"]="1 4 5"
    ["E6R"]="1 3 5"
)

# Add the updated loop here
# Loop through conditions and their specific replicates
for condition in "${!conditions_replicates[@]}"; do
    reps=(${conditions_replicates[$condition]})
    for rep in "${reps[@]}"; do
        echo "Processing condition: ${condition}, replicate: ${rep}"

        # Handle specific exceptions for TSA and E6R naming patterns
        if [[ "$condition" == "TSA" ]]; then
            R1_L001=$(find ${FASTQ_DIR} -type f -name "Yeast-20-TSA-${rep}_NGS19-*_L001_R1_001.fastq.gz")
            R1_L002=$(find ${FASTQ_DIR} -type f -name "Yeast-20-TSA-${rep}_NGS19-*_L002_R1_001.fastq.gz")
        elif [[ "$condition" == "E6R" ]]; then
            R1_L001=$(find ${FASTQ_DIR} -type f -name "Yeast-50-E6R-${rep}_NGS19-*_L001_R1_001.fastq.gz")
            R1_L002=$(find ${FASTQ_DIR} -type f -name "Yeast-50-E6R-${rep}_NGS19-*_L002_R1_001.fastq.gz")
        else
            R1_L001=$(find ${FASTQ_DIR} -type f -name "Yeast-${condition}-${rep}_NGS19-*_L001_R1_001.fastq.gz")
            R1_L002=$(find ${FASTQ_DIR} -type f -name "Yeast-${condition}-${rep}_NGS19-*_L002_R1_001.fastq.gz")
        fi

        # Check if both lane files are present
        if [[ -n "$R1_L001" && -n "$R1_L002" ]]; then
            # Concatenate the two files into a temporary file
            CONCAT_FILE="${OUTPUT_DIR}/${condition}_${rep}_concat_R1.fastq.gz"
            cat $R1_L001 $R1_L002 > $CONCAT_FILE
            R1_FILE=$CONCAT_FILE
            echo "Concatenated lanes: $R1_L001 + $R1_L002 -> $CONCAT_FILE"
        elif [[ -n "$R1_L001" ]]; then
            # Use the single lane file if only one exists
            R1_FILE=$R1_L001
            echo "Using single lane file: $R1_FILE"
        else
            echo "Missing FASTQ files for ${condition} replicate ${rep}. Skipping."
            continue
        fi

        # Run Salmon quantification
        salmon quant -i ${INDEX_DIR} \
            -l A \
            -r ${R1_FILE} \
            --validateMappings \
            --gcBias \
            -p 8 \
            -o ${OUTPUT_DIR}/${condition}_${rep}_quant

        # Remove the temporary concatenated file
        if [[ -f $CONCAT_FILE ]]; then
            rm $CONCAT_FILE
        fi

        echo "Salmon quantification completed for ${condition} replicate ${rep}."
        done
done

echo "Salmon quantification complete!"

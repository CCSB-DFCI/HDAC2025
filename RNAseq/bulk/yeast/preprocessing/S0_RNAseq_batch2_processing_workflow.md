# Bulk RNA-seq Data Processing Workflow

**Date:** Feb 21, 2025  
**System:** HMS O2 SLURM cluster  
**Paths:**  
- Scratch: `/n/scratch/users/o/old656/`  
- Home: `/home/old656/`

---


#Step-1: Download compressed data (zipped files):
wget --no-check-certificate -O 20241118-GEN1973-Jolivet.gz https://filesender.belnet.be/download.php?token=e0b82ce6-367a-485b-ac42-6d4e4a4337fd&files_ids=3359943 

*Note:* `--no-check-certificate` bypasses SSL verification (use only when necessary for trusted sources).

#Step-2: Submit Unzip job to SLURM:
sbatch -p short --job-name=unzip_fastq --output=/n/scratch/users/o/old656/unzip_fastq.log --error=/n/scratch/users/o/old656/unzip_fastq.err --time=2:00:00 --mem=64G --cpus-per-task=4 --wrap="unzip -o -j /home/old656/HDAC_drugome2025/bulk_RNAseq/yeast_analysis/raw_data/20241118-GEN1973-Jolivet.gz -d /n/scratch/users/o/old656/"

#Step-3: Organize yeast & human SKNBE2C data- moving & organizing fastq/QC files: 
#Move FASTQ files to /RNA_seq/human_SKNBE2C/fastq_bulk/ 
mkdir -p /n/scratch/users/o/old656/RNA_seq/human_SKNBE2C/fastq_bulk
mv /n/scratch/users/o/old656/DMSO-*.fastq.gz /n/scratch/users/o/old656/RNA_seq/human_SKNBE2C/fastq_bulk/
mv /n/scratch/users/o/old656/E6R-*.fastq.gz /n/scratch/users/o/old656/RNA_seq/human_SKNBE2C/fastq_bulk/
mv /n/scratch/users/o/old656/TSA-*.fastq.gz /n/scratch/users/o/old656/RNA_seq/human_SKNBE2C/fastq_bulk/


#Move FastQC results (.html & .zip) to /RNA_seq/human_SKNBE2C/QC_bulk/ 
mkdir -p /n/scratch/users/o/old656/RNA_seq/human_SKNBE2C/QC_bulk
find /n/scratch/users/o/old656/ -type f \( -name "DMSO-*_fastqc.html" -o -name "DMSO-*_fastqc.zip" -o \
                                       -name "E6R-*_fastqc.html" -o -name "E6R-*_fastqc.zip" -o \
                                       -name "TSA-*_fastqc.html" -o -name "TSA-*_fastqc.zip" \) \
     -exec mv {} /n/scratch/users/o/old656/RNA_seq/human_SKNBE2C/QC_bulk/ \;


ls -lh /n/scratch/users/o/old656/RNA_seq/human_SKNBE2C/fastq_bulk/
ls -lh /n/scratch/users/o/old656/RNA_seq/human_SKNBE2C/QC_bulk/ 



#Yeast:
mkdir -p /n/scratch/users/o/old656/RNA_seq/yeast/fastq_bulk
mv /n/scratch/users/o/old656/JO-M208-*-NGS24-*.fastq.gz /n/scratch/users/o/old656/RNA_seq/yeast/fastq_bulk/
mv /n/scratch/users/o/old656/JO-ume6-*-NGS24-*.fastq.gz /n/scratch/users/o/old656/RNA_seq/yeast/fastq_bulk/

#Move FastQC results to /RNA_seq/yeast/QC_bulk/ using grep pattern 
mkdir -p /n/scratch/users/o/old656/RNA_seq/yeast/QC_bulk
find /n/scratch/users/o/old656/ -type f \( -name "JO-M208-*_fastqc.html" -o -name "JO-M208-*_fastqc.zip" -o -name "JO-ume6-*_fastqc.html" -o -name "JO-ume6-*_fastqc.zip" \) -exec mv {} /n/scratch/users/o/old656/RNA_seq/yeast/QC_bulk/ \;
 


# RNA-seq Preprocessing and Quantification of Human Neuroblastoma Cell Lines: Dose- and Time-Dependent Effects Across Batches 1–4 

This repository documents the standardized preprocessing & quantification pipeline applied to four batches of human SKNBE2C neuroblastoma RNA-seq data. Here, we investigated the transcriptomic responses to varying doses of DMSO, TSA, & E6R across multiple experimental batches, using both single-end and paired-end sequencing strategies, and including both short-term (16 hours) and long-term (120 hours) treatments. All analyses employ the latest Salmon quantification (version 1.8.0) methods & a consistent GRCh38 reference to enable robust, cross-batch comparisons of gene expression & differential response to treatment.

1. **Batch 1**: Single-End, Dual-Lane Sequencing

## Cell Line: SKNBE2C neuroblastoma cells after 16 hours treatment 
## Treatments (16 hrs timepoint):
- DMSO (16h, vehicle control): 3 replicates
- TSA (1 µM, 16h): 3 replicates
- E6R (25 µM, 16h): 3 replicates

Batch 1 samples were sequenced as single-end reads and split across two sequencing lanes (L001 & L002). For each biological sample, the corresponding FASTQ files from both lanes were concatenated on the fly to ensure that all reads were included in quantification. Only the R1 files were used, as this batch did not include paired-end data.

- FASTQ pattern:
*_L001_R1_001.fastq.gz and *_L002_R1_001.fastq.gz

- Processing:
Reads from both lanes were merged using zcat during Salmon quantification, with no intermediate files created.

- Sample naming:
Output directories were named using the condition and replicate (e.g., DMSO1_batch1, TSA3_batch1).


2. **Batch 2**: Paired-End, 2025 Sequencing

## Cell Line: SKNBE2C neuroblastoma cells 
## Treatments (16 hrs timepoint):
- DMSO (16h, vehicle control): 3 replicates
- TSA (0.01 µM, 16h): 3 replicates
- E6R (0.25 µM, 16h): 3 replicates
   

4. **Batches 3 & 4**: Paired-End, 2025 Sequencing

##  Human Batch 3 (New Data)
## Cell Line: SKNBE2C neuroblastoma cells
## Treatments (16 hrs timepoint):
- DMSO (vehicle control): 3 replicates
- TSA (0.1 µM): 3 replicates
- E6R (1 µM): 3 replicates

## Analysis:
- Processed with the updated pipeline (Salmon v1.8.0, GRCh38 reference) for consistency with Batches 1–2.
- Focuses on short-term treatment effects (16 h) for comparison with earlier batches.

## Human Batch 4 (New Data)
## Cell Line: SKNBE2C neuroblastoma cells
## Treatments (120 hrs timepoint):
- DMSO (vehicle control): 3 replicates
- TSA (0.01 µM): 3 replicates
- E6R (0.25 µM): 3 replicates

## Analysis:
- Processed with the same updated pipeline as all other batches.
- Linked to cell invasion assays initiated after 120 h of compound incubation.
- Will be analyzed independently from 16 h timepoints but compared later for time-dependent effects.

   
Batches 3 & 4 were generated in 2025 and sequenced as paired-end reads (R1 and R2). Each sample consisted of matched R1 and R2 FASTQ files. No lane merging was required. The data was processed using Salmon’s paired-end mode.

- FASTQ pattern:
*_R1_001.fastq.gz and *_R2_001.fastq.gz

- Processing:
Each sample’s R1 and R2 files were provided to Salmon for quantification.

- Sample naming:
Output directories included the condition, replicate, and timepoint (e.g., DMSO1_16hrs, E6R2_16hrs, TSA3_5days).

# Quantification: 
All batches were quantified using Salmon v1.8.0 with the same human transcriptome index (GRCh38). The --validateMappings option was used for improved accuracy. Each sample was processed in its own subdirectory, & logs were retained for reproducibility.

# Summary: 
- Batch 1: Single-end, dual-lane; merged on the lanes. 
- Batch 3 & 4: Paired-end, single-lane; R1/R2 matched per sample.
- Consistent reference & parameters used across all batches.

# Output: For each sample, a quant.sf file with transcript-level quantification, ready for downstream analysis.

## Output Structure

salmon_output_batchX/
├── [SAMPLE]/
│ ├── quant.sf # Transcript-level quantification
│ ├── logs/ # Salmon processing logs
│ └── aux/ # Auxiliary files
├── salmon_quant.log # SLURM job log
└── salmon_quant.err # SLURM error log


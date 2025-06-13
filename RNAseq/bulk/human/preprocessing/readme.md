# RNA-seq Preprocessing and Quantification of Human Neuroblastoma Cell Lines: Dose- and Time-Dependent Effects Across Batches 1–4

This repository documents the standardized preprocessing & quantification pipeline applied to four batches of human SKNBE2C neuroblastoma RNA-seq data. Here, we investigated the transcriptomic responses to varying doses of DMSO, TSA, & E6R across multiple experimental batches, using both single-end and paired-end sequencing strategies, and including both short-term (16 hours) and long-term (120 hours) treatments. All analyses employ the latest Salmon quantification (version 1.8.0) methods & a consistent GRCh38 reference to enable robust, cross-batch comparisons of gene expression & differential response to treatment.

## Batch 1: Single-End, Dual-Lane Sequencing

### Cell Line
SKNBE2C neuroblastoma cells after 16 hours of treatment. 

### Treatments (16 hrs timepoint)
- DMSO (16h, vehicle control): 3 replicates
- TSA (1 µM, 16h): 3 replicates
- E6R (25 µM, 16h): 3 replicates

### Processing Details
- **Sequencing:** Single-end, split across two lanes (L001 & L002).
- **Input FASTQ pattern:**  
  `*_L001_R1_001.fastq.gz` and `*_L002_R1_001.fastq.gz`
- **Preprocessing:**  
  For each sample, lane-specific FASTQ files were merged on-the-fly using `zcat` to concatenate reads from both L001 and L002, ensuring all reads were included.
- **Quantification:**  
  Salmon was run in single-end mode (`-r`) on the merged reads. No intermediate files were generated; merging occurred in-memory.
- **Sample naming:**  
  Output directories were named as `[Condition][Replicate]_batch1` (e.g., `DMSO1_batch1`, `TSA3_batch1`).

## Batch 2: Paired-End; Single-Lane Sequencing 

### Cell Line
SKNBE2C neuroblastoma cells after 16 hours of treatment.

### Treatments (16 hrs timepoint)
- DMSO (16h, vehicle control): 3 replicates
- TSA (0.01 µM, 16h): 3 replicates
- E6R (0.25 µM, 16h): 3 replicates

### Processing Details
- **Sequencing:** Paired-end, single-lane (L002).
- **Input FASTQ pattern:**  
  `*_R1_001.fastq.gz` and `*_R2_001.fastq.gz`
- **Preprocessing:**  
  No merging required. For each sample, R1 & R2 files were paired by sample name.
- **Quantification:**  
  Salmon was run in paired-end mode (`-1` and `-2` flags).
- **Sample naming:**  
  Output directories were named as `[Condition][Replicate]_batch2`.

## Batch 3: Paired-End, 2025 Sequencing

### Cell Line
SKNBE2C neuroblastoma cells

### Treatments (16 hrs timepoint)
- DMSO (vehicle control): 3 replicates
- TSA (0.1 µM): 3 replicates
- E6R (1 µM): 3 replicates

### Processing Details
- **Sequencing:** Paired-end, single-lane.
- **Input FASTQ pattern:**  
  `*_R1_001.fastq.gz` and `*_R2_001.fastq.gz`
- **Preprocessing:**  
  For each sample, R1 and R2 files were paired by sample name.
- **Quantification:**  
  Salmon was run in paired-end mode (`-1` and `-2` flags).
- **Sample naming:**  
  Output directories included the condition, replicate, and timepoint (e.g., `DMSO1_16hrs`, `E6R2_16hrs`, `TSA3_16hrs`).

## Batch 4: Paired-End, 2025 Sequencing

### Cell Line
SKNBE2C neuroblastoma cells

### Treatments (120 hrs timepoint)
- DMSO (vehicle control): 3 replicates
- TSA (0.01 µM): 3 replicates
- E6R (0.25 µM): 3 replicates

### Processing Details
- **Sequencing:** Paired-end, single-lane.
- **Input FASTQ pattern:**  
  `*_R1_001.fastq.gz` and `*_R2_001.fastq.gz`
- **Preprocessing:**  
  For each sample, R1 and R2 files were paired by sample name.
- **Quantification:**  
  Salmon was run in paired-end mode (`-1` and `-2` flags).
- **Sample naming:**  
  Output directories included the condition, replicate, and timepoint (e.g., `DMSO1_5days`, `E6R2_5days`, `TSA3_5days`).

- **Note:**  
  Batch 4 is linked to cell invasion assays initiated after 120 h of compound incubation. It will be analyzed independently from 16 h timepoints but compared later for time-dependent effects.

## Quantification Pipeline

- **Tool:** Salmon v1.8.0
- **Reference:** Human transcriptome index (GRCh38)
- **Parameters:**
  - `--validateMappings` for selective alignment
  - `-p 8` threads per sample
  - `-l A` for automatic library type detection
- **Batch-specific options:**
  - **Batch 1:** Single-end mode (`-r <(zcat ...)`)
  - **Batches 2–4:** Paired-end mode (`-1` and `-2`)
- **Logging:**  
  Each sample was processed in its own subdirectory; logs were retained for reproducibility.

## Output Structure
salmon_output_batchX/
├── [SAMPLE]/
│ ├── quant.sf # Transcript-level quantification
│ ├── logs/ # Salmon processing logs
│ └── aux/ # Auxiliary files
├── salmon_quant.log # SLURM job log
└── salmon_quant.err # SLURM error log


## Summary

- **Batch 1:** Single-end, dual-lane; merged on the lanes.
- **Batch 2, 3 & 4:** Paired-end, single-lane; R1/R2 matched per sample.
- **Consistent reference & parameters** used across all batches.
- For each sample, a `quant.sf` file with transcript-level quantification, ready for downstream analysis.



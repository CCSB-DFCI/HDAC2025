Differential Expression Analysis with DESeq2: 
DESeq2 is used for identifying differentially expressed genes or transcripts from our neuroblastoma RNA-seq data. It models raw count data using negative binomial generalized linear models, estimates variance and dispersion across biological replicates, and corrects for differences in sequencing depth to provide robust statistical testing of expression changes between experimental conditions.

In this project, we:
- Quantified transcript abundances using Salmon.
- Performed batch-specific DESeq2 analyses, comparing each treatment group (TSA, E6R) to its respective DMSO control within each batch.
- Pre-filtered low-count features, normalized data & applied shrinkage estimation for log2 fold changes.
- Identified differentially expressed transcripts using an adjusted p-value (FDR) cutoff of 0.05 & an absolute log2 fold change cutoff of ≥1.
- Saved results & diagnostic plots in batch-specific output folders for downstream interpretation and reproducibility.

This approach ensures accurate detection of biologically meaningful expression changes while accounting for technical & biological variability inherent in RNA-seq experiments.

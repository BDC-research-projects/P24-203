# P24-203

Integration and analysis of publicly available single-nucleus RNA sequencing (snRNA-seq) datasets to investigate the expression of **GFRAL** and other genes within the brain's dorsal vagal complex (DVC).

**Zenodo:** https://zenodo.org/records/20040970

## Overview

This repository contains the code, metadata, and analysis workflows used to download, preprocess, integrate, and analyse publicly available snRNA-seq datasets from the mouse dorsal vagal complex. The workflow produces integrated Cell Ranger count matrices, performs downstream analyses in R, and generates the figures and tables used in the accompanying study.

## Workflow

The complete analysis consists of four major steps:

1. **Download raw sequencing data**

   Raw FASTQ files are downloaded from the NCBI Sequence Read Archive (SRA) using:

   ```bash
   scripts/run_fasterq_dump.sh
   ```

2. **Pre-processing**

   Raw sequencing data are processed using:

   - **nf-core/scrnaseq** v2.7.1
   - **Cell Ranger** v8.0.0

   This step generates filtered feature-barcode matrices for each sample.

3. **Data integration**

   Integrated analysis, quality control, clustering, annotation, and visualization are performed using:

   ```text
   scripts/snRNA-seq_integration.Rmd
   ```

4. **Co-expression analysis**

   Gene co-expression analyses for GFRAL and genes of interest are performed using:

   ```text
   scripts/snRNA_coexpression.Rmd
   ```

## Software

| Software | Version |
|----------|---------|
| SRA Toolkit | 3.1.0 |
| nf-core/scrnaseq | 2.7.1 |
| Cell Ranger | 8.0.0 |
| R | ≥4.4 (recommended) |

## Reproducibility

Detailed information required to reproduce the complete workflow, including software versions, reference files, and analysis parameters, is provided in **REPRODUCIBILITY.md**.

## Citation

If you use this repository, please cite:

- This GitHub repository
- The archived Zenodo release (DOI above)
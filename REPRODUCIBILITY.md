# Reproducibility guide for P24-203

This repository analyses publicly available single-nucleus RNA-seq datasets from the dorsal vagal complex to investigate expression of GFRAL and other genes of interest.

## Overview of Workflow

The analysis was performed in the following order:

1. Download raw sequencing data  
   → `scripts/run_fasterq_dump.sh`

2. Pre-process data (FASTQ → count matrices)  
   → nf-core/scrnaseq (v2.7.1) with Cell Ranger (v8.0.0)

3. Perform QC, integration, clustering, annotation, and visualisation  
   → `scripts/snRNA-seq_integration.Rmd`

4. Perform co-expression analysis, visualisation and statistical tests  
   → `scripts/snRNA_coexpression.Rmd`

## Installation guide 

### Required input data

See `data_manifest.tsv` for dataset accessions and expected files.

### Clone the repository 

```
git clone https://github.com/BDC-research-projects/P24-203.git
```

### Expected directory structure

```
P24-203/
├── data/
│   ├── fastq/
│   ├── count/
│   ├── csv/
│   │	├── neurons/
│   │	├── non_neurons/  
│
├── doc/
├── img/
└── scripts/
```

# Step 1: Download raw sequencing data

```
bash scripts/run_fasterq_dump_parallel.sh <SRR_FILE> <OUTPUT_DIR> <NUM_THREADS> <SINGULARITY>
```

SRR files can be found in `doc` and the singularity container can be built using command

```
apptainer build sratoolkit_3.1.0.sif docker://pegi3s/sratoolkit:3.1.0
``` 

# Step 2: Pre-processing using nf-core/scrnaseq

## Example samplesheet
```
sample,fastq_1,fastq_2,expected_cells
6V,data/fastq/SRR13694200_R1.fastq.gz,data/fastq/SRR13694200_R2.fastq.gz,4500
6V,data/fastq/SRR13694201_R1.fastq.gz,data/fastq/SRR13694201_R2.fastq.gz,4500
6V,data/fastq/SRR13694202_R1.fastq.gz,data/fastq/SRR13694202_R2.fastq.gz,4500
6V,data/fastq/SRR13694203_R1.fastq.gz,data/fastq/SRR13694203_R2.fastq.gz,4500
9V,data/fastq/SRR13694204_R1.fastq.gz,data/fastq/SRR13694204_R2.fastq.gz,4500
9V,data/fastq/SRR13694205_R1.fastq.gz,data/fastq/SRR13694205_R2.fastq.gz,4500
14V,data/fastq/SRR13694206_R1.fastq.gz,data/fastq/SRR13694206_R2.fastq.gz,4500
14V,data/fastq/SRR13694207_R1.fastq.gz,data/fastq/SRR13694207_R2.fastq.gz,4500
...
```
## Example nf-params.json 

```
{
    "input": "samplesheet.csv",
    "outdir": "results",
    "email": "your.email@org.com",
    "fasta": "reference/GRCm39.primary_assembly.genome.fa.gz",
    "gtf": "reference/gencode.vM28.annotation.gtf.gz",
    "aligner": "cellranger",
    "protocol": "10XV2"
}
```

The reference files can be downloaded from https://www.gencodegenes.org/mouse/release_M28.html

## Example command

```
nextflow run nf-core/scrnaseq -r 2.7.1 -profile your_profile,singularity -resume -params-file nf-params.json -c your.config
```

# Step 3: snRNA-seq integration and downstream analysis

Run `scripts/snRNA_integration.Rmd`

## This script performs:

- Quality control and filtering
- Ambient RNA removal (e.g. SoupX)
- Doublet detection (e.g. DoubletFinder)
- Normalisation and scaling
- Dataset integration
- Dimensionality reduction (e.g. PCA, UMAP)
- Clustering
- Cell-type annotation
- Visualisation

### Inputs
- Count matrices from data/count/
- Metadata and reference files from data/csv/
- AP/NTS marker genes from doc/

### Outputs
- Data objects → data/
- Non-neuronal subset outputs → data/marker_genes/
- Neuronal subset outputs → data/neuron_marker_genes/
- Figures (UMAPs, QC plots, etc.) → img/

# Step 4: Co-expression analysis and visualisation

Run `scripts/snRNA_plot_coexpression.Rmd`

## This script performs:

- Selection of genes of interest (GOIs)
- Co-expression network construction (e.g. scLink)
- Network analysis
- Statistical testing
- Visualisation of co-expression patterns

### Inputs
- Processed data generated in Step 3

### Outputs
- Figures → img/

# Notes on reproducibility
- The analysis is documented without modifying the original code (with a few minor exceptions, e.g. modified output directory names)
- Some variability may occur due to stochastic steps (e.g. clustering, UMAP)
- Setting random seeds is recommended where possible

# Software environment

## Tool and package versions
- SRA Toolkit: 3.1.0
- CELLRANGER: 8.0.0
- fastqc: 0.12.1
- gunzip: 1.1
- nf-core/scrnaseq: v2.7.1-g4171377
- Nextflow: 24.10.2
- R: 4.4.0
- Bioconductor: 3.19
- Seurat: 5.1.0
- SeuratObject: 5.0.2
- SoupX: 1.6.2
- scran: 1.32.0
- scuttle: 1.14.0
- Polychrome: 1.5.1
- readr: 2.1.5
- scLink: 1.0.1
- SAVER: 1.1.2
- reshape2: 1.4.4
- rlang: 1.1.4
- ggplot2: 3.5.1
- harmony: 1.2.2
- patchwork: 1.2.0
- cowplot: 1.1.3
- foreach: 1.5.2
- doParallel: 1.0.17
- GSEABase: 1.66.0
- AUCell: 1.26.0
- here: 1.0.1

## Operating system
- NAME="Red Hat Enterprise Linux"
- VERSION="8.10 (Ootpa)"
- ID="rhel"
- ID_LIKE="fedora"

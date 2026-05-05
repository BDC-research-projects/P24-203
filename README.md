# P24-203

Integratation and analysis of publicly available sn-RNA sequencing datasets to investigate the expression of the GFRAL and other genes in the brain's dorsal vagal complex region. 

## Workflow

The analysis consists of:

1. Downloading raw data using `scripts/run_fasterq_dump.sh`
2. Pre-processing using [nf-core/scrnaseq](https://nf-co.re/scrnaseq/2.7.1/docs/output/#warning-please-read-this-documentation-on-the-nf-core-website-httpsnf-corescrnaseqoutput) (v2.7.1) with Cell Ranger (v8.0.0)
3. Running `scripts/snRNA-seq_integration.Rmd`
4. Running `scripts/snRNA_coexpression.Rmd`

For full reproducibility details, see `REPRODUCIBILITY.md`.

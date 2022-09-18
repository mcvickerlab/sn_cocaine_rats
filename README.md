# sn_cocaine_rats

This repository contains relevant code/pipelines for work done in our paper ["Cocaine addiction-like behaviors are associated with long-term changes in gene regulation, energy metabolism, and GABAergic inhibition within the amygdala"](https://www.biorxiv.org/content/10.1101/2022.09.08.506493v1).

## [`rna_analysis`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis)
This directory contains code/pipelines for analysis of snRNA-seq data used in this project. 

### [`make_seurat_rna`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/make_seurat_rna)
This pipeline runs [`cellranger-count`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count) on FASTQ files for each sample and creates Seurat objects for each sample along with QC plots.

### ['rpca'](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/dge_labeled)
This pipeline takes the output of the [`make_seurat_rna`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/make_seurat_rna) pipline and subsets the cells based on QC metrics; normalizes the data; and integrates the data from all samples using [RPCA](https://satijalab.org/seurat/articles/integration_rpca.html). Cells are then clustered and marker genes of clusters are identified. All of this is performed with Seurat. Outputs of this pipeline are used for manual cell type annotation. 

### ['dge_labeled'](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/dge_labeled)
This pipeline takes **annotated** snRNA-seq data and performs differential gene expression analysis between high and low AI rats (as described in our paper) using MAST in Seurat. The `config.yml` file can be modified to change the parameters of the test for differential expression. 



## [`atac_analysis`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis)
This directory contains code/pipelines for analysis of snATAC-seq data used in this project.

# sn_cocaine_rats

This repository contains relevant code/pipelines for work done in our paper ["Cocaine addiction-like behaviors are associated with long-term changes in gene regulation, energy metabolism, and GABAergic inhibition within the amygdala"](https://www.biorxiv.org/content/10.1101/2022.09.08.506493v1) (bioRxiv).

## [`rna_analysis`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis)
This directory contains code/pipelines for analysis of snRNA-seq data used in this project. 

### [`make_seurat_rna`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/make_seurat_rna)
This pipeline runs [`cellranger-count`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count) on FASTQ files for each sample to align them to an rn6 reference genome and creates Seurat objects for each sample along with QC plots.

### [`rpca`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/dge_labeled)
This pipeline takes the output of the [`make_seurat_rna`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/make_seurat_rna) pipline and subsets the cells based on QC metrics; normalizes the data; and integrates the data from all samples using [RPCA](https://satijalab.org/seurat/articles/integration_rpca.html). Cells are then clustered and marker genes of clusters are identified. All of this is performed with Seurat. Outputs of this pipeline are used for manual cell type annotation. 

### [`dge_labeled`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/dge_labeled)
This pipeline takes **annotated** snRNA-seq data and performs differential gene expression analysis between high and low AI rats (as described in our paper) using MAST in Seurat. The `config.yml` file can be modified to change the parameters of the test for differential expression. 

### [`clusterProfiler`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/rna_analysis/clusterProfiler/scripts)
This directory contains scripts used to run clusterProfiler on the annotated snRNA-seq data to perform cell type-specific GSEA. 

## [`atac_analysis`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis)
This directory contains code/pipelines for analysis of snATAC-seq data used in this project.

### [`filter_rn6.sh`](https://github.com/mcvickerlab/sn_cocaine_rats/blob/master/atac_analysis/filter_rn6.sh)
This script contains code used to create a custom rn6 reference genome.

### [`realign_snatac_fastqs`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis/realign_snatac_fastqs)
This directory contains scripts used to realign snATAC-seq FASTQs to the custom rn6 reference genome using [`cellranger-atac count`](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count).

### [`snatac_integration`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis/snatac_integration)
This pipeline takes the [`cellranger-atac count`](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count) outputs from the previous step and generates a combined peak matrix which is then used to generate new peak by cell matrices for each sample. These are then loaded into Signac and preprocessed accordingly to generate an integrated dataset across all samples (using [Harmony](https://portals.broadinstitute.org/harmony/articles/quickstart.html)). This pipeline also requires an annotated integrated snRNA-seq Seurat object for label transfer (for annotation of cell types in the snATAC-seq dataset).

### [`dapeak_analysis_rn6`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis/dapeak_analysis_rn6)
This pipeline takes annotated snATAC-seq data as a Signac object and performs cell type-specific peak accessibility analysis between high and low AI rats. It also performs a permutation test and pulls out foreground/background sequences based on significance values, which may be useful for downstream enrichment analyses (e.g. motif finding).

### [`coembed`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis)
This directory contains scripts used to coembed the snRNA-seq with the snATAC-seq.

### [`dapeaks_partition_h2`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis/dapeaks_partition_h2)
This pipeline contains code used for performing a partitioned heritability analysis for enrichment of GWAS risk variants in cell type-specific open chromatin regions.

### [`call_peaks_macs`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis/call_peaks_macs)
This pipeline uses MACS to call peaks on the BAM outputs of `cellranger-atac count` in a cell type-specific manner. Requires a mapping of cell barcodes to cell type annotations. May be useful for downstream analyses. 

### [`chromvar_analysis`](https://github.com/mcvickerlab/sn_cocaine_rats/tree/master/atac_analysis/chromvar_analysis)
This directory contains code used to perform differential analysis of transcription factors. Notebooks include pipelines for running chromVar analysis, as well as data visualization.


**Questions? Contact jlz014 [at] eng [dot] ucsd [dot] edu**

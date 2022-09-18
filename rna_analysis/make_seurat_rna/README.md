# make_seurat

This pipeline runs [cellranger-count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) with a custom reference genome on FASTQ files for scRNA-seq experiments and generates Seurat objects (`.rds`) from the resulting count data with corresponding QC plots. 

This pipeline was designed for analyzing snRNA-seq from outbred rats classified by addiction index following self-administration of oxycodone or cocaine, in addition to rats that were never exposed to any drug (naive). 

## Usage
To run this pipeline, fill out the fields in the [config file](https://github.com/zrcjessica/make_seurat/blob/main/config.yml) with your respective paths to the specified files and directories. 

### [`config.yml`](https://github.com/zrcjessica/make_seurat/blob/main/config.yml) 

#### `out`
Specify the relative path to the directory where output files will be written

#### `data`
Specify path where you want outputs of `cellranger count` to be written - can specify the same directory as `out`

#### `fastqs_ref`
This specifies a .tsv file with four columns:
| RFID                              | Treatment                      | Addiction index      | FASTQ prefix                                      |
|-----------------------------------|--------------------------------|----------------------|---------------------------------------------------|
| Should correspond to sample names | e.g. cocaine, oxycodone, naive | e.g. high, low, none | Sample prefix for all FASTQ files for this sample; see [10x docs](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input) for details |

The tsv file should keep the row of header names. Header names can be changed.

#### `fastqs_dir`
Path to directory containing subdirectories each named after a sample (RFID). Each subdirectory contains all the FASTQ files for each sample. 

e.g. 
```bash
fastqs_dir
├── SampleA
│   ├── SampleA_1_S41_L002_R1_001.fastq.gz
│   ├── SampleA_1_S41_L002_R2_001.fastq.gz
│   ├── SampleA_2_S42_L002_R1_001.fastq.gz
│   ├── SampleA_2_S42_L002_R2_001.fastq.gz
│   ├── SampleA_3_S43_L002_R1_001.fastq.gz
│   ├── SampleA_3_S43_L002_R2_001.fastq.gz
│   ├── SampleA_4_S44_L002_R1_001.fastq.gz
│   ├── SampleA_4_S44_L002_R2_001.fastq.gz
│   ├── SampleA_S2_L001_R1_001.fastq.gz
│   └── SampleA_S2_L001_R2_001.fastq.gz
├── SampleB
│   ├── SampleB_S4_L002_R1_001.fastq.gz
│   └── SampleB_S4_L002_R2_001.fastq.gz
├── SampleC
│   ├── SampleC_S3_L002_R1_001.fastq.gz
│   └── SampleC_S3_L002_R2_001.fastq.gz
```
#### `ref_transcriptome`
Path to 10x custom reference genome built with [cellranger mkref](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr) for your organism. 

#### cellranger count input params
Can specify `expect_cells`, `chemistry`, `localmem`, and `localcores`. See 10x [Command-Line Argument Reference](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#args) for details. 

#### Seurat `CreateSeuratObject()` arguments
Can specify `min.cells` and `min.features`. See docs [here](https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/CreateSeuratObject). 

#### `samples`
List of sample names corresponding to first column of tsv file referenced by `fastqs_ref`.


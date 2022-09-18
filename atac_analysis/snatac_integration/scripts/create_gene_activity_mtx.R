library(argparse)
suppressMessages(library(Seurat))
suppressMessages(library(Signac))

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rds", action = "store",
	help = "path to rds object to compute gene activity for")
parser$add_argument("--out", action = "store",
	help = "where to save rds object with gene activities")

args <- parser$parse_args()

dat <- readRDS(args$rds)

gene.activities <- GeneActivity(dat)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
dat[['RNA']] <- CreateAssayObject(counts = gene.activities)
dat <- NormalizeData(
  object = dat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(dat$nCount_RNA)
)

DefaultAssay(dat) <- "RNA"

saveRDS(dat, args$out)
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
library(ggplot2)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--atac', action = "store", type = "character", 
	help = "path to snATAC-seq rds file")
parser$add_argument("--rna", action = "store", type = "character",
	help = "path to snRNA-seq rds file")
parser$add_argument("--out", action = "store", type = 'character',
	help = "path to save RDS file with predicted.id")
# parser$add_argument("--downsample", action = "store", type = "character", nargs= = "+", default = NULL,
# 	help = "if provided, downsample scRNA-seq dataset to this many samples")

args <- parser$parse_args()

print("loading RDS files")
atac <- readRDS(args$atac)
rna <- readRDS(args$rna)

print("setting default assay")
DefaultAssay(rna) <- "integrated"
DefaultAssay(atac) <- "RNA"

rna$celltype <- Idents(rna)

print("finding transfer anchors")
transfer.anchors <- FindTransferAnchors(
  reference = rna, 
  query = atac,
  reduction = 'cca'
)

print("transferring anchors")
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna$celltype,
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

print("adding predicted labels")
atac <- AddMetaData(object = atac, metadata = predicted.labels)

print("saving RDS")
saveRDS(atac, file = args$out)

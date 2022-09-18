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
parser$add_argument("--activity", action = "store_true", 
  help = "if true, GeneActivity assay already exists under name RNA")

args <- parser$parse_args()

print("loading RDS files")
atac <- readRDS(args$atac)
rna <- readRDS(args$rna)

if (args$activity) {
  print("renaming existing RNA assay to GeneActivity")
  atac <- RenameAssays(atac, RNA = "GeneActivity")
}

print("setting default assay")
DefaultAssay(rna) <- "integrated"
DefaultAssay(atac) <- "GeneActivity"

rna$celltype <- Idents(rna)

print("finding transfer anchors")
transfer.anchors <- FindTransferAnchors(
  reference = rna, 
  query = atac,
  reduction = 'cca'
)

print("imputing RNA expression")
imputed.rna <- TransferData(
  anchorset = transfer.anchors,
  refdata = GetAssayData(rna, assay = "RNA", slot = "data"),
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

print("adding imputed expression data")
atac[["RNA"]] <- imputed.rna

print("saving RDS")
saveRDS(atac, file = args$out)

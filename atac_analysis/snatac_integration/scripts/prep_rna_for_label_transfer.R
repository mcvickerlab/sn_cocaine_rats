suppressMessages(library(Seurat))
suppressMessages(library(Signac))
library(ggplot2)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rna", action = "store", type = "character",
	help = "path to snRNA-seq rds file")
parser$add_argument("--out", action = "store", type = 'character',
	help = "path to save RDS file with predicted.id")
args <- parser$parse_args()

dat <- readRDS(args$rna)

rats.list <- SplitObject(dat, split.by = "sample")

# normalize and identify variable features for each dataset independently
rats.list <- lapply(X = rats.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- ScaleData(x, vars.to.regress = "percent.mt")
})

features <- SelectIntegrationFeatures(object.list = rats.list)

anchors <- FindIntegrationAnchors(object.list = rats.list, anchor.features = features)

combined <- IntegrateData(anchorset = anchors)

combined <- ScaleData(combined)

saveRDS(combined, args$out)
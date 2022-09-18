library(Seurat)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--rds', action = "store", type = "character", 
	help = "path to rds file for integrated object")
parser$add_argument("--test", action = "store", type = "character",
	help = "test to use for FindAllMarkers")
parser$add_argument("--out", action = "store", type = "character",
	help = "where to save table of cluster markers")
parser$add_argument("--prefilter", action = "store_true",
	help = "if true, prefilter genes with Seurat defaults")
parser$add_argument("--covars", action = "store", nargs = "+", type = "character",
	help = "list of covars for model")

args <- parser$parse_args()

data <- readRDS(args$rds)
# DefaultAssay(data) <- "SCT"
DefaultAssay(data) <- "RNA"
data <- NormalizeData(data)
data <- ScaleData(data)

if (args$prefilter) {
	markers <- FindAllMarkers(data, test.use = args$test, 
		latent.vars = args$covars)
	write.csv(markers, file = args$out, quote = FALSE, row.names = FALSE)
} else {
	markers <- FindAllMarkers(data, test.use = args$test, 
		logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, 
		min.cells.group=0, latent.vars = args$covars)
	write.csv(markers, file = args$out, quote = FALSE, row.names = FALSE)
}


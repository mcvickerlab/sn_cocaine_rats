library(Seurat)
library(argparse)
library(patchwork)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--rds', action = "store", type = "character", 
	help = "path to rds file for integrated object")
parser$add_argument("--marker_genes", action = "store", type = "character",
	nargs = "+", help = "path to tsv file with cell type marker genes")
parser$add_argument("--featureplot_fh", action = "store", type = "character",
	help = "path to save feature plot")
parser$add_argument("--dotplot_fh", action = "store", type = "character",
	help = "path to save dotplot")

args <- parser$parse_args()

data <- readRDS(args$rds)

DefaultAssay(data) <- "RNA"
data <- NormalizeData(data)

if (sum(args$marker_genes %in% rownames(data))>0) {
	fp <- FeaturePlot(data, features = args$marker_genes)

	png(args$featureplot_fh, units = "in", width = 12, height = 12, res = 300)
	print(fp)
	dev.off()

	dp <- DotPlot(data, features = args$marker_genes)

	png(args$dotplot_fh, units = "in", width = 12, height = 12, res = 300)
	print(dp)
	dev.off()

} else {
	print("None of the requested features were found for celltype")
}


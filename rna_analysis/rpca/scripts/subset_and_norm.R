library(Seurat)
library(argparse)
library(sctransform)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--rds', action = "store", type = "character", 
	help = "path to rds file")
parser$add_argument('--subset', action = "store_true",
	help = "if true, filter out low-quality cells")
parser$add_argument("--sctransform", action = "store_true",
	help = "if true, use sctransform in place of normalization + scaling steps")
parser$add_argument("--out", action = "store", type = "character",
	help = "where to write output file (rds)")

args <- parser$parse_args()

data <- readRDS(args$rds)

# subset
# if (!is.null(args$percent_mt) & !is.null(args$nFeature_RNA)) {
# 	cat("subsetting data\n")
# 	data.nFeature_RNA <- FetchData(object = data, vars = 'nFeature_RNA')
# 	data.percent_mt <- FetchData(object = data, vars = 'percent.mt')
# 	cat(sprintf("keeping cells where %d < nFeature_RNA < %d and %.1f < percent.mt < %.1f\n", 
# 		args$nFeature_RNA[1], args$nFeature_RNA[2], args$percent_mt[1], args$percent_mt[2]))
# 	data <- data[, which(x = data.nFeature_RNA > args$nFeature_RNA[1] & 
# 		data.nFeature_RNA < args$nFeature_RNA[2] & 
# 		data.percent_mt > args$percent_mt[1] &
# 		data.percent_mt < args$percent_mt[2])]
# }

if (args$subset) {
	cat("subsetting data\n")
	features.data.mean <- mean(data@meta.data$nFeature_RNA)
	features.data.ds <- sd(data@meta.data$nFeature_RNA)

	mito.data.mean <- mean(data@meta.data$percent.mt)
	mito.data.ds <- sd(data@meta.data$percent.mt)

	data <- subset(data, subset = nFeature_RNA>(features.data.mean-3*features.data.ds) &
		nFeature_RNA<(features.data.mean + 3*features.data.ds) &
		percent.mt < mito.data.mean + 3*mito.data.ds)
}

if (args$sctransform) {
	cat("normalizing data with sctransform\n")
	data <- SCTransform(data, vars.to.regress = "percent.mt",
		verbose = FALSE, return.only.var.genes = FALSE)
	saveRDS(data, file = args$out)
} else {
	cat("normalizing data\n")
	data <- NormalizeData(data)
	cat("finding variable features\n")
	data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

	# # write variable features to table
	# cat("writing variable features to table\n")
	# write.csv(HVFInfo(object = data), 
	# 	file = file.path(var6, 'variable_features', 
	# 		paste(sample, 'variable_features', sprintf('%d_varfeatures.csv', var5), sep = '_')), 
	# 	quote = F)

	saveRDS(data, file = args$out)
}
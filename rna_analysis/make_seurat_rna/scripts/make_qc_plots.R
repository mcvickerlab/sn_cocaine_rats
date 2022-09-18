library(Seurat)
library(optparse)
library(ggplot2)

# load arguments
option_list <- list(
	make_option("--input_dir", action = "store", type = "character", 
		help = "path to directory where cell ranger filtered feature barcode matrix is located"),
	make_option("--sample", action = "store", type = "character",
		help = "sample name"),
	make_option("--treatment", action = "store", type = "character",
		help = "drug treatment (cocaine, oxycodone, naive)"),
	make_option("--addiction_index", action = "store", type = "character",
		help = "addiction_index (high, low, none)"),
	make_option("--min_cells", action = "store", default = 3,
		type = "integer", help = "value for min.cells argument of CreateSeuratObject(); default = 3"),
	make_option("--min_features", action = "store", default = 200,
		type = "integer", help = "value for min.features argument of CreateSeuratObject(), default = 200"),
	make_option("--vln_plot", action = "store", type = "character",
		help = "full path to file to save QC violin plots"),
	make_option("--rds", action = "store", type = "character",
		help = "full path to file to save .rds"))

plotViolin <- function(data, out, sample) {
	# Visualize QC metrics as a violin plot
	qcmetrics <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	png(out, units = "in", width = 12, height = 8, res = 300)
	print(qcmetrics)
	dev.off()
}

opt <- parse_args(OptionParser(option_list=option_list))
print('args parsed')

cat(sprintf("processing sample %s\n", opt$sample))
data <- Read10X(opt$input_dir)
data <- CreateSeuratObject(counts = data, 
	project = opt$sample, 
	min.cells = opt$min_cells, 
	min.features = opt$min_features)

cat("Seurat object created\n")

data[["sample"]] <- opt$sample
data[["treatment"]] <- opt$treatment
data[["addiction.index"]] <- opt$addiction_index
if (opt$treatment == "naive") {
	data[["label"]] <- "naive"
} else {
	data[["label"]] <- sprintf("%s_%s", opt$treatment, opt$addiction_index)
}
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^Mt-")

plotViolin(data, opt$vln_plot, opt$sample)
saveRDS(data, file = opt$rds)


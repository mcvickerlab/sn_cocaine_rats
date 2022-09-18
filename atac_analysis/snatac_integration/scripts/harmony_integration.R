library(argparse)
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
library(ggplot2)
suppressMessages(library(harmony))

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--rds', action = "store", type = "character", 
	nargs = "+", help = "paths to rds files to integrate")
parser$add_argument("--merged_out", action = "store", type = "character",
	help = "where to save RDS for merged and unintegrated data")
parser$add_argument("--out", action = "store", type = "character",
	help = "where to write output file (rds)")
parser$add_argument("--umap", action = "store", type = "character",
	help = "path to save UMAP (png)")
parser$add_argument("--conditions_umap", action = "store", type = "character",
	help = "path to save UMAP split by conditions (png)")
parser$add_argument("--samples_umap", action = "store", type = "character",
	help = "path to save UMAP split by samples")
parser$add_argument("--reference", action = "store", default = NULL,
	help = "name of sample to use as reference for merging")
parser$add_argument("--vars", action = "store", nargs = "+", type = "character",
	help = "variables to integrate on with Harmony")

args <- parser$parse_args()

print(args$vars)

# load samples
samples.list <- lapply(args$rds, function(x) {
	print(x)
	readRDS(x)
	})

sample.names <- lapply(samples.list, function(x) {
	return(x$sample[1])
	})

names(samples.list) <- sample.names

print(samples.list)
# merge samples
combined <- NULL
if (!is.null(args$reference)) {
	print("reference given")
	ref <- samples.list[[args$reference]]
	combined <- merge(x = ref,
		y = samples.list[sample.names != args$reference],
		add.cell.ids = c(args$reference, sample.names[sample.names != args$reference]))
} else {
	print("no reference given")
	combined <- merge(x = samples.list[[1]],
		y = samples.list[2:length(samples.list)],
		add.cell.ids = sample.names)
}

# compute LSI for merged object
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
# combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
# combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)
saveRDS(combined , args$merged_out)

# integrated with harmony
hm.integrated <- RunHarmony(
  object = combined,
  group.by.vars = args$vars,
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
# hm.integrated <- FindNeighbors(object = hm.integrated, reduction = 'lsi', dims = 2:30)
hm.integrated <- FindNeighbors(object = hm.integrated, reduction = 'harmony', dims = 2:30)
hm.integrated <- FindClusters(object = hm.integrated, verbose = FALSE, algorithm = 3)

# save RDS
saveRDS(hm.integrated, args$out)

# plot harmony UMAP
harmony.umap <- DimPlot(hm.integrated, pt.size = 0.1) + 
	ggtitle("Harmony integration") + NoAxes() +  
	theme(legend.text = element_text(size = 20),
		title = element_text(size = 26))

# save UMAP
png(args$umap, units = "in", width = 16, height = 8, res = 300)
print(harmony.umap)
dev.off()

# plot UMAP split by condition
conditions.umap <- DimPlot(hm.integrated, group.by = 'condition', pt.size = 0.1) + 
	ggtitle("Harmony integration") + NoAxes() +  
	theme(legend.text = element_text(size = 20),
		title = element_text(size = 26))

png(args$conditions_umap, units = "in", width = 16, height = 8, res = 300)
conditions.umap
dev.off()

# plot UMAP split by samples
samples.umap <- DimPlot(hm.integrated, group.by = 'sample', pt.size = 0.1) + 
	ggtitle("Harmony integration") + NoAxes() +  
	theme(legend.text = element_text(size = 20),
		title = element_text(size = 26))

png(args$samples_umap, units = "in", width = 16, height = 8, res = 300)
samples.umap
dev.off()

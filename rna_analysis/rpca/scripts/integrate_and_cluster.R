library(Seurat)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--rds', action = "store", type = "character", 
	nargs = "+", help = "paths to rds files to integrate")
parser$add_argument('--sample', action = "store", type = "character",
	nargs = "+", help = "sample names corresponding to rds files")
parser$add_argument("--out_rds", action = "store", type = "character",
	help = "where to write output file (rds)")
parser$add_argument("--clusters_umap", action = "store", type = "character",
	help = "where to save UMAP plot colored by cluster (png)")
parser$add_argument("--conditions_umap", action = "store", type = "character",
	help = "where to save UMAP plot colored by sample (png)")
parser$add_argument("--split_umap", action = "store", type = "character",
	help = "where to save side by side UMAP plot")
parser$add_argument("--rpca", action = "store_true",
	help = "use rpca for integration")
parser$add_argument("--reference", action = "store", type = "character",
	nargs = "+", help = "samples to use as reference for rpca")

args <- parser$parse_args()

samples.list <- lapply(args$rds, function(x) {
	print(x)
	readRDS(x)})

if (args$rpca) {
	print("integrating with rpca")
	names(samples.list) <- args$sample
	refs <- sapply(args$reference, function(x) {which(names(samples.list)==x)})
	features <- SelectIntegrationFeatures(samples.list)
	samples.list <- lapply(samples.list, function(x) {RunPCA(x, features = features)})
	samples.list <- PrepSCTIntegration(object.list = samples.list, anchor.features = features)
	anchors <- FindIntegrationAnchors(samples.list, reference = refs, 
		reduction = "rpca", normalization.method="SCT", anchor.features = features, dims = 1:30)
	integrated <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "SCT")
	integrated <- RunPCA(integrated)
	integrated <- RunUMAP(integrated, dims = 1:30)
	integrated <- FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
	integrated <- FindClusters(integrated, verbose = F)

	# Visualization 
	p1 <- DimPlot(integrated, reduction = "umap", label = TRUE)
	p2 <- DimPlot(integrated, reduction = "umap", group.by = "label")
	p3 <- DimPlot(integrated, reduction = "umap", split.by = "label", label = FALSE)

	# write 
	png(args$clusters_umap, units = "in", width = 9, height = 6, res = 300)
	print(p1)
	dev.off()

	png(args$conditions_umap, units = "in", width = 9, height = 6, res = 300)
	print(p2)
	dev.off()

	png(args$split_umap, units = "in", width = 16, height = 8, res = 300)
	print(p3)
	dev.off()

	saveRDS(integrated, file = args$out_rds)

} else {
	print("integrating with SCT")
	features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
	samples.list <- PrepSCTIntegration(object.list = samples.list, anchor.features = features)
	anchors <- FindIntegrationAnchors(object.list = samples.list, normalization.method = "SCT", 
	anchor.features = features)
	integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

	integrated <- RunPCA(integrated, verbose = FALSE)
	integrated <- RunUMAP(integrated, dims = 1:30)
	integrated <- FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
	integrated <- FindClusters(integrated, verbose = FALSE)

	# Visualization 
	p1 <- DimPlot(integrated, reduction = "umap", label = TRUE)
	p2 <- DimPlot(integrated, reduction = "umap", group.by = "label")
	p3 <- DimPlot(integrated, reduction = "umap", split.by = "label", label = FALSE)

	# write 
	png(args$clusters_umap, units = "in", width = 9, height = 6, res = 300)
	print(p1)
	dev.off()

	png(args$conditions_umap, units = "in", width = 9, height = 6, res = 300)
	print(p2)
	dev.off()

	png(args$split_umap, units = "in", width = 16, height = 8, res = 300)
	print(p3)
	dev.off()

	saveRDS(integrated, file = args$out_rds)
}



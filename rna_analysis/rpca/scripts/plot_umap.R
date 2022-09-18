library(Seurat)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument('--rds', action = "store", type = "character", 
	help = "path to rds file for integrated object")
parser$add_argument("--clusters_umap", action = "store", type = "character",
	help = "where to save UMAP plot colored by cluster (png)")
parser$add_argument("--conditions_umap", action = "store", type = "character",
	help = "where to save UMAP plot colored by condition (png)")
parser$add_argument("--split_umap", action = "store", type = "character",
	help = "where to save side by side UMAP plot")

args <- parser$parse_args()

data <- readRDS(args$rds)

# Visualization 
p1 <- DimPlot(data, reduction = "umap", label = TRUE)
p2 <- DimPlot(data, reduction = "umap", group.by = "label")
p3 <- DimPlot(data, reduction = "umap", split.by = "label", label = FALSE)


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

####################################################################################
# this script is flexible for getting barcodes by rat, by celltype, or by both
####################################################################################


library(argparse)
library(Signac)
library(Seurat)
library(stringr)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rds", action = "store", type = "character",
	help = "path to labeled rds file for integrated object")
parser$add_argument("--celltype", action = "store", type = "character",
	default = NULL,
	help = "name of celltype")
parser$add_argument("--sample", action = "store", type = "character",
	default = NULL,
    help = "which rat to get barcodes for")
parser$add_argument("--out", action = "store", type = "character",
	help = "output directory")

print("parsing args")
args <- parser$parse_args()
print(args)

atac <- readRDS(args$rds)
DefaultAssay(atac) <- "peaks"
Idents(atac) <- atac$predicted.id

if (!(is.null(args$celltype)) & !(is.null(args$celltype))) {
	cells <- WhichCells(atac, expression = (sample == args$sample) & (predicted.id == args$celltype))
	cells <- str_remove(cells, paste0(args$sample, "_"))
	writeLines(cells, file.path(args$out, sprintf("%s-%s_barcodes.tsv", 
		args$sample, gsub("/", "-", args$celltype)))))
} else if (!(is.null(args$celltype))) {
	cells <- WhichCells(atac, expression = predicted.id == args$celltype)
	writeLines(cells, file.path(args$out, sprintf("%s_barcodes.tsv", gsub("/", "-", args$celltype)))))
} else if (!(is.null(args$sample))) {
	cells <- WhichCells(atac, expression = sample == args$sample)
	cells <- str_remove(cells, paste0(args$sample, "_"))
	writeLines(cells, file.path(args$out, sprintf("%s_barcodes.tsv", args$sample)))
}

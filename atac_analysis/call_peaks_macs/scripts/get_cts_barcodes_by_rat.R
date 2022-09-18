####################################################################################
# this script subsets by rat but iterates through all celltypes
####################################################################################
library(argparse)
library(Signac)
library(Seurat)
library(stringr)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rds", action = "store", type = "character",
	help = "path to labeled rds file for integrated object")
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

atac <- subset(atac, subset = sample == args$sample)

if (!dir.exists(args$out)) {
	dir.create(args$out, recursive = TRUE)
}
for (ct in Idents(atac)) {
	cat(sprintf("getting barcodes for %s\n", ct))
	cells <- WhichCells(atac, expression = predicted.id == ct)
	cells <- str_remove(cells, paste0(args$sample, "_"))
	writeLines(cells, file.path(args$out, sprintf("%s_barcodes.tsv", gsub("/", "-", ct))))
}

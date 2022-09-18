library(Seurat)
library(Signac)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rds", action = "store", type = "character",
	help = "path to labeled rds file for integrated object")
parser$add_argument("--prefilter", action = "store_true",
    help = "if true, prefilter genes for DGE analysis")
parser$add_argument("--covars", action = "store", nargs = "+",
    help = "covariates to include in DEG test")
parser$add_argument("--name", action = "store", 
    help = "name of analysis")
parser$add_argument("--celltype", action = "store", 
    help = "name of celltype")
parser$add_argument("--groupby", action = "store", default = NULL,
    help = "argument for group.by argument of Seurat FindMarkers")
parser$add_argument("--ident1", action = "store", default = "cocaine_high",
    help = "identity class to define markers for (corresponds to group.by argument)")
parser$add_argument("--ident2", action = "store", default = NULL,
    help = "(optional) second identity class for comparison; if NULL use all other identity classes")
parser$add_argument("--outfh", action = "store", type = "character",
    help = "where to save outputs")
parser$add_argument("--test", action = "store", type = "character", default = "wilcox",
	help = "test to use for FindMarkers() function")
parser$add_argument("--permute", action = "store_true", 
    help = "if true, run permuted test")

args <- parser$parse_args()
print(args)

# load RDS
data <- readRDS(args$rds)

DefaultAssay(data) <- "peaks"
Idents(data) <- data$predicted.id

# data$condition <- as.factor(data$condition)
data$library <- as.factor(data$library)
data$sample <- as.factor(data$sample)

celltype <- args$celltype
if (args$celltype == "Cck+-Vip+") {
    celltype <- "Cck+/Vip+"
}

if (!dir.exists(dirname(args$outfh))) {
    dir.create(dirname(args$outfh), recursive = T)
}

if (args$permute) {
    print("shuffling labels")
    data$shuf <- as.character(sample(data$condition, length(data$condition)))
}

cat(sprintf("computing DA peaks for %s\n", celltype))
if (args$prefilter) {
    if (args$permute) {
        print("permuting")
        dapeaks <- FindMarkers(data,
            test.use = args$test, ident.1 = args$ident1, ident.2 = args$ident2,
            latent.vars = args$covars, group.by = "shuf", subset.ident = celltype,
            min.pct = 0, logfc.threshold = 0)
        dapeaks$celltype <- celltype
        print(head(dapeaks))
        cat(sprintf("saving to %s\n", args$outfh))
        write.csv(dapeaks, file = args$outfh, quote = FALSE, row.names = TRUE, col.names = TRUE)
    } else {
        print("finding DA peaks")
        dapeaks <- FindMarkers(data,
            test.use = args$test, ident.1 = args$ident1, ident.2 = args$ident2,
            latent.vars = args$covars, group.by = args$groupby, subset.ident = celltype,
            min.pct = 0, logfc.threshold = 0)
        dapeaks$celltype <- celltype
        print(head(dapeaks))
        cat(sprintf("saving to %s\n", args$outfh))
        write.csv(dapeaks, file = args$outfh, quote = FALSE, row.names = TRUE, col.names = TRUE)
    }
} else {
    if (args$permute) {
        print('permuting')
        dapeaks <- FindMarkers(data,
            test.use = args$test, ident.1 = args$ident1, ident.2 = args$ident2,
            latent.vars = args$covars, group.by = "shuf", subset.ident = celltype,
            logfc.threshold = 0, min.pct = 0)
        dapeaks$celltype <- celltype
        cat(sprintf("saving to %s\n", args$outfh))
        write.csv(dapeaks, file = args$outfh, quote = FALSE, row.names = TRUE, col.names = TRUE)
    } else {
        print("finding DA peaks")
        dapeaks <- FindMarkers(data,
            test.use = args$test, ident.1 = args$ident1, ident.2 = args$ident2,
            latent.vars = args$covars, group.by = args$groupby, subset.ident = celltype,
            logfc.threshold = 0, min.pct = 0)
        dapeaks$celltype <- celltype
        cat(sprintf("saving to %s\n", args$outfh))
        write.csv(dapeaks, file = args$outfh, quote = FALSE, row.names = TRUE, col.names = TRUE)
    }
}

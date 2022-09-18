library(Seurat)
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
parser$add_argument("--groupby", action = "store", default = NULL,
    help = "argument for group.by argument of Seurat FindMarkers")
parser$add_argument("--ident1", action = "store", default = "cocaine_high",
    help = "identity class to define markers for (corresponds to group.by argument)")
parser$add_argument("--ident2", action = "store", default = NULL,
    help = "(optional) second identity class for comparison; if NULL use all other identity classes")
parser$add_argument("--outdir", action = "store", type = "character",
    help = "where to save outputs")
parser$add_argument("--test", action = "store", type = "character", default = "wilcox",
	help = "test to use for FindMarkers() function")

args <- parser$parse_args()
print(args)

data <- readRDS(args$rds)
# DefaultAssay(data) <- "SCT"
DefaultAssay(data) <- "RNA"
data <- NormalizeData(data)
data <- ScaleData(data)
# if (!(args$ident1 %in% data@meta.data[,args$groupby])) {
#     stop(sprintf("ident.1 %s not a value of groupby argument %s" % args$ident1, args$groupby))
# }
data$batch.code <- as.factor(sapply(data$batch, function(x) { return (letters[x]) }))
data$sample <- as.factor(data$sample)

cat(sprintf("outdir = %s\n", args$outdir))

if (!dir.exists(args$outdir)) {
    cat(sprintf("creating dir %s\n", args$outdir))
    dir.create(args$outdir, recursive = TRUE)
} else {
    cat(sprintf("!dir.exists(%s) returns FALSE\n", args$outdir))
}

for (type in names(table(Idents(data)))) {
    if (!file.exists(file.path(args$outdir, paste0(gsub(" ", "_", type), "_", args$name, ".csv")))) {
        cat(sprintf("computing DEG for %s\n", type))
        if (args$prefilter) {
            print("running analysis with prefiltering")
            degs <- FindMarkers(data,
                test.use = args$test, ident.1 = args$ident1, ident.2 = args$ident2,
                latent.vars = args$covars, group.by = args$groupby, subset.ident = type,
                logfc.threshold = 0, min.pct= 0.1)
            degs$celltype <- type
            outfh <- file.path(args$outdir, paste0(gsub("/", "-", type), "_", args$name, ".csv"))
            cat(sprintf("saving to %s\n", outfh))
            print(head(degs))
            write.csv(degs, file = outfh, quote = FALSE, row.names = TRUE, col.names = TRUE)
        } else {
            print("running analysis without prefiltering")
            degs <- FindMarkers(data,
                test.use = args$test, ident.1 = args$ident1, ident.2 = args$ident2,
                latent.vars = args$covars, group.by = args$groupby, subset.ident = type,
                logfc.threshold = 0, min.pct = 0, assay = "RNA")
            degs$celltype <- type
            outfh <- file.path(args$outdir, paste0(gsub("/", "-", type), "_", args$name, ".csv"))
            cat(sprintf("saving to %s\n", outfh))
            print(head(degs))
            write.csv(degs, file = outfh, quote = FALSE, row.names = TRUE, col.names = TRUE)
        }
    } else {
        cat(sprintf("%s already exists",
            file.path(args$outdir, paste0(gsub(" ", "_", type), "_", args$name, ".csv"))))
    }    
}

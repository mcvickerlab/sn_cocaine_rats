library(Seurat)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rds", action = "store", type = "character",
	help = "path to labeled rds file for integrated object")
parser$add_argument("--prefilter", action = "store_true",
    help = "if true, prefilter genes for DGE analysis")
parser$add_argument("--covars", action = "store", nargs = "+",
    help = "covariates to include in DEG test")
parser$add_argument("--outdir", action = "store", type = "character",
    help = "where to save outputs")
parser$add_argument("--test", action = "store", type = "character", 
	default = "wilcox",
	help = "test to use for FindMarkers() function")

args <- parser$parse_args()
print("args parsed")
print(args)

data <- readRDS(args$rds)
DefaultAssay(data) <- "RNA"
data <- NormalizeData(data)

for (type in names(table(Idents(data)))) {
    cat(sprintf("computing DEG for %s\n", type))
    
    if (args$prefilter) {
        # naive vs. cocaine treatment
        naive.vs.cocaine <- FindMarkers(data, 
            test.use = args$test, ident.1 = "naive", latent.vars = c("percent.mt", "sample", "label"),
            group.by = "label", subset.ident = type,
            logfc.threshold = 0)
        outfh <- file.path(args$outdir,
            # paste0(gsub(" ", "_", gsub(pattern = "[[:punct:]]", replacement = "", type)), 
            #     "_naive_vs_cocaine.csv"))
                paste0(gsub(" ", "_", type), "_naive_vs_cocaine.csv"))
        cat(sprintf("saving to %s\n", outfh))
        write.csv(naive.vs.cocaine, file = outfh, quote = F, row.names = TRUE, col.names = TRUE)

        # naive vs. cocaine low
        naive.vs.low <- FindMarkers(data,
            test.use = args$test, ident.1 = "naive", ident.2 = "cocaine_low", group.by = "label",
            latent.vars = c("percent.mt", "sample"), logfc.threshold = 0, subset.ident = type)
        outfh <- file.path(args$outdir,
            # paste0(gsub(" ", "_", gsub(pattern = "[[:punct:]]", replacement = "", type)), 
            paste0(gsub(" ", "_", type), "_naive_vs_low.csv"))
        cat(sprintf("saving to %s\n", outfh))
        write.csv(naive.vs.low, file = outfh, quote = F, row.names = TRUE, col.names = TRUE)

        # naive vs. cocaine high
        naive.vs.high <- FindMarkers(data,
            test.use = args$test, ident.1 = "naive", ident.2 = "cocaine_high", group.by = "label",
            latent.vars = c("percent.mt", "sample"),logfc.threshold = 0, subset.ident = type)
        outfh <- file.path(args$outdir,
            # paste0(gsub(" ", "_", gsub(pattern = "[[:punct:]]", replacement = "", type)), 
            paste0(gsub(" ", "_", type), "_naive_vs_high.csv"))
        cat(sprintf("saving to %s\n", outfh))
        write.csv(naive.vs.high, file = outfh, quote = F, row.names = TRUE, col.names = TRUE)

        # cocaine low vs. high
        low.vs.high <- FindMarkers(data, 
            test.use = args$test, ident.1 = "cocaine_low", ident.2 = "cocaine_high", group.by = "label", 
            latent.vars = c("percent.mt", "sample"),subset.ident = type, logfc.threshold = 0)
        outfh <- file.path(args$out,
            # paste0(gsub(" ", "_", gsub(pattern = "[[:punct:]]", replacement = "", type)), 
            paste0(gsub(" ", "_", type), "_cocaine_low_vs_high.csv"))
        cat(sprintf("saving to %s\n", outfh))
        write.csv(low.vs.high, file = outfh, quote = F, row.names = TRUE, col.names = TRUE)
    } else {
        # naive vs. cocaine treatment
        naive.vs.cocaine <- FindMarkers(data, 
            test.use = args$test, latent.vars = "percent.mt",
                         ident.1 = "naive", group.by = "label", subset.ident = type, 
                         latent.vars = c("percent.mt", "sample", "label"),
                         logfc.threshold = 0, min.pct = 0)
        outfh <- file.path(args$outdir,
            # paste0(gsub(" ", "_", gsub(pattern = "[[:punct:]]", replacement = "", type)), 
            paste0(gsub(" ", "_", type), "_naive_vs_cocaine.csv"))
        cat(sprintf("saving to %s\n", outfh))
        write.csv(naive.vs.cocaine, file = outfh, quote = F, row.names = TRUE, col.names = TRUE)

        # naive vs. cocaine low
        naive.vs.low <- FindMarkers(data,
            test.use = args$test, latent.vars = "percent.mt",
            ident.1 = "naive", ident.2 = "cocaine_low", group.by = "label",
            latent.vars = c("percent.mt", "sample"),
            subset.ident = type, logfc.threshold = 0, min.pct = 0)
        outfh <- file.path(args$outdir,
            # paste0(gsub(" ", "_", gsub(pattern = "[[:punct:]]", replacement = "", type)), 
            paste0(gsub(" ", "_", type), "_naive_vs_low.csv"))
        cat(sprintf("saving to %s\n", outfh))
        write.csv(naive.vs.low, file = outfh, quote = F, row.names = TRUE, col.names = TRUE)

        # naive vs. cocaine high
        naive.vs.high <- FindMarkers(data,
            test.use = args$test, latent.vars = "percent.mt",
            ident.1 = "naive", ident.2 = "cocaine_high", group.by = "label",
            latent.vars = c("percent.mt", "sample"),
            subset.ident = type, logfc.threshold = 0, min.pct = 0)
        outfh <- file.path(args$outdir,
            # paste0(gsub(" ", "_", gsub(pattern = "[[:punct:]]", replacement = "", type)), 
            paste0(gsub(" ", "_", type), "_naive_vs_high.csv"))
        cat(sprintf("saving to %s\n", outfh))
        write.csv(naive.vs.high, file = outfh, quote = F, row.names = TRUE, col.names = TRUE)

        # cocaine low vs. high
        low.vs.high <- FindMarkers(data, test.use = args$test, latent.vars = "percent.mt",
            ident.1 = "cocaine_low", ident.2 = "cocaine_high", group.by = "label", 
            latent.vars = c("percent.mt", "sample"),
            subset.ident = type, logfc.threshold = 0, min.pct = 0)
        outfh <- file.path(args$outdir,
            paste0(gsub(" ", "_", gsub(pattern = "[[:punct:]]", replacement = "", type)), 
                "_cocaine_low_vs_high.csv"))
        cat(sprintf("saving to %s\n", outfh))
        write.csv(low.vs.high, file = outfh, quote = F, row.names = TRUE, col.names = TRUE)
    }

}
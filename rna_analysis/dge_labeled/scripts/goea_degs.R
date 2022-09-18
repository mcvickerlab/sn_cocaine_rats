library(msigdbr)
library(fgsea)
library(ggplot2)
library(dplyr)
library(argparse)

parser <- ArgumentParser(help = "parse input args")

# list DEG files
low.vs.high.files <- list.files(path = "../rpca/out/grant_figure/dge",
                               pattern = "*cocaine_low_vs_high_dge.csv",
                               full.names = TRUE)

naive.vs.cocaine.files <- list.files(path = "../rpca/out/grant_figure/dge",
                                    pattern = "*naive_vs_cocaine_dge.csv",
                                    full.names = TRUE)

names(low.vs.high.files) <- sapply(low.vs.high.files, function(x) {
    fh <- strsplit(x, "/")[[1]][6]
    substr(fh, 1, (regexpr("_cocaine_low_vs_high_dge.csv", fh)-1))
})

names(naive.vs.cocaine.files) <- sapply(naive.vs.cocaine.files, function(x) {
    fh <- strsplit(x, "/")[[1]][6]
    substr(fh, 1, (regexpr("_naive_vs_cocaine_dge.csv", fh)-1))
})

# retrieve gene sets for Rattus norvegicus
# get GO:BP
go.bp <- msigdbr(species = "Rattus norvegicus", subcategory = "GO:BP")
go.bp_list <- split(x = go.bp$gene_symbol, f = go.bp$gs_name)
# get KEGG
kegg <- msigdbr(species = "Rattus norvegicus", subcategory = "CP:KEGG")
kegg_list <- split(x = kegg$gene_symbol, f = kegg$gs_name)

cat("analyzing low vs. high\n")
for (i in 1:length(low.vs.high.files)) {
	file <- low.vs.high.files[i]
	print(file)
	celltype <- names(low.vs.high.files)[i]
	cat(sprintf("analyzing %s\n", celltype))
	df <- read.csv(file, row.names = 1)
	ranks <- df$avg_logFC
	names(ranks) <- rownames(df)
	fgsea.go <- fgsea(pathways = go.bp_list, 
		stats = ranks,
		minSize=10,
		maxSize=500,
		nperm=1000000)
	fgsea.kegg <- fgsea(pathways = kegg_list,
		stats = ranks,
		minSize = 10,
		maxSize = 500,
		nperm = 1000000)
	fgsea.go <- fgsea.go[order(pval),]
	# print(head(fgsea.go))
	# print(nrow(fgsea.go))
	fgsea.kegg <- fgsea.kegg[order(pval),]
	if (nrow(fgsea.go)>0) {
		print("writing go")
		# print(head(fgsea.go))
		fgsea.go$leadingEdge <- vapply(fgsea.go$leadingEdge, paste, collapse = ", ", character(1L))
		print(head(fgsea.go))
		# print(class(fgsea.go$leadingEdge))
		fh <- sprintf("/iblm/netapp/home/jezhou/rat_snrnaseq_pipeline/rpca/out/grant_figure/goea/low_vs_high/go/%s_low_vs_high_go.txt", celltype)
		print(fh)
		write.table(fgsea.go, 
			file = file.path(fh),
			quote = FALSE, sep = "\t", row.names = FALSE)
		# print("done writing go")
	}
	if (nrow(fgsea.kegg)>0) {
		print("writing kegg")
		fgsea.kegg$leadingEdge <- vapply(fgsea.kegg$leadingEdge, paste, collapse = ", ", character(1L))
		fh <- sprintf("/iblm/netapp/home/jezhou/rat_snrnaseq_pipeline/rpca/out/grant_figure/goea/low_vs_high/kegg/%s_low_vs_high_kegg.txt", celltype)
		write.table(fgsea.kegg, 
			file = file.path(fh),
			quote = FALSE, sep = "\t", row.names = FALSE)
	}


}

cat("analyzing naive vs. cocaine\n")
for (i in 1:length(naive.vs.cocaine.files)) {
	file <- naive.vs.cocaine.files[i]
	celltype <- names(naive.vs.cocaine.files)[i]
	cat(sprintf("analyzing %s\n", celltype))
	df <- read.csv(file, row.names = 1)
	ranks <- df$avg_logFC
	names(ranks) <- rownames(df)
	fgsea.go <- fgsea(pathways = go.bp_list, 
		stats = ranks,
		minSize=10,
		maxSize=500,
		nperm=1000000)
	fgsea.kegg <- fgsea(pathways = kegg_list,
		stats = ranks,
		minSize = 10,
		maxSize = 500,
		nperm = 1000000)
	fgsea.go <- fgsea.go[order(pval),]
	fgsea.kegg <- fgsea.kegg[order(pval),]
	if (nrow(fgsea.go)>0) {
		
		fgsea.go$leadingEdge <- vapply(fgsea.go$leadingEdge, paste, collapse = ", ", character(1L))
		sprintf("/iblm/netapp/home/jezhou/rat_snrnaseq_pipeline/rpca/out/grant_figure/goea/naive_vs_cocaine/go/%s_naive_vs_cocaine_go.txt", celltype)
		write.table(fgsea.go, 
			file = file.path(fh),
			quote = FALSE, sep = "\t", row.names = FALSE)
	}
	if (nrow(fgsea.kegg)>0) {
		fgsea.kegg$leadingEdge <- vapply(fgsea.kegg$leadingEdge, paste, collapse = ", ", character(1L))
		fh <- sprintf("/iblm/netapp/home/jezhou/rat_snrnaseq_pipeline/rpca/out/grant_figure/goea/naive_vs_cocaine/kegg/%s_naive_vs_cocaine_kegg.txt", celltype)
		write.table(fgsea.kegg, 
			file = file.path(fh),
			quote = FALSE, sep = "\t", row.names = FALSE)
	}


}

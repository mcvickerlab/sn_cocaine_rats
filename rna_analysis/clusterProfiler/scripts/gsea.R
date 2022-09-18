suppressMessages(library(dplyr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Rn.eg.db))
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--markers", action = "store", type = "character",
	help = "path to file with output of FindAllMarkers")
parser$add_argument("--zscore", action = "store_true",
	help = "if set, compute GSEA using z-score of avg_logFC")
parser$add_argument("--celltype", action = "store",
	help = "which celltype to run GSEA on")
parser$add_argument("--kegg", action = "store_true",
	help = "if set, run GSEA with KEGG")
parser$add_argument("--go", action = "store_true",
	help = "if set, run GSEA with GO:BP")
parser$add_argument("--cutoff", action = "store", type = "double",
	default = 0.1,
	help = "p-value cutoff for reporting results")
parser$add_argument("--padjust", action = "store",
	default = "fdr",
	help = "padjust method")
parser$add_argument("--outdir", action = "store", type = "character",
    help = "where to save outputs")

args <- parser$parse_args()
print(args)

# add genes as column and q_val
markers <- read.csv(args$markers, row.names = 1)
markers <- tibble::rownames_to_column(markers, "gene")
markers$q_val <- p.adjust(markers$p_val, method = "fdr")

if (args$zscore) {
	# add z-scores for avg_logFC
	markers <- markers %>% mutate(zscore = (avg_log2FC - mean(avg_log2FC))/sd(avg_log2FC))
}

print("markers")
print(head(markers))

gene.df <- bitr(markers$gene, fromType = "SYMBOL",
               # toType = c("ENSEMBL", "ENTREZID"),
               toType = "ENTREZID",
               OrgDb = org.Rn.eg.db)

gene.list <- NULL

if (args$zscore) {
	gene.list <- markers$zscore
} else {
	gene.list <- markers$avg_log2FC
}
names(gene.list) <- markers$gene
gene.list <- sort(gene.list, decreasing = TRUE)


if (args$go) {
	print("running GSEA for GO:BP")
	gobp <- gseGO(geneList = gene.list,
	     ont = "BP",
	     OrgDb = org.Rn.eg.db,
	     keyType = "SYMBOL",
	     pvalueCutoff = args$cutoff,
	     pAdjustMethod = args$padjust)
	outfh <- NULL
	if (args$zscore) {
		outfh <- file.path(args$outdir, sprintf("%s_gobp_zscore.rds", args$celltype))
	} else {
		outfh <- file.path(args$outdir, sprintf("%s_gobp.rds", args$celltype))
	}
	cat(sprintf("saving gobp results to %s\n", outfh))
	saveRDS(gobp, outfh)
	table <- NULL
	if (args$zscore) {
		table <- file.path(args$outdir, sprintf("%s_gobp_zscore.txt", args$celltype))
	} else {
		table <- file.path(args$outdir, sprintf("%s_gobp.txt", args$celltype))
	}
	cat(sprintf("saving gobp table to %s\n", table))
	write.table(data.frame(gobp), file = table,
		row.names = FALSE, col.names = TRUE, sep = '\t')

}

if (args$kegg) {
	print("running GSEA for KEGG")
	# get ENTREZID
	entrez.ids<-bitr(names(gene.list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Rn.eg.db)
	entrez.ids <- entrez.ids[!duplicated(entrez.ids[c("SYMBOL")]),]
	markers.with.ids <- markers %>% left_join(entrez.ids, by = c("gene" = "SYMBOL"))
	markers.with.ids <- na.omit(markers.with.ids)
	kegg.ranks <- c()
	if (args$zscore) {
		kegg.ranks <- markers.with.ids$zscore
	} else {
		kegg.ranks <- markers.with.ids$avg_log2FC
	}
	names(kegg.ranks) <- markers.with.ids$ENTREZID
	kegg.ranks <- sort(kegg.ranks, decreasing = TRUE)
	kegg <- gseKEGG(geneList     = kegg.ranks,
               organism     = "rno",
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = args$cutoff,
               pAdjustMethod = args$padjust,
               keyType       = "ncbi-geneid")
	outfh <- NULL
	if (args$zscore) {
		outfh <- file.path(args$outdir, sprintf("%s_kegg_zscore.rds", args$celltype))
	} else {
		outfh <- file.path(args$outdir, sprintf("%s_kegg.rds", args$celltype))
	}
	cat(sprintf("saving kegg results to %s\n", outfh))
	saveRDS(kegg, outfh)
	table <- NULL
	if (args$zscore) {
		table <- file.path(args$outdir, sprintf("%s_kegg_zscore.txt", args$celltype))
	} else {
		table <- file.path(args$outdir, sprintf("%s_kegg.txt", args$celltype))
	}
	cat(sprintf("saving kegg table to %s\n", table))
	write.table(data.frame(kegg), file = table,
		row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
}



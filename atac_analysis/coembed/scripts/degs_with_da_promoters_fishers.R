suppressMessages(library(dplyr))
suppressMessages(require(TxDb.Rnorvegicus.UCSC.rn6.refGene))
require(org.Rn.eg.db)
library(annotate)
library(GenomicRanges)
library(argparse)
library(broom)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--celltype", action = "store",
	help = "name of celltype")
parser$add_argument("--degs", action = "store",
	help = "path to file with DEGs")
parser$add_argument("--dapeaks", action = "store",
	help = "path to file with DA peaks")
parser$add_argument("--upstream", action = "store",
	default = 1500,
	help = "value to pass to upstream argument of promoters function")
parser$add_argument("--downstream", action = "store",
	default = 500,
	help = "value to pass to downstream argument of promoters function")
parser$add_argument("--outfh", action = "store",
	help = "where to save Fishers test results (csv)")

args <- parser$parse_args()

###  load DEGs
degs <- read.csv(args$degs, row.names = 1)

print(head(degs))

if (!("q_val" %in% colnames(degs))) {
	degs$q_val <- p.adjust(degs$p_val, method = "fdr")
}

### load DA peaks
dapeaks <- read.csv(args$dapeaks, row.names = 1)

if (!("q_val" %in% colnames(dapeaks))) {
	dapeaks$q_val <- p.adjust(dapeaks$p_val, method = "fdr")
}

### promoter info
txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
promoters.gr <- promoters(genes(txdb), upstream = 1500, downstream = 500)
promoters.gr$symbol <- lookUp(promoters.gr$gene_id, 'org.Rn.eg.db', 'SYMBOL') 


### get list of significant DEGs
significant.degs <- row.names(degs %>% filter(q_val<0.1))

### get promoter regions of significant DEGs
significant.degs.pr <- promoters.gr[promoters.gr$symbol %in% significant.degs]

### get GRanges for significant DA peaks
significant.peaks <- rownames(dapeaks %>% filter(q_val<0.1))
significant.peaks <- as.data.frame(t(as.data.frame(strsplit(significant.peaks, '-'))))
colnames(significant.peaks) <- c("seqnames", "start", "end")
significant.peaks.gr <- GRanges(significant.peaks)
significant.peaks.gr

### `findOverlaps` between promoters of significant DEGs and DA peaks to identify DEGs with DA promoters
deg.dap.ols <- findOverlaps(significant.degs.pr, significant.peaks.gr)

### get values for DEGs row of confusion matrix
n.deg.dap.ols <- length(unique(queryHits(deg.dap.ols)))
# n.deg.nondap.ols <- length(significant.degs)-n.deg.dap.ols
n.deg.nondap.ols <- queryLength(deg.dap.ols) - n.deg.dap.ols

### find non-DEGs and repeat process for non-DEGs
non.degs <- rownames(filter(degs, q_val>0.1))

### get promoter regions of non-DEGs
non.degs.pr <- promoters.gr[promoters.gr$symbol %in% non.degs]

### `findOverlaps` between promoters of non-significant DEGs and DA peaks to identify non-DEGs with DA promoters
nondeg.dap.ols <- findOverlaps(non.degs.pr, significant.peaks.gr)

### get values for non-DEGs row of confusion matrix
n.nondeg.dap.ols <- length(unique(queryHits(nondeg.dap.ols)))
# n.deg.nondap.ols <- length(significant.degs)-n.deg.dap.ols
n.nondeg.nondap.ols <- queryLength(nondeg.dap.ols) - n.nondeg.dap.ols

### get matrix for chisq test
mtx <- data.frame("DA promoter" = c(n.deg.dap.ols,n.deg.nondap.ols), 
	"non-DA promoter" = c(n.nondeg.dap.ols,n.nondeg.nondap.ols),
	row.names = c("DEG", "non-DEG"))

print(args$celltype)
print(mtx)

res <- fisher.test(mtx)
res.tidy <- tidy(res)
res.tidy$celltype <- args$celltype

write.csv(res.tidy, args$outfh, row.names = FALSE, quote = FALSE)

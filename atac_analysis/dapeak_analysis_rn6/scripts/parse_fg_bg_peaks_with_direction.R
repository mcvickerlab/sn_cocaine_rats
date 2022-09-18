suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))
suppressMessages(library(rtracklayer))
suppressMessages(library(Signac))
library(GenomeInfoDb)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--dapeaks", action = "store", type = "character",
	help = "path to csv file for dapeaks")
parser$add_argument("--celltype", action = "store", 
    help = "celltype")
parser$add_argument("--outdir", action = "store", 
    help = "where to save outputs")

args <- parser$parse_args()
print(args)

df <- read.csv(args$dapeaks, row.names = 1)
df <- tibble::rownames_to_column(df, "peaks")
df$q_val <- p.adjust(df$p_val, method = "fdr")
significant <- filter(df, q_val<0.1)
pos <- filter(significant, avg_log2FC > 0)
neg <- filter(significant, avg_log2FC < 0)
pos.gr <- StringToGRanges(pos$peaks)
neg.gr <- StringToGRanges(neg$peaks)
seqlevelsStyle(pos.gr) <- "UCSC"
seqlevelsStyle(neg.gr) <- "UCSC"
# save pos
export.bed(pos.gr, con = file.path(args$outdir, sprintf("%s_fg_pos.bed", args$celltype)))
# save neg
export.bed(neg.gr, con = file.path(args$outdir, sprintf("%s_fg_neg.bed", args$celltype)))
# save background
bg <- (filter(df, q_val>0.5) %>% select(peaks))$peaks
bg.gr <- StringToGRanges(bg)
seqlevelsStyle(bg.gr) <- "UCSC"
export.bed(bg.gr, con = file.path(args$outdir, sprintf("%s_bg.bed", args$celltype)))

library(argparse)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Rnorvegicus.UCSC.rn6)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rds", action = "store", type = "character",
	help = "path to labeled rds file for integrated object")
parser$add_argument("--celltype", action = "store", type = "character",
	help = "name of celltype")
parser$add_argument("--outfh", action = "store", type = "character",
    help = "where to save outputs")
parser$add_argument("--bg", action = "store", type = "character",
    help = "path to where to save matched background sequences")
parser$add_argument("--dapeaks", action = "store", type = "character",
	help = "file with differential accessibility analysis output")
parser$add_argument("--alpha", action = "store", default = 0.1,
	help = "significance cutoff for selecting FG peaks")
parser$add_argument("--npeaks", action = "store", default = 10000,
    help = "number of peaks for background matching")

print("parsing args")
args <- parser$parse_args()
print(args)

# load RDS
print("loading ATAC")
atac <- readRDS(args$rds)
DefaultAssay(atac) <- "peaks"
Idents(atac) <- atac$predicted.id
table(Idents(atac))

# add RegionStats
print("adding RegionStats")
atac <- RegionStats(
    object = atac,
    genome = BSgenome.Rnorvegicus.UCSC.rn6,
    assay = "peaks"
)

# add motifs 
print("getting MatrixSet")
pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

print("creating motif matrix")
motif.matrix <- CreateMotifMatrix(
    features = granges(atac),
    pwm = pfm,
    genome = BSgenome.Rnorvegicus.UCSC.rn6
)

print("creating motif object")
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

print("adding motifs to object")
Motifs(atac) <- motif

# load DA peaks
print("loading da peaks")
dapeaks <- read.csv(args$dapeaks, row.names = 1)
if (!("q_val" %in% colnames(dapeaks))) {
	dapeaks$q_val <- p.adjust(dapeaks$p_val, "fdr")
}

# get top differentially accessible peaks (FG peaks)
print("getting top DA peaks")
top.da.peak <- rownames(dapeaks[dapeaks$q_val <= args$alpha, ])

# get BG peaks
print("getting matched BG peaks")
celltype <- args$celltype
if (args$celltype == "Cck+-Vip+") {
    celltype <- "Cck+/Vip+"
}
open.peaks <- AccessiblePeaks(atac, idents = celltype)
meta.feature <- GetAssayData(atac, assay = "peaks", slot = "meta.features")

peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = args$npeaks,
  features.match = c("GC.percent", "percentile")
)

cat(sprintf("%d peaks matched for %s\n", length(peaks.matched), celltype))
writeLines(peaks.matched, args$bg)

# find enriched motifs
print("finding enriched motifs")
enriched.motifs <- FindMotifs(
  object = atac,
  features = top.da.peak,
    background = peaks.matched
)

write.csv(enriched.motifs, args$outfh, quote = F, row.names = F)

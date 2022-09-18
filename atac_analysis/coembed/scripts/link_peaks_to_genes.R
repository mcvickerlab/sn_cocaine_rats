library(Signac)
library(Seurat)
library(argparse)
library(BSgenome.Rnorvegicus.UCSC.rn6)


parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rds", action = "store",
	help = "path to rds object to compute gene activity for")
parser$add_argument("--out", action = "store",
	help = "where to save rds object with gene activities")

args <- parser$parse_args()

dat <- readRDS(args$rds)
Idents(dat) <- dat$predicted.id 

dat <- RegionStats(
    object = dat,
    genome = BSgenome.Rnorvegicus.UCSC.rn6,
    assay = "peaks"
)
# linkage
links <- LinkPeaks(dat,
	peak.assay = "peaks",
	expression.assay = "RNA",
	expression.slot = "data",
	gene.coords = NULL,
	distance = 5e+05,
	min.distance = NULL,
	min.cells = 10,
	# method = "pearson",
	genes.use = NULL,
	n_sample = 200,
	pvalue_cutoff = 0.05,
	score_cutoff = 0.05,
	gene.id = FALSE,
	verbose = TRUE
)

# save
saveRDS(links, args$out)


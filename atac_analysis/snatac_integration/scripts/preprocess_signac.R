suppressMessages(library(Signac))
suppressMessages(library(Seurat))
library(argparse)
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(AnnotationHub))
suppressMessages(library(BSgenome.Rnorvegicus.UCSC.rn6))
library(rtracklayer)
library(ggplot2)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--mtx", action = "store", type = "character",
	help = "path to peak barcode matrix (h5)")
parser$add_argument("--metadata", action = "store", type = "character",
	help = "path to singlecell.csv output file (cellranger-atac)")
parser$add_argument("--fragments", action = "store", type = "character",
	help = "path to fragments.tsv.gz output file (cellranger-atac)")
parser$add_argument("--mincells", action = "store", default = 10,
	help = "min.cells argument for CreateChromatinAssay")
parser$add_argument("--minfeatures", action = "store", default = 200,
	help = "min.features argument for CreateChromatinAssay")
parser$add_argument("--vlnplt", action = "store", type = "character", 
	default = NULL,
	help = "where to save QC plots")
parser$add_argument("--rds", action = "store", type = "character",
	help = "where to save rds file for processed Signac object")
parser$add_argument("--umap", action = "store", type = "character",
	default = NULL,
	help = "where to save UMAP plot")
parser$add_argument("--sample", action = "store",
	help = "sample name")
parser$add_argument("--condition", action = "store",
	help = "drug condition, e.g. cocaine.high, cocaine.low, naive")
parser$add_argument("--library", action = "store", type = "character",
	help = "library construction date")
parser$add_argument("--ann", action = "store", type = "character",
	help = "path to gtf file with genome annotation")

args <- parser$parse_args()

# create Signac object
counts <- Read10X_h5(filename = args$mtx, use.names = F)

# filter chromosomes
seqs.keep <- c(paste0('chr', 1:20), "chrX", 'chrY')
counts.seqs <- lapply(rownames(counts), function(x) {
    strsplit(x, "-")[[1]][1]
})

# get 1-based coordinates
peaks.1based <- GRangesToString(grange = StringToGRanges(regions = rownames(counts)), 
	sep = c("-","-"), starts.in.df.are.0based = TRUE)
rownames(counts) <- peaks.1based

metadata <- read.csv(
  file = args$metadata,
  header = TRUE,
  row.names = 1
)

# get SeqInfo() object for rn6
rn6.seqinfo <- Seqinfo(seqnames(BSgenome.Rnorvegicus.UCSC.rn6), 
	seqlengths(BSgenome.Rnorvegicus.UCSC.rn6), 
	isCircular(BSgenome.Rnorvegicus.UCSC.rn6), 
	genome(BSgenome.Rnorvegicus.UCSC.rn6))

# get annotation
# ann <- rtracklayer::import("/iblm/netapp/data1/jezhou/Rattus_norvegicus.Rnor_6.0.98.filtered.gtf")
ann <- rtracklayer::import(args$ann)

chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = rn6.seqinfo,
    fragments = args$fragments,
    min.cells = args$mincells,
    min.features = args$minfeatures,
    annotation = ann
)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(chrom_assay),
  pwm = pfm,
  genome = BSgenome.Rnorvegicus.UCSC.rn6
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# add motif information
Motifs(chrom_assay) <- motif


# create Seurat object from ChromatinAssay
dat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

dat$sample <- args$sample
dat$condition <- args$condition
dat$library <- args$library

### QC
# get RegionStats
dat <- RegionStats(
  object = dat,
  genome = BSgenome.Rnorvegicus.UCSC.rn6,
    assay = "peaks"
)

# compute nucleosome signal score per cell
dat <- NucleosomeSignal(object = dat)

# compute TSS enrichment score per cell
dat <- TSSEnrichment(object = dat, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
dat$pct_reads_in_peaks <- dat$peak_region_fragments / dat$passed_filters * 100
# dat$blacklist_ratio <- dat$blacklist_region_fragments / dat$peak_region_fragments

# plot QC metrics
if (!is.null(args$vlnplt)) {
	vln <- VlnPlot(
		object = dat,
		features = c('pct_reads_in_peaks', 'peak_region_fragments',
	               'TSS.enrichment', 'nucleosome_signal'),
		pt.size = 0.1,
		ncol = 4
	)
	png(args$vlnplt, units = "in", width = 12, height = 6, res = 300)
	print(vln)
	dev.off()
}

# filter
print('filtering')
dat <- subset(dat,
	subset = peak_region_fragments>max(0,(mean(dat[[]]$peak_region_fragments) - 2*sd(dat[[]]$peak_region_fragments))) &
	peak_region_fragments<(mean(dat[[]]$peak_region_fragments) + 2*sd(dat[[]]$peak_region_fragments)) &
	pct_reads_in_peaks>max(0, (mean(dat[[]]$pct_reads_in_peaks) - 2*sd(dat[[]]$pct_reads_in_peaks))) &
	nucleosome_signal<(mean(dat[[]]$nucleosome_signal) + 2*sd(dat[[]]$nucleosome_signal)) &
	TSS.enrichment>(mean(dat[[]]$TSS.enrichment) - 2*sd(dat[[]]$TSS.enrichment))
	)

print(dat)

# normalization and linear dimensional reduction
print("TFIDF")
dat <- RunTFIDF(dat)
print("FindTopFeatures")
dat <- FindTopFeatures(dat, min.cutoff = 'q0')
print("RunSVD")
dat <- RunSVD(dat)
print("RunUMAP")
dat <- RunUMAP(object = dat, reduction = 'lsi', dims = 2:30)

# cluster plot UMAP
dat <- FindNeighbors(object = dat, reduction = 'lsi', dims = 2:30)
dat <- FindClusters(object = dat, verbose = FALSE, algorithm = 3)
umap <- DimPlot(object = dat, label = TRUE, repel = TRUE, label.size = 8) + 
NoAxes() + theme(legend.text=element_text(size=15))

if (!is.null(args$umap)) {
	png(args$umap, units = "in", width = 8, height = 6, res = 300)
	print(umap)
	dev.off()
}

saveRDS(dat, file = args$rds)
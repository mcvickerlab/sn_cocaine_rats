suppressMessages(library(Signac))
# suppressMessages(library(Seurat))
library(GenomicRanges)
library(Matrix)
suppressMessages(library(HDF5Array))
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--fragments", action = "store", type = "character",
	help = "path to fragments file for sample (cellranger output)")
parser$add_argument("--peaks", action = "store", type = "character",
	help = "path to peaks file for sample")
parser$add_argument("--barcodes", action = "store", type = "character", 
	default = NULL,
	help = "path to file containing barcodes to include in output mtx")
parser$add_argument("--out", action = "store", type = "character",
	help = "where to write output h5 matrix (must have .h5 extension)")
parser$add_argument("--ucsc", action = "store_true",
	help = "if set, convert chromosome names to UCSC format")

args <- parser$parse_args()

# load barcodes (if provided) and fragments
bcs <- NULL
if (!is.null(args$barcodes)) {
	print("barcodes provided")
	bcs <- readLines(args$barcodes)
}
frags <- CreateFragmentObject(args$fragments, cells = bcs)

# load MACS2 peaks
print("loading MACS2 peaks")
peaks.df <- read.table(args$peaks)

if (args$ucsc) {
	chroms <- c(as.character(1:20), "X", "Y", "MT")
	peaks.df$V1 <- sapply(peaks.df$V1, function(x) {
		if (x %in% chroms) {
			return(paste0("chr",x))
			} else {
				return(x)
			}
		})
}

peaks.gr <- GRanges(seqnames = peaks.df$V1,
	ranges = IRanges(start = peaks.df$V2, end = peaks.df$V3, names = peaks.df$V4),
	score = peaks.df$V5)

# generate new feature-bc mtx from peaks 
print("generating new feature barcode mtx")
mtx <- FeatureMatrix(fragments = frags,
                    features = peaks.gr,
                    cells = bcs)
print(dim(mtx))

# write mtx to h5 file
cat(sprintf("writing new matrix to h5 file at %s\n", args$out))
writeTENxMatrix(mtx, filepath = args$out, verbose = TRUE)
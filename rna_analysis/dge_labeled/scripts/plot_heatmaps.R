suppressMessages(library(Seurat))
library(ggplot2)
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
library(stringr)
library(colorspace)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--rds", action = "store", type = "character",
	help = "path to labeled rds file for integrated object")
parser$add_argument("--degs", action = "store", nargs = "+",
	help = "files with dge analysis results")

args <- parser$parse_args()


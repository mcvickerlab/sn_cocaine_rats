suppressMessages(library(dplyr))
library(ggplot2)
library(argparse)
library(patchwork)
library(argparse)

parser <- ArgumentParser(description = "process input arguments")
parser$add_argument("--obs", action = "store", type = "character",
    help = "path to file with observed results")
parser$add_argument("--perm", action = "store", type = "character",
    help = "path to file with permuted results")
parser$add_argument("--celltype", action = "store", 
    help = "name of celltype")
parser$add_argument("--outfh", action = "store", type = "character",
    help = "where to save outputs")

args <- parser$parse_args()
print(args)

print("loading data")
obs <- read.csv(args$obs, row.names = 1)
perm <- read.csv(args$perm, row.names = 1)


print("obs - hist")
obs.hist <- ggplot(obs, aes(x=p_val)) + 
geom_histogram(color="black", fill="white") + theme_bw(24) +   
theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
  ) + ggtitle("Observed")

print("perm - hist")
perm.hist <- ggplot(perm, aes(x=p_val)) + 
geom_histogram(color="black", fill="white") + theme_bw(24) +   
theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
  ) + ggtitle("Permuted")

gg_qqplot <- function(ps.obs, ps.perm, ci = 0.95) {
    n  <-  max(length(ps.obs), length(ps.perm))
    
    observed <- -log10(sort(ps.obs))
    permuted <- -log10(sort(ps.perm))
    expected <- -log10(ppoints(n))
    
    observed.points = c(observed, rep(NA, n - length(observed)))
    permuted.points = c(permuted, rep(NA, n - length(permuted)))
    
    df <- data.frame(observed = observed.points,
                permuted = permuted.points,
                expected = expected,
                clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
                cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)))
    
    plotdf <- tidyr::pivot_longer(df, cols = c("observed", "permuted"), 
      names_to = "data", values_to = "-log10(p)")

    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("-log"[10], plain(P)))

    ggplot(plotdf) +
        geom_ribbon(
          mapping = aes(x = expected, ymin = clower, ymax = cupper),
          alpha = 0.1
        ) +
        geom_point(aes(x = expected, y = `-log10(p)`, color = data), size = 3) +
        geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
        geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
        geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
        xlab(log10Pe) + 
        ylab(log10Po) +   theme_bw(base_size = 24) +
      theme(
        axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank(),
          legend.title = element_blank()
      )

}

print("plotting qqplot")
p.gg <- gg_qqplot(obs$p_val, perm$p_val)

p.out <- p.gg | (obs.hist/perm.hist) + plot_annotation(title = args$celltype)

cat(sprintf("saving plot to %s\n", args$outfh))
png(args$outfh, res = 300, units = "in", width = 20, height = 10)
print(p.out)
dev.off()

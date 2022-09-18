library(Seurat)

filename <- '../out/grant_figure/rpca_integrated.rds'
data <- readRDS(filename)

DefaultAssay(data) <- "RNA"

cat("computing cluster 2\n")
deg.2 <- FindMarkers(data, ident.1 = 2, test.use = "negbinom",
                    latent.vars = c("label","percent.mt"))

cat("computing lcuster 35\n")
deg.35 <- FindMarkers(data, ident.1 = 35, test.use = "negbinom",
                     latent.vars = c("label","percent.mt"))

write.table(deg.2, file = "../out/grant_figure/cluster2_degs.txt", row.names = T, quote = F)

write.table(deg.35, file = "../out/grant_figure/cluster35_degs.txt", row.names = T, quote = F)

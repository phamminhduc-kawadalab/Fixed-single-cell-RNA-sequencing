library(Seurat)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tidyr)
library(piano)
library(msigdbr)
library(monocle)
library(SingleCellExperiment)

install.packages('piano')

req_subset <- readRDS("/Volumes/LaCie/For\ articles/NKT\ cell\ /NKT.rds")

req_subset <- readRDS("/Volumes/LaCie/20240530_HSC_subset/20240530_HSC_subset.rds")

req_subset <- readRDS("/Volumes/LaCie/20240530_Human\ MALFD/20240530_HSC_Human_MALFD.rds")

req_subset <- readRDS("/Volumes/LaCie/20251005_Data\ analysis\ /20250105_Data\ analysis\ /20250105_NKT.rds")

req_subset <- readRDS("/Users/minhduc/Desktop/20241027_MANUSCRIPT_Update/20241101_Core\ data/Supplement\ Table\ /Script\ for\ github\ publishing/20240102_Macrophages/Macrophages.rds")

DimPlot(req_subset)
req_subset@meta.data$group


SO.traj1 <- req_subset[,req_subset$group %in% c("TAA", "TAA_RE")]

req_subset@meta.data

pd <- new('AnnotatedDataFrame', data = SO.traj1@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(SO.traj1), row.names = row.names(SO.traj1)))

cds <- newCellDataSet(as(SO.traj1@assays$RNA@data, "matrix"), phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds$clusters <- SO.traj1$group

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

set.seed(123)
cds <- detectGenes(cds, min_expr = 0.1)

SO.traj1$orig.ident <- factor(SO.traj1$group, levels = c("TAA", "TAA_RE"))

SO.traj1 <- SetIdent(SO.traj1, value = "group")

markers <- FindAllMarkers(object = SO.traj1, min.pct = 0.25, thresh.use = 0.25)
markers <- subset(markers, p_val_adj < 0.05)

order.genes <- unique(as.character(markers$gene))

cds <- setOrderingFilter(cds, order.genes)
cds <- reduceDimension(cds = cds, max_components = 3,method = 'DDRTree')
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = T, cell_size = 0.5, label=T) + scale_x_reverse() + theme(axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), legend.text = element_text(size = 14), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20)) +  guides(color = guide_legend(override.aes = list(size=4)) ) 

plot_cell_trajectory(cds, color_by = "group", show_branch_points = T, cell_size = 0.5, label=T) + scale_x_reverse() + theme(axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), legend.text = element_text(size = 14), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20)) +  guides(color = guide_legend(override.aes = list(size=4)) ) 

plot_genes_branched_heatmap(cds, color_by = "group", branch_points = 1, branch_states = 1, cell_size = 0.5, label=T)


order.genes <- order.genes[!grepl("ENSMPUG", order.genes)]


plot_pseudotime_heatmap(cds[order.genes,], num_clusters = 4, cores = 1, show_rownames = T, return_heatmap = T)
dev.off()

diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res <- differentialGeneTest(cds[order.genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")


sig_gene_names <- row.names(subset(diff_test_res, qval < 1.596244e-10))

plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 2,
                        cores = 1,
                        show_rownames = T)

dev.off()

gene_short_name %in% c("Fasl")

plot_genes_in_pseudotime(cds, color_by = "group")

my_genes <- row.names(subset(fData(cds),
                             gene_short_name %in% c("Fasl")))
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Cell.types") +  guides(color = guide_legend(override.aes = list(size=4))) 

##

SO.traj1$group <- factor(SO.traj1$group, levels = c("CONT", "TAA", "TAARE"))

SO.traj1 <- SetIdent(SO.traj1, value = "group")

markers <- FindAllMarkers(object = SO.traj1, min.pct = 0.25, thresh.use = 0.25)
markers <- subset(markers, p_val_adj < 0.05)

order.genes <- unique(as.character(markers$gene))

cds <- setOrderingFilter(cds, order.genes)
cds <- reduceDimension(cds = cds, max_components = 3,method = 'DDRTree')
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = T, cell_size = 0.5, label=T) + scale_x_reverse() + theme(axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), legend.text = element_text(size = 14), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20)) +  guides(color = guide_legend(override.aes = list(size=4)) ) 

plot_cell_trajectory(cds, color_by = "group", show_branch_points = T, cell_size = 0.5, label=T) + scale_x_reverse() + theme(axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), legend.text = element_text(size = 14), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20)) +  guides(color = guide_legend(override.aes = list(size=4)) ) 

order.genes <- order.genes[!grepl("ENSMPUG", order.genes)]


plot_pseudotime_heatmap(cds[order.genes,], num_clusters = 4, cores = 1, show_rownames = T, return_heatmap = T)
dev.off()

diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res <- differentialGeneTest(cds[order.genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")



sig_gene_names <- row.names(subset(diff_test_res, qval < 4.345613e-250))


plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 1,
                        cores = 1,
                        show_rownames = T) 


?plot_pseudotime_heatmap

dev.off()

gene_short_name %in% c("LMCD1")

plot_genes_in_pseudotime(cds, color_by = "group")

my_genes <- row.names(subset(fData(cds),
                             gene_short_name %in% c("COL1A1")))
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "group") + theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 16, color = "black"), legend.text = element_text(size = 16), title = element_text(size = 20), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20)) +  guides(color = guide_legend(override.aes = list(size=6), ncol=1)) 


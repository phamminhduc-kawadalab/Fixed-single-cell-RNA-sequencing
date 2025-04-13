# Run Nichnet

install.packages("devtools")
devtools::install_github("saeyslab/nichenetr", force = T)
devtools::install_github("saeyslab/multinichenetr", force = T)
install.packages("glmmTMB",
                 repos = c("https://glmmTMB.github.io/glmmTMB/repos",
                           getOption("repos")["CRAN"]))

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
install.packages("glmmTMB",repos="https://glmmTMB.github.io/glmmTMB/repos",type="binary")

organism = "mouse"
if(organism == "human"){
  lr_network = readRDS("/Users/lethithanhthuy/Desktop/lr_network_human_21122021.rds")
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS("/Users/lethithanhthuy/Desktop/ligand_target_matrix_nsga2r_final.rds")
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
} else if(organism == "mouse"){
  lr_network = readRDS("/Users/lethithanhthuy/Desktop/lr_network_mouse_21122021.rds")
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS("/Users/lethithanhthuy/Desktop/ligand_target_matrix_nsga2r_final_mouse.rds")
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
}

req_subset <- subset(data.combined, idents = c("HSC"))

sce = Seurat::as.SingleCellExperiment(req_subset, assay = "RNA")

sce = alias_to_symbol_SCE(sce, "mouse") %>% makenames_SCE()

SummarizedExperiment::colData(sce)$orig.ident = SummarizedExperiment::colData(sce)$orig.ident %>% make.names()
SummarizedExperiment::colData(sce)$group = SummarizedExperiment::colData(sce)$group %>% make.names()
SummarizedExperiment::colData(sce)$Cell.types = SummarizedExperiment::colData(sce)$Cell.types %>% make.names()


sample_id = "orig.ident"
group_id = "group"
celltype_id = "Cell.types"
covariates = NA
batches = NA

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()


senders_oi = "Mesothelial cell"

receivers_oi = "HSC"

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]

min_cells = 10

abundance_expression_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)

abundance_expression_info$abund_plot_sample

contrasts_oi = c("'TAARE-(CONT+TAA)/2','CONT-(TAA+TAARE)/2','TAA-(CONT+TAARE)/2'")


contrast_tbl = tibble(contrast = 
                        c("TAARE-(CONT+TAA)/2","CONT-(TAA+TAARE)/2","TAA-(CONT+TAARE)/2"), 
                      group = c("TAARE","CONT","TAA"))


DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)

DE_info$celltype_de$de_output_tidy %>% arrange(p_adj) %>% head()

celltype_de = DE_info$celltype_de$de_output_tidy

sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

sender_receiver_de %>% head(20)

logFC_threshold = 0.250
p_val_threshold = 0.05
fraction_cutoff = 0.05

p_val_adj = FALSE

top_n_target = 500

verbose = TRUE
cores_system = 6
n.cores = min(cores_system, sender_receiver_de$receiver %>% unique() %>% length())

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose, 
  n.cores = n.cores
)))

ligand_activities_targets_DEgenes$de_genes_df %>% head(20)

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))

prioritization_tables$group_prioritization_tbl %>% head(20)

lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

path = "./Users/lethithanhthuy/Desktop/multinichnet/"

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

save = FALSE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output.rds"))
  
}

prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 200, rank_per_group = FALSE)

prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)

circos_list$TAA
circos_list$legend

group_oi = "TAA"

prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 200, groups_oi = group_oi)

plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi

library(Seurat)
FeaturePlot(data.combined, "Il18r1", split.by = "group")

group_oi = "TAA"
receiver_oi = "HSC"
lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi, receiver == receiver_oi)
lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50)) # downregulation -- negative correlation
lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)

prioritized_tbl_oi = get_top_n_lr_pairs(prioritization_tables, 50, groups_oi = group_oi, receivers_oi = receiver_oi)

lr_target_correlation_plot = make_lr_target_correlation_plot(multinichenet_output$prioritization_tables, prioritized_tbl_oi,  lr_target_prior_cor_filtered , multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, receiver_oi,plot_legend = FALSE)
lr_target_correlation_plot$combined_plot

ligand_oi = "Wnt5a"
receptor_oi = "Fzd1"
sender_oi = "HSC"
receiver_oi = "HSC"
lr_target_scatter_plot = make_lr_target_scatter_plot(multinichenet_output$prioritization_tables, ligand_oi, receptor_oi, sender_oi, receiver_oi, multinichenet_output$celltype_info, multinichenet_output$grouping_tbl, lr_target_prior_cor_filtered)
lr_target_scatter_plot

combined_plot$combined_plot

write.csv(prioritized_tbl_oi_M_50, "/Users/lethithanhthuy/Desktop/prioritized_tbl_oi_M_50.csv")

data.combined@meta.data

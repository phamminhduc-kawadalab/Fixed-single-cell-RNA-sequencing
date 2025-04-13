# RUN CELLCHAT

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
cellchat <- createCellChat(object = data.combined_cont, group.by = "Cell.types", assay = "RNA")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
cellchat@netP[["pathways"]]
pathways.show <- c("COLLAGEN")
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = seq(1,4), layout = "hierarchy")
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
par(mfrow = c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show <- c("COLLAGEN"), title = "Contribution of each L-R/all")
netAnalysis_contribution(cellchat, signaling = pathways.show <- c(cellchat@netP[["pathways"]]), title = "Contribution of each L-R/all")
pairLR1 <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR1[2,]
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, vertex.receiver = seq(1,4), layout = "hierarchy")
netVisual_bubble(cellchat, sources.use = c(1:11), targets.use = 9, remove.isolate = FALSE, font.size.title = 15)
netVisual_bubble(cellchat, sources.use = c("CAF1", "CAF2", "Endothelial cell", "HCC", "B cell", "T cell", "NK-T cell", "NK cell", "Monocytes", "Plasmocytes"), targets.use = "Endothelial cell", signaling = cellchat@netP[["pathways"]], remove.isolate = T, font.size = 10)
pairLR.use <- extractEnrichedLR(cellchat, signaling = cellchat@netP[["pathways"]])
netVisual_bubble(cellchat, sources.use = c(1:10), targets.use = 9, pairLR.use = pairLR.use, remove.isolate = TRUE)
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}

for (i in 1:nrow(mat)) {  
  mat2 <- matrix(1, nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}
pathways.show <- c("TNFa") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
pairLR.use <- extractEnrichedLR(cellchat, signaling = cellchat@netP[["pathways"]])
netVisual_bubble(cellchat, sources.use = c(1:5), targets.use = c(1:11), pairLR.use = pairLR.use, remove.isolate = TRUE)
netVisual_chord_gene(cellchat, sources.use = c(1:10), targets.use = 9, legend.pos.x = 15)
netVisual_chord_gene(cellchat, sources.use=c(1:11), targets.use = 9, signaling = c("COLLAGEN"), legend.pos.x = 8)
plotGeneExpression(cellchat, signaling = cellchat@netP[["pathways"]])
plotGeneExpression(cellchat, signaling = "COLLAGEN", enriched.only = FALSE)
pathways.show <- cellchat@netP[["pathways"]]
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, signaling = "CD45", width = 8, height = 2.5, font.size = 15, font.size.title = 12)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("COLLAGEN", "FN1"))
gg1 + gg2
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("COLLAGEN", "FN1", "APP", "LAMININ", "COMPLEMENT", "TENASCIN", "CD45", "VTN", "GALECTIN", "SEMA4", "GAS", "PARs", "SPP1", "ANGPTL", "JAM", "PECAM1", "ESAM", "VEGF", "RELN", "CD22", "CDH", "SELPLG", "HSPG", "ICAM", "NOTCH", "VCAM", "AGT", "CSF", "NGF", "ARGN", "EPHA", "SEMA3", "SN", "ITGAL-ITGB2", "PROS", "CADM", "CD52", "PDGF", "FGF", "CD23", "THY1","CCL"), pattern = "outgoing")
ht1
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("COLLAGEN", "FN1", "APP", "LAMININ", "COMPLEMENT", "TENASCIN", "CD45", "VTN", "GALECTIN", "SEMA4", "GAS", "PARs", "SPP1", "ANGPTL", "JAM", "PECAM1", "ESAM", "VEGF", "RELN", "CD22", "CDH", "SELPLG", "HSPG", "ICAM", "NOTCH", "VCAM", "AGT", "CSF", "NGF", "ARGN", "EPHA", "SEMA3", "SN", "ITGAL-ITGB2", "PROS", "CADM", "CD52", "PDGF", "FGF", "CD23", "THY1","CCL"), pattern = "incoming")
ht2
ht1+ht2

gg1 <- netAnalysis_signalingRole_scatter(cellchat)


cont_data <- subset(data.combined, idents = "CONT")
taa_data <- subset(data.combined, idents = "TAA")
taare_data <- subset(data.combined, idents = "TAARE")

unique_identities <- levels(data.combined$Idents)
print(unique_identities)

slotNames(data.combined)

cont_data <- subset(data.combined, subset = group == "CONT")

taa_data <- subset(data.combined, subset = group == "TAA")

taare_data <- subset(data.combined, subset = group == "TAARE")

# Run CellChat for control group
library(CellChat)
cont_data@meta.data

cellchat <- createCellChat(object = cont_data, group.by = "Cell.types", assay = "RNA")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
cellchat@netP[["pathways"]]

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("COLLAGEN", "FN1"))
gg1 + gg2
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("COLLAGEN"))
gg2

# Run CellChat for TAA group

taa_data@meta.data

cellchat_TAA <- createCellChat(object = taa_data, group.by = "Cell.types", assay = "RNA")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat_TAA@DB <- CellChatDB.use
cellchat_TAA <- subsetData(cellchat_TAA)
cellchat_TAA <- identifyOverExpressedGenes(cellchat_TAA)
cellchat_TAA <- identifyOverExpressedInteractions(cellchat_TAA)
cellchat_TAA <- projectData(cellchat_TAA, PPI.mouse)
cellchat_TAA <- computeCommunProb(cellchat_TAA, raw.use = TRUE)
cellchat_TAA <- filterCommunication(cellchat_TAA, min.cells = 10)
cellchat_TAA <- computeCommunProbPathway(cellchat_TAA)
cellchat_TAA <- aggregateNet(cellchat_TAA)
groupSize <- as.numeric(table(cellchat_TAA@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_TAA@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_TAA@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat_TAA@net$weight
cellchat_TAA@netP[["pathways"]]

cellchat_TAA <- netAnalysis_computeCentrality(cellchat_TAA, slot.name = "netP") 

gg1 <- netAnalysis_signalingRole_scatter(cellchat_TAA)
gg2 <- netAnalysis_signalingRole_scatter(cellchat_TAA, signaling = c("COLLAGEN", "FN1"))
gg1 + gg2

# Run CellChat for recovery group

taare_data@meta.data

cellchat_TAARE <- createCellChat(object = taare_data, group.by = "Cell.types", assay = "RNA")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat_TAARE@DB <- CellChatDB.use
cellchat_TAARE <- subsetData(cellchat_TAARE)
cellchat_TAARE <- identifyOverExpressedGenes(cellchat_TAARE)
cellchat_TAARE <- identifyOverExpressedInteractions(cellchat_TAARE)
cellchat_TAARE <- projectData(cellchat_TAARE, PPI.mouse)
cellchat_TAARE <- computeCommunProb(cellchat_TAARE, raw.use = TRUE)
cellchat_TAARE <- filterCommunication(cellchat_TAARE, min.cells = 10)
cellchat_TAARE <- computeCommunProbPathway(cellchat_TAARE)
cellchat_TAARE <- aggregateNet(cellchat_TAARE)
groupSize <- as.numeric(table(cellchat_TAARE@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_TAARE@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_TAARE@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat_TAARE@net$weight
cellchat_TAARE@netP[["pathways"]]

netVisual_bubble(cellchat_TAARE, sources.use = "NK-T cell", targets.use = c("HSC"), remove.isolate = FALSE, font.size.title = 15)

cellchat_TAARE <- netAnalysis_computeCentrality(cellchat_TAARE, slot.name = "netP") 

gg1 <- netAnalysis_signalingRole_scatter(cellchat_TAARE)
gg2 <- netAnalysis_signalingRole_scatter(cellchat_TAARE, signaling = c("COLLAGEN", "FN1"))
gg1 + gg2

gg1 <- netAnalysis_signalingRole_scatter(cellchat_TAA, signaling = c("FN1"))
gg2 <- netAnalysis_signalingRole_scatter(cellchat_TAARE, signaling = c("FN1"))
gg1 + gg2

# Merge cellchat

object.list <- list(TAA = cellchat_TAA, TAARE = cellchat_TAARE)
cellchat_merge <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- rankNet(cellchat_merge, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merge, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

netVisual_bubble(cellchat_merge, sources.use = c("Hepatocytes", "Macrophages","Cholangiocytes", "NK-T cell", "LSEC", "HSC", "Mesothelial cell", "Dendritic cell", "Neutrophils"), targets.use = "HSC",  comparison = c(1, 2, 3), angle.x = 45, font.size.x = 6.7) 

pos.dataset = "TAA"

features.name = pos.dataset

cellchat_merge <- identifyOverExpressedGenes(cellchat_merge, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

net <- netMappingDEG(cellchat_merge, features.name = features.name)

net.up <- subsetCommunication(cellchat_merge, net = net, datasets = "TAA",ligand.logFC = 0.2, receptor.logFC = NULL)

net.down <- subsetCommunication(cellchat_merge, net = net, datasets = "TAARE",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merge)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merge)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_merge, pairLR.use = pairLR.use.up, sources.use = c("Hepatocytes", "Macrophages","Cholangiocytes", "NK-T cell", "LSEC", "HSC", "Mesothelial cell", "Neutrophils", "Dendritic cell"), targets.use = "HSC", comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]), font.size = 10)
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat_merge, pairLR.use = pairLR.use.down, sources.use = c("Hepatocytes", "Macrophages","Cholangiocytes", "NK-T cell", "LSEC", "HSC", "Mesothelial cell", "Neutrophils", "Dendritic cell"), targets.use = "HSC", comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), font.size = 10)

gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c("Hepatocytes", "Macrophages","Cholangiocytes", "NK-T cell", "LSEC", "HSC", ""), targets.use = "HSC", slot.name = 'net', net = net.up, lab.cex = 0.5, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
netVisual_chord_gene(object.list[[1]], sources.use = c("Hepatocytes", "Macrophages","Cholangiocytes", "NK-T cell", "LSEC", "HSC"), targets.use = "HSC", slot.name = 'net', net = net.down, lab.cex = 0.5, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))

future::plan("multisession", workers = 2)

cellchat_merge <- computeNetSimilarityPairwise(cellchat_merge, type = "functional")

cellchat_merge <- netEmbedding(cellchat_merge, type = "functional")

cellchat_merge <- netClustering(cellchat_merge, type = "functional", do.parallel=FALSE)

netVisual_embeddingPairwise(cellchat_merge, type = "functional", label.size = 3.5)


future::plan("multisession", workers = 2)

cellchat_merge <- computeNetSimilarityPairwise(cellchat_merge, type = "structural")

cellchat_merge <- netEmbedding(cellchat_merge, type = "structural")

cellchat_merge <- netClustering(cellchat_merge, type = "structural", do.parallel = FALSE)

netVisual_embeddingPairwise(cellchat_merge, type = "structural", label.size = 3.5)

library(ComplexHeatmap)

i = 1

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 16, font.size = 6.3)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 16, font.size = 6.3)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 16, font.size = 6.3)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 16, font.size = 6.3)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
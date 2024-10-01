### this script will process single cell RNA seq from ADARKO teratoma (n=4), using non-perturbed teratoma as a reference for cell annotation, 
### as well as call perturbation (which KO vs. wild type) for each cell

library(Seurat)
library(cowplot)
library(ggplot2)
library(DoubletFinder)
library(cellMapper)

input_dir = '/media/NAS1/Sammi_NAS1/Cellranger_Output/ADAR_20220614/dual_ref/'
ter_list <- c("ter1","ter2","ter3","ter4")

# -----------------------------------------------------------
# Preprocessing ADAR-KO teratoma, removing low quality cells and doublets

#Load dataset
ter_KO <- list(Read10X(data.dir = paste0(input_dir, "/ter1/outs/filtered_feature_bc_matrix")),
               Read10X(data.dir = paste0(input_dir, "/ter2/outs/filtered_feature_bc_matrix")),
               Read10X(data.dir = paste0(input_dir, "/ter3/outs/filtered_feature_bc_matrix")),
               Read10X(data.dir = paste0(input_dir, "/ter4/outs/filtered_feature_bc_matrix")))
names(ter_KO) <- c("ter1","ter2","ter3","ter4")

# Start Seurat Processing
QC_table <- data.frame(row.names =  names(ter_KO))
QC_table$Cell_Count_from_Cellranger <- sapply(sapply(ter_KO, "[[", 1), dim)[2,]
QC_table

# SeruatObject of RNA datasets
ter_KO.SO.list <- list()
for (i in 1:length(ter_KO)){
  ter_KO.SO.list[i][[1]] <- CreateSeuratObject(counts = ter_KO[i][[1]]$`Gene Expression`, project = "KO_H1",
                                               min.cells = round(0.001*ncol(ter_KO[i][[1]]$`Gene Expression`)), 
                                               min.features = 50)
}

QC_table$SeuObj_min_cell_0.1per_min_features_200 <- sapply(ter_KO.SO.list, dim)[2,]
QC_table

head(ter_KO.SO.list[[1]]@meta.data)
names(ter_KO.SO.list) <- names(ter_KO)


# putting in some metadata
info_exp <- c("ter1","ter2","ter3","ter4")

for (i in 1:length(ter_KO.SO.list)) {
  ter_KO.SO.list[i][[1]][['teratoma']] <- info_exp[i]
  
  ter_KO.SO.list[i][[1]] <- RenameCells(object = ter_KO.SO.list[i][[1]], 
                                        new.names = paste0(info_exp[i], "_",
                                                           colnames(ter_KO.SO.list[i][[1]])))
}

head(ter_KO.SO.list[1][[1]]@meta.data)

rm(info_exp)

# note mito/mouse content cells
for (i in 1:length(ter_KO.SO.list)){
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  ter_KO.SO.list[i][[1]][["percent.mt"]] <- PercentageFeatureSet(ter_KO.SO.list[i][[1]], pattern = "^GRCh38-MT-")
  ter_KO.SO.list[i][[1]][["mito.content"]] <- "High"
  ter_KO.SO.list[i][[1]]@meta.data[which(ter_KO.SO.list[i][[1]][["percent.mt"]]<10),"mito.content"] <- "Low"
}

head(ter_KO.SO.list[1][[1]]@meta.data)

p <- list(list())
for (i in 1:length(ter_KO.SO.list)){
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  p[i][[1]] <- plot_grid(VlnPlot(ter_KO.SO.list[i][[1]], features = c("nFeature_RNA")) + NoLegend(),
                         VlnPlot(ter_KO.SO.list[i][[1]], features = c("nCount_RNA")) + NoLegend(),
                         VlnPlot(ter_KO.SO.list[i][[1]], features = c("percent.mt")) + NoLegend(),
                         ncol = 3, labels = names(ter_KO.SO.list)[i], scale = 0.9)
}
pdf("Figure3_ADARKO_R1/KO_4ter_Vln_postcr_prefilter_v1.pdf", width = 11, height = 6)
p
dev.off()

### Filtering for features
p <- list(list())
for (i in 1:length(ter_KO.SO.list)){
  ter_KO.SO.list[i][[1]] <- subset(ter_KO.SO.list[i][[1]], subset = nFeature_RNA > 200 & 
                                     nFeature_RNA < 7500 & percent.mt < 10) #& percent.mm10 < 10)
  p[i][[1]] <- plot_grid(VlnPlot(ter_KO.SO.list[i][[1]], features = c("nFeature_RNA")) + NoLegend(),
                         VlnPlot(ter_KO.SO.list[i][[1]], features = c("nCount_RNA")) + NoLegend(),
                         VlnPlot(ter_KO.SO.list[i][[1]], features = c("percent.mt")) + NoLegend(),
                         ncol = 3, labels = names(ter_KO.SO.list)[i], scale = 0.9)
}
pdf("Figure3_ADARKO_R1/KO_4ter_Vln_postcr_postfilter_v1.pdf", width = 11, height = 6)
p
dev.off()

QC_table$Count_nFeature_200_7500_mt_10 <- sapply(ter_KO.SO.list, dim)[2,]
QC_table

# preprocessing
for (i in 1:length(ter_KO.SO.list)){
  ter_KO.SO.list[i][[1]] <- NormalizeData(ter_KO.SO.list[i][[1]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  ter_KO.SO.list[i][[1]] <- FindVariableFeatures(ter_KO.SO.list[i][[1]], selection.method = "vst", nfeatures = 4000, verbose = F)
  ter_KO.SO.list[i][[1]] <- ScaleData(ter_KO.SO.list[i][[1]]) # this step could include vars.to.regress
}
hist(log10(ter_KO.SO.list[1][[1]]@assays$RNA@counts@x)) # raw
hist(ter_KO.SO.list[1][[1]]@assays$RNA@data@x) # normalized data
hist(ter_KO.SO.list[1][[1]]@assays$RNA@scale.data) # scaled

p <- list(list())
for (i in 1:length(ter_KO.SO.list)){
  ter_KO.SO.list[i][[1]] <- RunPCA(ter_KO.SO.list[i][[1]], npcs = 30, verbose = F)
  p[i][[1]] <- ElbowPlot(ter_KO.SO.list[i][[1]], ndims = 30)
}
pdf("Figure3_ADARKO_R1/KO_4ter_Elbow_v1.pdf", width = 11, height = 6)
p
dev.off()

for (i in 1:length(ter_KO.SO.list)){
  ter_KO.SO.list[i][[1]] <- RunUMAP(ter_KO.SO.list[i][[1]], dims = 1:30, verbose = F)
  ter_KO.SO.list[i][[1]] <- FindNeighbors(ter_KO.SO.list[i][[1]], k.param = 10, reduction = "pca", dims = 1:30, verbose = F)
  ter_KO.SO.list[i][[1]] <- FindClusters(ter_KO.SO.list[i][[1]], resolution = 0.6, 
                                         algorithm = 1, print.output = F)
}

save(ter_KO.SO.list, file ='Figure3_ADARKO_R1/KO_4ter_preprocessed_for_DF_20230201.Rda')
write.csv(QC_table, "Figure3_ADARKO_R1/KO_4ter_QC_table_20230201.csv", row.names = TRUE, quote = FALSE)

load("Figure3_ADARKO_R1/KO_4ter_preprocessed_for_DF_20230201.Rda")
QC_table <- read.table("Figure3_ADARKO_R1/KO_4ter_QC_table_20230201.csv", sep = ",", header = T)

### doublet finder
ListTest <- ter_KO.SO.list

DF_File_Names <- paste0(ter_list, "_DFed_SO.rds")
npcs = 75

# model for linear fit between multiplate rates & recovered nuclei)
#doubleRateDF <- data.frame(multiplate_rate = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6)/100,
#                           recovered_cell = c(0.5e3, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3))
#lmDR <- lm(multiplate_rate~recovered_cell, data = doubleRateDF)
#estimated_doublet_rate = sapply(ListTest,dim)[2,]*lmDR$coefficients[2]+lmDR$coefficients[1]
estimated_doublet_rate = 0.076 # loaded 10k cells for all datasets
QC_table$estimated_doublet_rate <- rep(estimated_doublet_rate, length(ListTest))

length(ListTest)
for (i in 1:length(ListTest)) {
  seuratTemp <- ListTest[[i]]
  # Run pK indentification
  sweep.res.list_seurat <- paramSweep_v3(seuratTemp, PCs = 1:npcs, sct = T)
  sweep.stats_seurat <- summarizeSweep(sweep.res.list_seurat, GT = F)
  bcmvn <- find.pK(sweep.stats_seurat)
  pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message(pK)
  # Homotypic Doublet Proportion Estimate
  annotations <- seuratTemp@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  # Run DoubletFinder
  nExp_poi <- round(estimated_doublet_rate * nrow(seuratTemp@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seuratTemp_DF <- doubletFinder_v3(seuratTemp, PCs = 1:npcs, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
  # Check # of Doublets and Save
  head(seuratTemp_DF@meta.data)
  table(seuratTemp_DF@meta.data[,dim(seuratTemp_DF@meta.data)[2]-1])
  
  saveRDS(seuratTemp_DF, file = paste0("Figure3_ADARKO_R1/preprocess_20230201/DF/",DF_File_Names[i]))
}

# Load DFed datasets
DF_File_Names <- paste0(ter_list, "_DFed_SO.rds")

KO <- list()
for (i in 1:length(ter_list)){
  KO[[i]] <- readRDS(paste0("Figure3_ADARKO_R1/preprocess_20230201/DF/", DF_File_Names[[i]]))
}
names(KO) <- ter_list

rm(ListTest)
rm(sweep.res.list_seurat,sweep.stats_seurat,seuratTemp,seuratTemp_DF,bcmvn)
rm(npcs,pK,nExp_poi,nExp_poi.adj)
gc()

# brute force modifying the colnames of DF results
for (i in 1:length(KO)) {
  colnames(KO[i][[1]]@meta.data)[ncol(KO[i][[1]]@meta.data)-1] <- "DF_pANN"
  colnames(KO[i][[1]]@meta.data)[ncol(KO[i][[1]]@meta.data)] <- "DF_Classification"
  QC_table$DF_Singlet[i] <- length(which(KO[i][[1]][["DF_Classification"]] == "Singlet"))
}
QC_table
write.csv(QC_table, "Figure3_ADARKO_R1/KO_4ter_QC_table_20230201.csv", row.names = TRUE, quote = FALSE)

save(KO, file = 'Figure3_ADARKO_R1/preprocess_20230201/KO_4ter_DFed_Seurat_list_20230201.rda')
load('Figure3_ADARKO_R1/preprocess_20230201/KO_4ter_DFed_Seurat_list_20230201.rda')

### distinguish the human and the mouse cells using cellranger output
cn_prefix <- c("ter1", "ter2", "ter3", "ter4")
ter_gem_cla <- list()
for (i in 1:length(ter_list)){
  ter_gem_cla[[i]] <- read.table(paste0(input_dir, ter_list[[i]], "/outs/analysis/gem_classification.csv"), sep = ",", header = TRUE)
  names(ter_gem_cla)[[i]] <- cn_prefix[[i]]
  row.names(ter_gem_cla[[i]]) <- paste0(names(ter_gem_cla)[[i]], "_", ter_gem_cla[[i]]$barcode)
}
length(ter_gem_cla)
head(ter_gem_cla[[1]],6)

QC_table$human_count_cr <- sapply(ter_gem_cla, function(x) table(x$call))["GRCh38",]
QC_table$mouse_count_cr <- sapply(ter_gem_cla, function(x) table(x$call))["mm10",]
QC_table$multiplet_count_cr <- sapply(ter_gem_cla, function(x) table(x$call))["Multiplet",]
QC_table

for (i in 1:length(ter_list)){
  cr_call_temp <- ter_gem_cla[[i]]$call
  names(cr_call_temp) <- row.names(ter_gem_cla[[i]])
  KO[[i]] <- AddMetaData(KO[[i]], metadata = cr_call_temp, col.name = "cr_call")
  print(table(KO[[i]]$cr_call))
  print(length(which(KO[[i]]$cr_call == "GRCh38" & KO[[i]]$DF_Classification == "Singlet")))
}

table(KO[[2]]$cr_call)

### Filtering for Singlet Human cells
for (i in 1:length(KO)) {
  Idents(KO[i][[1]]) <- KO[i][[1]]$DF_Classification
  KO[i][[1]] <- subset(KO[i][[1]], ident = "Singlet")
  
  Idents(KO[i][[1]]) <- KO[i][[1]]$cr_call
  KO[i][[1]] <- subset(KO[i][[1]], ident = "GRCh38")
  
  Idents(KO[i][[1]]) <- KO[i][[1]]$orig.ident
}
QC_table$Singlet_Human_Cells <- sapply(KO, dim)[2,]
sum(sapply(KO, dim)[2,]) # 33710

### generate dataset that only contain human genes
KO.human <- list()
for (i in 1:length(KO)){
  genes.to.keep <- grep(pattern = "^GRCh38-", x = row.names(x = KO[i][[1]]), value = TRUE)
  KO[i][[1]] <- subset(KO[i][[1]], features = genes.to.keep)
  
  counts.temp <- GetAssayData(KO[i][[1]], assay = "RNA", slot = "counts")
  counts.temp@Dimnames[[1]] <- sub('GRCh38-', '', counts.temp@Dimnames[[1]])
  
  KO.human[i][[1]] <- CreateSeuratObject(counts.temp,
                                                       meta.data = KO[i][[1]]@meta.data)
  rm(counts.temp,genes.to.keep)
}
row.names(KO.human[i][[1]])
QC_table$human_genes_to_keep <- sapply(KO.human, dim)[1,]
QC_table

KO <- KO.human
rm(KO.human)

# Remove noncoding RNAs and mitochondrial genes (RP pseudogenes not functional)
KO.cleanup <- list()
for (i in 1:length(KO)){
  genes.keep <- rownames(KO[i][[1]])[!grepl("^RP11|^MT-", rownames(KO[i][[1]]))]
  KO[i][[1]] <- subset(KO[i][[1]], features = genes.keep)
  
  counts.temp <- GetAssayData(KO[i][[1]], assay = "RNA", slot = "counts")
  KO.cleanup[i][[1]] <- CreateSeuratObject(counts.temp,
                                                         meta.data = KO[i][[1]]@meta.data)
  rm(counts.temp,genes.keep)
}
row.names(KO.cleanup[i][[1]])
QC_table$relevant_genes_to_keep <- sapply(KO.cleanup, dim)[1,]
QC_table

KO <- KO.cleanup
rm(KO.cleanup)

names(KO) <- ter_list
KO

head(KO[1][[1]]@meta.data)

grep(pattern = "^mm10-", x = row.names(x = KO[i][[1]]), value = TRUE)
grep(pattern = "^GRCh38-", x = row.names(x = KO[i][[1]]), value = TRUE)
grep(pattern = "^MT-", x = row.names(x = KO[i][[1]]), value = TRUE)

save(KO, file = 'Figure3_ADARKO_R1/KO_4ter_DFed_Seurat_humanCells_list_20230201.rda')
load('Figure3_ADARKO_R1/KO_4ter_DFed_Seurat_humanCells_list_20230201.rda')

# -----------------------------------------------------------
# map QC'ed cells in each ADAR-KO experiment onto the reference
load('/media/Scratch_SSD_Voyager/sammi/RNA_editing/WT4H1_YWsubset_integratedObj_umapmodel_20230115.Robj')

anchors <- list()
predictions <- list()
KO.query <- list()
for (i in 1:length(KO)){
  KO[[i]] <- NormalizeData(object = KO[[i]], verbose = FALSE)
  anchors[[i]] <- FindTransferAnchors(reference = WT_4H1, query = KO[[i]],
                                      dims = 1:30, reference.reduction = "pca")#, features = intersect(row.names(terA.human), row.names(KO[[1]])))
  predictions[[i]] <- TransferData(anchorset = anchors[[i]], refdata = WT_4H1$cluster,
                                   dims = 1:30)
  KO.query[[i]] <- KO[[i]]
  KO.query[[i]] <- AddMetaData(KO[[i]], metadata = predictions[[i]])
  KO.query[[i]] <- NormalizeData(KO.query[[i]], verbose = FALSE)
  KO.query[[i]] <- FindVariableFeatures(KO.query[[i]], selection.method = "vst", nfeatures = 2000,
                                                      verbose = FALSE)
  KO.query[[i]] <- MapQuery(anchorset = anchors[[i]], reference = WT_4H1, query = KO.query[[i]],
                                          refdata = list(celltype = "cluster"), reference.reduction = "pca", reduction.model = "umap", verbose = FALSE) 
}
sum(sapply(KO.query, dim)[2,])
# 33710

hist(predictions[[4]]$prediction.score.max)
DimPlot(KO.query[[1]], group.by = "predicted.celltype")
table(KO.query[[1]]$predicted.celltype)
table(KO.query[[2]]$predicted.celltype)
table(KO.query[[3]]$predicted.celltype)
table(KO.query[[4]]$predicted.celltype)
names(KO.query) <- names(KO)
min(KO.query[[1]]$prediction.score.max)
min(KO.query[[2]]$prediction.score.max)
min(KO.query[[3]]$prediction.score.max)
min(KO.query[[4]]$prediction.score.max)


save(KO.query, file = "Figure3_ADARKO_R1/KO_4ter_project_SO_list_20230201.rda")
load("Figure3_ADARKO_R1/KO_4ter_project_SO_list_20230201.rda")

KO <- merge(KO.query[[1]], y = c(KO.query[[2]],KO.query[[3]],KO.query[[4]]), merge.data = T, merge.dr = c("ref.pca","ref.umap"))

KO <- FindNeighbors(KO, k.param = 10, reduction = "ref.pca", dims = 1:30, verbose = F)
KO <- FindClusters(KO, resolution = 0.6, algorithm = 2, print.output = F)

pdf("Figure3_ADARKO_R1/KO_4ter_project_refUMAP_QC_20230201.pdf", height = 12, width = 10)
plot_grid(
  VlnPlot(KO, features = c("nCount_RNA","nFeature_RNA","percent.mt"), group.by = "RNA_snn_res.0.6", ncol = 1, pt.size = -1),
  DimPlot(KO, group.by = "RNA_snn_res.0.6", label = T, repel = T) + NoLegend()
)
dev.off()

plot_grid(
  DimPlot(KO, group.by = "teratoma", label = T, repel = T) + NoLegend(),
  DimPlot(KO, group.by = "predicted.id", label = T, repel = T) + NoLegend()
)

cluster_levels <- c("Radial Glia","CycProg","Early Neurons","Retinal Neurons","Retinal Epi","Schwann Cells","Melanoblasts",
                    "Foregut Epi","Mid/Hindgut Epi","Airway Epi","HSC","Immune","Erythrocyte",
                    "Adipogenic MSC/Fib","Chondrogenic MSC/Fib","Cycling MSC/Fib","MSC/Fib","MyoFib",
                    "Muscle Prog","Cardiac/Skeletal Muscle","Pericytes","Smooth Muscle",
                    "Kidney Prog")
cluster_col <- c("#4795F5","#0B62CD","#7474FF","#DADAFF","#F2B77C","#9EEA00","#4FF4A2",
                 "#002700","#005A00","#00A700","#000032","#04409B","#4545FF",
                 "#0000FF","#FF0000","#00C000","#AD07E3","#000000",
                 "#B35A00","#FF9933","#A00000","#94641F","#610051")
names(cluster_col) <- cluster_levels

pdf("Figure3_ADARKO_R1/KO_4ter_project_refUMAP_20230201.pdf")
DimPlot(KO, group.by = "predicted.id", label = T, repel = T, cols = cluster_col) + NoLegend()
dev.off()

pdf("Figure3_ADARKO_R1/KO_4ter_project_refUMAP_postGC_20230212.pdf", height = 8, width = 10)
DimPlot(KO_sub, group.by = "predicted.id", label = T, repel = T, cols = cluster_col) + NoLegend()
dev.off()





### get a correlation plot with the reference that I am using

YanMarkersPub <- c("VIM","SOX2","HES5","HMGB2","DCX","MAP2",
                   "OTX2","NRL","MITF","TTR",
                   "MPZ","SOX10","MLANA",
                   "ELF3","PAX9","KRT4","FOXJ1","CDHR3","CDX2",
                   "CD74","CD34","HHEX","THY1","ITM2A","SHOX2","COL2A1","SOX9","COL15A1",
                   "MYOD1","PAX7","TNNT1","TNNI2","FOXC1","CYP1B1","ACTA2","RGS5","WT1")
KO <- FindVariableFeatures(KO, selection.method = "vst", nfeatures = 4000, verbose = FALSE)

genes.use <- unique(c(as.character(YanMarkersPub), VariableFeatures(WT_4H1), VariableFeatures(KO)))
#genes.use <- Reduce(intersect, list(genes.use, rownames(GetAssayData(WT_4H1, assay = "integrated", slot = "scale.data")), 
#                                    rownames(GetAssayData(KO, assay = "RNA"))))
length(genes.use) # 3994
Tera.r.assay <- GetAssayData(KO, assay = "RNA")
cor.result <- MapClustersCor(query.data = Tera.r.assay, query.clusters = KO$predicted.id, 
                             train.data = GetAssayData(WT_4H1, assay = "integrated", slot = "scale.data"),
                             train.clusters = WT_4H1$cluster,
                             genes.use = genes.use)
pdf("Figure3_ADARKO_R1/KO_4ter_project_correlation_with_WT_20230201.pdf",height = 8, width = 8.5)
ggHeat(t(cor.result$cluster.correlations), x.lab.size = 11, y.lab.size = 11) + ggtitle("Row: Wild Type Cell Types, Column: ADAR KO Cell Types")
dev.off()

save(KO, file = "Figure3_ADARKO_R1/KO_4ter_project_SO_annot_20230201.rda")

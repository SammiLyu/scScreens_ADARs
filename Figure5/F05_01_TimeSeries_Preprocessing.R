### script pre-process H1 Teratoma Time Series snRNAseq data
### doublet identification method: DoubletFinder
### Annotation method: scanvi, manual
### Integration method: log normalized seurat integration

### following Blue's advice from presentation 20230118, just transfer all labels using SCANVI

library(Seurat)
library(ggplot2)
library(cowplot)
library(DoubletFinder)
library(SeuratDisk)
library(cellMapper)
library(swne)
library(RColorBrewer)
library(stringr)

setwd('/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/')

# define markers
YanMarkersPub <- c("VIM","SOX2","HES5","HMGB2","DCX","MAP2",
                   "OTX2","NRL","MITF","TTR",
                   "MPZ","SOX10","MLANA",
                   "ELF3","PAX9","KRT4","FOXJ1","CDHR3","CDX2",
                   "CD74","CD34","HHEX","THY1","ITM2A","SHOX2","COL2A1","SOX9","COL15A1",
                   "MYOD1","PAX7","TNNT1","TNNI2","FOXC1","CYP1B1","ACTA2","RGS5","WT1")

# -----------------------------------------------------------
# Preprocessing RNA seurat object, removing low quality cells and doublets

### load RNA
input_dir = '/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/Teratoma_TS_RNA_intronT/multispecies/'
ter_list <- c("wk4_M1_1","wk6_M15_1","wk8_M8_1","wk10_A_2","wk4_M14_1","wk6_M9_1","wk8_M10_1","wk10_B_1")

TS.RNA.list <- list()
for (ter in ter_list){
  TS.RNA.list[[ter]] <- Read10X(data.dir = paste0(input_dir, ter, "/outs/filtered_feature_bc_matrix"))
}

rm(input_dir)

TS.RNA.list
length(TS.RNA.list)

QC_table <- data.frame(row.names =  names(TS.RNA.list))
QC_table$Cell_Count_from_Cellranger <- sapply(TS.RNA.list, dim)[2,]
QC_table

# SeruatObject of RNA datasets
TS.RNA.SO.list <- list()
for (i in 1:length(TS.RNA.list)){
  TS.RNA.SO.list[i][[1]] <- CreateSeuratObject(counts = TS.RNA.list[[i]], project = "TS_RNA",
                                               min.cells = round(0.001*ncol(TS.RNA.list[[i]])), 
                                               min.features = 200)
}

QC_table$SeuObj_min_cell_0.1per_min_features_200 <- sapply(TS.RNA.SO.list, dim)[2,]
QC_table

head(TS.RNA.SO.list[[1]]@meta.data)
names(TS.RNA.SO.list) <- names(TS.RNA.list)

# putting in some metadata
info_exp <- c("20210205","20210205","20210303","20210714","20210714","20210714","20210721","20210721")
info_age <- c(4,6,8,10,4,6,8,10)
info_mouse <- c('M1','M15','M8','A','M14','M9','M10','B')
info_replicate <- c(1,1,1,2,1,1,1,1)

for (i in 1:length(TS.RNA.SO.list)) {
  TS.RNA.SO.list[i][[1]][['teratoma']] <- info_exp[i]
  TS.RNA.SO.list[i][[1]][['age']] <- info_age[i]
  TS.RNA.SO.list[i][[1]][['mouse']] <- info_mouse[i]
  TS.RNA.SO.list[i][[1]][['replicate']] <- info_replicate[i]
  
  TS.RNA.SO.list[i][[1]] <- RenameCells(object = TS.RNA.SO.list[i][[1]], 
                                        new.names = paste0(info_exp[i], "_wk",
                                                           info_age[i], "_",
                                                           info_mouse[i],"_",
                                                           info_replicate[i],"_", 
                                                           colnames(TS.RNA.SO.list[i][[1]])))
}

head(TS.RNA.SO.list[1][[1]]@meta.data)

rm(info_age,info_exp,info_mouse,info_replicate)

# note mito/mouse content cells
for (i in 1:length(TS.RNA.SO.list)){
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  TS.RNA.SO.list[i][[1]][["percent.mt"]] <- PercentageFeatureSet(TS.RNA.SO.list[i][[1]], pattern = "^GRCh38-MT-")
  TS.RNA.SO.list[i][[1]][["mito.content"]] <- "High"
  TS.RNA.SO.list[i][[1]]@meta.data[which(TS.RNA.SO.list[i][[1]][["percent.mt"]]<10),"mito.content"] <- "Low"
}

p <- list(list())
for (i in 1:length(TS.RNA.SO.list)){
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  p[i][[1]] <- plot_grid(VlnPlot(TS.RNA.SO.list[i][[1]], features = c("nFeature_RNA")) + NoLegend(),
                         VlnPlot(TS.RNA.SO.list[i][[1]], features = c("nCount_RNA")) + NoLegend(),
                         VlnPlot(TS.RNA.SO.list[i][[1]], features = c("percent.mt")) + NoLegend(),
                         ncol = 3, labels = names(TS.RNA.SO.list)[i], scale = 0.9)
}
pdf("Figure/preprocess_20221002/TS_8RNA_Vln_postcr_prefilter_v1.pdf", width = 11, height = 6)
p
dev.off()


### Filtering for features
p <- list(list())
for (i in 1:length(TS.RNA.SO.list)){
  TS.RNA.SO.list[i][[1]] <- subset(TS.RNA.SO.list[i][[1]], subset = nFeature_RNA > 200 & 
                                     nFeature_RNA < 5000 & percent.mt < 10) #& percent.mm10 < 10)
  p[i][[1]] <- plot_grid(VlnPlot(TS.RNA.SO.list[i][[1]], features = c("nFeature_RNA")) + NoLegend(),
                         VlnPlot(TS.RNA.SO.list[i][[1]], features = c("nCount_RNA")) + NoLegend(),
                         VlnPlot(TS.RNA.SO.list[i][[1]], features = c("percent.mt")) + NoLegend(),
                         ncol = 3, labels = names(TS.RNA.SO.list)[i], scale = 0.9)
}
pdf("Figure/preprocess_20221002/TS_8RNA_Vln_postcr_postfilter_v1.pdf", width = 11, height = 6)
p
dev.off()

QC_table$Count_nFeature_200_5000_mt_10 <- sapply(TS.RNA.SO.list, dim)[2,]
QC_table

# preprocessing
for (i in 1:length(TS.RNA.SO.list)){
  TS.RNA.SO.list[i][[1]] <- NormalizeData(TS.RNA.SO.list[i][[1]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  TS.RNA.SO.list[i][[1]] <- FindVariableFeatures(TS.RNA.SO.list[i][[1]], selection.method = "vst", nfeatures = 4000, verbose = F)
  TS.RNA.SO.list[i][[1]] <- ScaleData(TS.RNA.SO.list[i][[1]]) # this step could include vars.to.regress
}
hist(log10(TS.RNA.SO.list[1][[1]]@assays$RNA@counts@x)) # raw
hist(TS.RNA.SO.list[1][[1]]@assays$RNA@data@x) # normalized data
hist(TS.RNA.SO.list[1][[1]]@assays$RNA@scale.data) # scaled

p <- list(list())
for (i in 1:length(TS.RNA.SO.list)){
  TS.RNA.SO.list[i][[1]] <- RunPCA(TS.RNA.SO.list[i][[1]], npcs = 30, verbose = F)
  p[i][[1]] <- ElbowPlot(TS.RNA.SO.list[i][[1]], ndims = 30)
}
pdf("Figure/preprocess_20221002/TS_8RNA_Elbow_v1.pdf", width = 11, height = 6)
p
dev.off()

# seurat moved umap from reticulate to uwot - shouldn't affect results too much
# https://jlmelville.github.io/uwot/pycompare
for (i in 1:length(TS.RNA.SO.list)){
  TS.RNA.SO.list[i][[1]] <- RunUMAP(TS.RNA.SO.list[i][[1]], dims = 1:30, verbose = F)
  TS.RNA.SO.list[i][[1]] <- FindNeighbors(TS.RNA.SO.list[i][[1]], k.param = 10, reduction = "pca", dims = 1:30, verbose = F)
  TS.RNA.SO.list[i][[1]] <- FindClusters(TS.RNA.SO.list[i][[1]], resolution = 0.6, 
                                         algorithm = 1, print.output = F)
}

DimPlot(TS.RNA.SO.list[4][[1]])
FeaturePlot(TS.RNA.SO.list[2][[1]], features = "percent.mt")
#FeaturePlot(TS.RNA.SO.list[2][[1]], features = "percent.mm10")

save(TS.RNA.SO.list, file ='SeuratObj/preprocess_20221002/TS_8RNA_preprocessed_for_DF_20221002.Rda')
write.csv(QC_table, "TS_8RNA_QC_table_20221002.csv", row.names = TRUE, quote = FALSE)

### doublet finder
ListTest <- TS.RNA.SO.list

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
  
  saveRDS(seuratTemp_DF, file = paste0("SeuratObj/preprocess_20221002/DF/",DF_File_Names[i]))
}

# Load DFed datasets
DF_File_Names <- paste0(ter_list, "_DFed_SO.rds")

seuratObj.TS.RNA <- list()
for (i in 1:length(ter_list)){
  seuratObj.TS.RNA[[i]] <- readRDS(paste0("SeuratObj/preprocess_20221002/DF/", DF_File_Names[[i]]))
}

rm(ListTest)
rm(sweep.res.list_seurat,sweep.stats_seurat,seuratTemp,seuratTemp_DF,bcmvn)
rm(npcs,pK,nExp_poi,nExp_poi.adj)
gc()

# brute force modifying the colnames of DF results
for (i in 1:length(seuratObj.TS.RNA)) {
  colnames(seuratObj.TS.RNA[i][[1]]@meta.data)[ncol(seuratObj.TS.RNA[i][[1]]@meta.data)-1] <- "DF_pANN"
  colnames(seuratObj.TS.RNA[i][[1]]@meta.data)[ncol(seuratObj.TS.RNA[i][[1]]@meta.data)] <- "DF_Classification"
  QC_table$DF_Singlet[i] <- length(which(seuratObj.TS.RNA[i][[1]][["DF_Classification"]] == "Singlet"))
}
QC_table
write.csv(QC_table, "TS_8RNA_QC_table_20221002.csv", row.names = TRUE, quote = FALSE)

save(seuratObj.TS.RNA, file = 'SeuratObj/preprocess_20221002/TS_8RNA_DFed_Seurat_list_20221002.rda')
load('SeuratObj/preprocess_20221002/TS_8RNA_DFed_Seurat_list_20221002.rda')

### distinguish the human and the mouse cells using cellranger output

info_exp <- c("20210205","20210205","20210303","20210714","20210714","20210714","20210721","20210721")
info_age <- paste0("wk", c(4,6,8,10,4,6,8,10))
info_mouse <- c('M1','M15','M8','A','M14','M9','M10','B')
info_replicate <- c(1,1,1,2,1,1,1,1)

cn_prefix <- c("20210205_wk4_M1_1", "20210205_wk6_M15_1", "20210303_wk8_M8_1", "20210714_wk10_A_2",
               "20210714_wk4_M14_1","20210714_wk6_M9_1","20210721_wk8_M10_1","20210721_wk10_B_1")
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
  seuratObj.TS.RNA[[i]] <- AddMetaData(seuratObj.TS.RNA[[i]], metadata = cr_call_temp, col.name = "cr_call")
  print(table(seuratObj.TS.RNA[[i]]$cr_call))
  print(length(which(seuratObj.TS.RNA[[i]]$cr_call == "GRCh38" & seuratObj.TS.RNA[[i]]$DF_Classification == "Singlet")))
}

table(seuratObj.TS.RNA[[1]]$cr_call)

### Filtering for Singlet Human cells
for (i in 1:length(seuratObj.TS.RNA)) {
  Idents(seuratObj.TS.RNA[i][[1]]) <- seuratObj.TS.RNA[i][[1]]$DF_Classification
  seuratObj.TS.RNA[i][[1]] <- subset(seuratObj.TS.RNA[i][[1]], ident = "Singlet")
  
  Idents(seuratObj.TS.RNA[i][[1]]) <- seuratObj.TS.RNA[i][[1]]$cr_call
  seuratObj.TS.RNA[i][[1]] <- subset(seuratObj.TS.RNA[i][[1]], ident = "GRCh38")
  
  Idents(seuratObj.TS.RNA[i][[1]]) <- seuratObj.TS.RNA[i][[1]]$orig.ident
}
QC_table$Singlet_Human_Cells <- sapply(seuratObj.TS.RNA, dim)[2,]
sum(sapply(seuratObj.TS.RNA, dim)[2,]) # 64646

### generate dataset that only contain human genes
seuratObj.TS.RNA.human <- list()
for (i in 1:length(seuratObj.TS.RNA)){
  genes.to.keep <- grep(pattern = "^GRCh38-", x = row.names(x = seuratObj.TS.RNA[i][[1]]), value = TRUE)
  seuratObj.TS.RNA[i][[1]] <- subset(seuratObj.TS.RNA[i][[1]], features = genes.to.keep)
  
  counts.temp <- GetAssayData(seuratObj.TS.RNA[i][[1]], assay = "RNA", slot = "counts")
  counts.temp@Dimnames[[1]] <- sub('GRCh38-', '', counts.temp@Dimnames[[1]])
  
  seuratObj.TS.RNA.human[i][[1]] <- CreateSeuratObject(counts.temp,
                                                       meta.data = seuratObj.TS.RNA[i][[1]]@meta.data)
  rm(counts.temp,genes.to.keep)
}
row.names(seuratObj.TS.RNA.human[i][[1]])
QC_table$human_genes_to_keep <- sapply(seuratObj.TS.RNA.human, dim)[1,]
QC_table

seuratObj.TS.RNA <- seuratObj.TS.RNA.human
rm(seuratObj.TS.RNA.human)

# Remove noncoding RNAs and mitochondrial genes (RP pseudogenes not functional)
seuratObj.TS.RNA.cleanup <- list()
for (i in 1:length(seuratObj.TS.RNA)){
  genes.keep <- rownames(seuratObj.TS.RNA[i][[1]])[!grepl("^RP11|^MT-", rownames(seuratObj.TS.RNA[i][[1]]))]
  seuratObj.TS.RNA[i][[1]] <- subset(seuratObj.TS.RNA[i][[1]], features = genes.keep)
  
  counts.temp <- GetAssayData(seuratObj.TS.RNA[i][[1]], assay = "RNA", slot = "counts")
  seuratObj.TS.RNA.cleanup[i][[1]] <- CreateSeuratObject(counts.temp,
                                                         meta.data = seuratObj.TS.RNA[i][[1]]@meta.data)
  rm(counts.temp,genes.keep)
}
row.names(seuratObj.TS.RNA.cleanup[i][[1]])
QC_table$relavant_genes_to_keep <- sapply(seuratObj.TS.RNA.cleanup, dim)[1,]
QC_table

seuratObj.TS.RNA <- seuratObj.TS.RNA.cleanup
rm(seuratObj.TS.RNA.cleanup)

names(seuratObj.TS.RNA) <- cn_prefix
seuratObj.TS.RNA

head(seuratObj.TS.RNA[1][[1]]@meta.data)

grep(pattern = "^mm10-", x = row.names(x = seuratObj.TS.RNA[i][[1]]), value = TRUE)
grep(pattern = "^GRCh38-", x = row.names(x = seuratObj.TS.RNA[i][[1]]), value = TRUE)
grep(pattern = "^MT-", x = row.names(x = seuratObj.TS.RNA[i][[1]]), value = TRUE)


save(seuratObj.TS.RNA, file = 'SeuratObj/key_output/TS_8RNA_DFed_Seurat_humanCells_list_20230115.rda')
load('SeuratObj/key_output/TS_8RNA_DFed_Seurat_humanCells_list_20230115.rda')


# -----------------------------------------------------------
# Preprocessing Multiome RNA seurat object, removing low quality cells and doublets

### load Multiome RNA
input.hg <- list(Read10X_h5("/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/wk4_M1_2/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/wk6_M15_2/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/wk8_M8_2/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/wk10_A_3/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/wk10_A_4/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/wk10_B_3/outs/filtered_feature_bc_matrix.h5"))

input.mm <- list(Read10X_h5("/media/NAS2/Users/Sammi_NAS2/cellranger_Output/MLT_20220103_NovaSeq/mm10ref/wk4_M1_2/MLT_20211104_wk4_M1_2/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS2/Users/Sammi_NAS2/cellranger_Output/MLT_20220103_NovaSeq/mm10ref/wk6_M15_2/MLT_20211104_wk6_M15_2/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS2/Users/Sammi_NAS2/cellranger_Output/MLT_20220103_NovaSeq/mm10ref/wk8_M8_2/MLT_20211104_wk8_M8_2/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS2/Users/Sammi_NAS2/cellranger_Output/MLT_20220103_NovaSeq/mm10ref/wk10_A_3/MLT_20211104_wk10_A_3/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS2/Users/Sammi_NAS2/cellranger_Output/MLT_20220103_NovaSeq/mm10ref/MLT_20211115_wk10_A_4/outs/filtered_feature_bc_matrix.h5"),
                 Read10X_h5("/media/NAS2/Users/Sammi_NAS2/cellranger_Output/MLT_20220103_NovaSeq/mm10ref/MLT_20211115_wk10_B_3/outs/filtered_feature_bc_matrix.h5"))

so.hg <- list()
for (i in 1:length(input.hg)){
  so.hg[i][[1]] <- CreateSeuratObject(counts = input.hg[i][[1]]$`Gene Expression`, project = 'TS_Multiome')
  so.hg[i][[1]][["percent.mt"]] <- PercentageFeatureSet(so.hg[i][[1]], pattern = "^MT-")
  so.hg[i][[1]][["percent.mm10"]] <- 0
  
  for (j in 1:dim(so.hg[i][[1]])[2]){
    if (colnames(so.hg[i][[1]])[j] %in% colnames(input.mm[i][[1]]$`Gene Expression`)){
      mmCounts <- sum(input.mm[i][[1]]$`Gene Expression`[,colnames(so.hg[i][[1]])[j]])
      hgCounts <- sum(input.hg[i][[1]]$`Gene Expression`[,colnames(so.hg[i][[1]])[j]])
      so.hg[i][[1]]$percent.mm10[j] = mmCounts/(mmCounts+hgCounts)*100
    }
  }
}
head(so.hg[1][[1]]@meta.data)
pdf("illustration_multiome_hg_vs_mm.pdf")
hist(so.hg[1][[1]]$percent.mm10)
dev.off()

save(so.hg, file = 'SeuratObj/preprocess_20221002/TS_6Multi_pre-preprocessed_for_DF_20230116.Rda')
load("SeuratObj/preprocess_20221002/TS_6Multi_pre-preprocessed_for_DF_20230116.Rda")

file_name <- paste0("SeuratObj/preprocess_20221002/DF/", ter_multi_names)
npcs = 75
estimated_doublet_rate = 0.076

SO_Process <- function(seuratObj, file_name, estimated_doublet_rate){
  print("Creating Seurat Object from Counts...")
  seuratObj <- CreateSeuratObject(counts = seuratObj@assays$RNA@counts, project = "TS_multiome", min.cells = round(0.001*ncol(seuratObj)), min.features = 200,
                                  meta.data = seuratObj@meta.data)
  
  print("Real quick preprocess...")
  seuratObj <- subset(seuratObj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10) #& percent.mm10 < 10)
  seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 4000, verbose = F)
  seuratObj <- ScaleData(seuratObj) # this step could include vars.to.regress
  seuratObj <- RunPCA(seuratObj, npcs = 30, verbose = F)
  seuratObj <- RunUMAP(seuratObj, dims = 1:30, verbose = F)
  seuratObj <- FindNeighbors(seuratObj, k.param = 10, reduction = "pca", dims = 1:30, verbose = F)
  seuratObj <- FindClusters(seuratObj, resolution = 0.6, algorithm = 1, print.output = F)
  
  print("Running Doublet Finder...")
  # Run pK indentification
  sweep.res.list_seurat <- paramSweep_v3(seuratObj, PCs = 1:npcs, sct = T)
  sweep.stats_seurat <- summarizeSweep(sweep.res.list_seurat, GT = F)
  bcmvn <- find.pK(sweep.stats_seurat)
  pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message(pK)
  # Homotypic Doublet Proportion Estimate
  annotations <- seuratObj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  # Run DoubletFinder
  nExp_poi <- round(estimated_doublet_rate * nrow(seuratObj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seuratTemp_DF <- doubletFinder_v3(seuratObj, PCs = 1:npcs, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
  # Check # of Doublets and Save
  head(seuratTemp_DF@meta.data)
  table(seuratTemp_DF@meta.data[,dim(seuratTemp_DF@meta.data)[2]-1])
  
  saveRDS(seuratTemp_DF, file = paste0(file_name, "_DFed_SO.rds"))
}

for (i in 3:length(so.hg)){
  SO_Process(so.hg[i][[1]], file_name[i], estimated_doublet_rate)
}

seuratObj.TS.multi <- list()
for (i in 1:length(ter_multi_names)){
  seuratObj.TS.multi[[i]] <- readRDS(paste0(file_name[i],"_DFed_SO.rds"))
}

### add metadata
for (i in 1:length(seuratObj.TS.multi)){
  colnames(seuratObj.TS.multi[[i]]@meta.data)[ncol(seuratObj.TS.multi[[i]]@meta.data)-1] <- "DF_pANN"
  colnames(seuratObj.TS.multi[[i]]@meta.data)[ncol(seuratObj.TS.multi[[i]]@meta.data)] <- "DF_Classification"
  
  Idents(seuratObj.TS.multi[[i]]) <- seuratObj.TS.multi[[i]]$DF_Classification
  seuratObj.TS.multi[[i]] <- subset(seuratObj.TS.multi[[i]], ident = "Singlet")
  Idents(seuratObj.TS.multi[[i]]) <-seuratObj.TS.multi[[i]]$orig.ident
}

info_exp <- c("20211104","20211104","20211104","20211104","20211115","20211115")
info_age <- c(4,6,8,10,10,10)
info_mouse <- c('M1','M15','M8','A','A','B')
info_replicate <- c(2,2,2,3,4,3)

for (i in 1:length(seuratObj.TS.multi)) {
  seuratObj.TS.multi[i][[1]][['teratoma']] <- info_exp[i]
  seuratObj.TS.multi[i][[1]][['age']] <- info_age[i]
  seuratObj.TS.multi[i][[1]][['mouse']] <- info_mouse[i]
  seuratObj.TS.multi[i][[1]][['replicate']] <- info_replicate[i]
  
  seuratObj.TS.multi[i][[1]] <- RenameCells(object = seuratObj.TS.multi[i][[1]], 
                                            new.names = paste0(info_exp[i], "_wk",
                                                               info_age[i], "_",
                                                               info_mouse[i],"_",
                                                               info_replicate[i],"_", 
                                                               colnames(seuratObj.TS.multi[i][[1]])))
}

for (i in 1:length(seuratObj.TS.multi)){
  for (j in 1:ncol(seuratObj.TS.multi[[i]])){
    if (seuratObj.TS.multi[i][[1]]@meta.data[j,"percent.mm10"] > 10)  {
      seuratObj.TS.multi[i][[1]]@meta.data[j,'cr_call'] <- "mm10"
    }
    else{
      seuratObj.TS.multi[i][[1]]@meta.data[j,'cr_call'] <- "GRCh38"
    }
  }
}

table(seuratObj.TS.multi[1][[1]][["cr_call"]])
table(seuratObj.TS.multi[2][[1]][["cr_call"]])
table(seuratObj.TS.multi[3][[1]][["cr_call"]])
table(seuratObj.TS.multi[4][[1]][["cr_call"]])
table(seuratObj.TS.multi[5][[1]][["cr_call"]])
table(seuratObj.TS.multi[6][[1]][["cr_call"]])

### subset

for (i in 1:length(seuratObj.TS.multi)){
  seuratObj.TS.multi[i][[1]] <- subset(seuratObj.TS.multi[i][[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & percent.mm10 < 10)
  Idents(seuratObj.TS.multi[i][[1]]) <- seuratObj.TS.multi[i][[1]]$cr_call
  seuratObj.TS.multi[[i]] <- subset(seuratObj.TS.multi[[i]], ident = "GRCh38")
  Idents(seuratObj.TS.multi[[i]]) <-seuratObj.TS.multi[[i]]$orig.ident
}

### remove mito genes
# Remove noncoding RNAs and mitochondrial genes (RP pseudogenes not functional)
seuratObj.TS.multi.cleanup <- list()
for (i in 1:length(seuratObj.TS.multi)){
  genes.keep <- rownames(seuratObj.TS.multi[i][[1]])[!grepl("^RP11|^MT-", rownames(seuratObj.TS.multi[i][[1]]))]
  seuratObj.TS.multi[i][[1]] <- subset(seuratObj.TS.multi[i][[1]], features = genes.keep)
  
  counts.temp <- GetAssayData(seuratObj.TS.multi[i][[1]], assay = "RNA", slot = "counts")
  seuratObj.TS.multi.cleanup[i][[1]] <- CreateSeuratObject(counts.temp,
                                                           meta.data = seuratObj.TS.multi[i][[1]]@meta.data)
  rm(counts.temp,genes.keep)
}
grep(pattern = "^MT-", x = row.names(x = seuratObj.TS.multi.cleanup[i][[1]]), value = TRUE)
seuratObj.TS.multi <- seuratObj.TS.multi.cleanup
rm(seuratObj.TS.multi.cleanup)

save(seuratObj.TS.multi, file = 'SeuratObj/key_output/TS_6multi_RNA_DFed_Seurat_humanCells_list_20230115.rda')
load('SeuratObj/key_output/TS_6multi_RNA_DFed_Seurat_humanCells_list_20230115.rda')

# -----------------------------------------------------------
# put TS-RNA and TS-Multiome-RNA together (integration within Seurat to construct a larger seuratobj)

obj.list <- c(seuratObj.TS.RNA, seuratObj.TS.multi)
for(i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], selection.method = "vst", 
                                              nfeatures = 4000, verbose = F)
}

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 4000)
# remove RP features during integration
grep(pattern = "^RP[SL][[:digit:]]|^RPSA",
     features,
     value=TRUE, invert=F)
features <- grep(pattern = "^RP[SL][[:digit:]]|^RPSA",features,value=TRUE, invert=T);length(features)
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE)#, vars.to.regress = "percent.rb")
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = obj.list, reference = c(1,2,3,4), reduction = "rpca",
                                  dims = 1:30, anchor.features = 4000)
obj.int.allrna <- IntegrateData(anchorset = anchors, dims = 1:30)

obj.int.allrna <- ScaleData(obj.int.allrna, verbose = TRUE, features = rownames(GetAssayData(obj.int.allrna)))
obj.int.allrna <- RunPCA(obj.int.allrna, verbose = FALSE, npcs = 50)
ElbowPlot(obj.int.allrna, ndims = 50)
obj.int.allrna <- RunUMAP(obj.int.allrna, dims = 1:30, reduction = "pca", metric = "cosine")
DimPlot(obj.int.allrna, group.by = "age") 

obj.int.allrna <- FindNeighbors(obj.int.allrna, k.param = 10, reduction = "pca", dims = 1:30, verbose = F)
obj.int.allrna <- FindClusters(obj.int.allrna, resolution = 0.6, algorithm = 2, print.output = F)


# quickly fix age if needed
obj.int.allrna$age[1:64646] <- sapply(strsplit(obj.int.allrna$age[1:64646], split = "wk"),"[[",2)
table(obj.int.allrna$age)

pdf("TS_all_RNA_lognorm_int_presubset_umap.pdf")
DimPlot(obj.int.allrna, label = T, repel = T) + NoLegend()
dev.off()
pdf("TS_all_RNA_lognorm_int_presubset_vln.pdf", width = 13, height = 10)
VlnPlot(obj.int.allrna, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1, pt.size = -1) + NoLegend()
dev.off()
# low quality cluster: 1

obj.int.allrna.sub <- subset(obj.int.allrna, cells = row.names(obj.int.allrna@meta.data[!obj.int.allrna$integrated_snn_res.0.6 %in% c("1"),]))
obj.int.allrna.sub # 84640
obj.int.allrna.sub <- RunPCA(object = obj.int.allrna.sub, npcs = 50, verbose = F, assay = "integrated")
obj.int.allrna.sub <- RunUMAP(object = obj.int.allrna.sub, reduction = "pca", assay = "integrated",
                       dims = 1:50, metric = "cosine")
obj.int.allrna.sub <- FindNeighbors(obj.int.allrna.sub, k.param = 50, reduction = "pca", assay = "integrated", dims = 1:50)
cluster.res <- 0.5
obj.int.allrna.sub <- FindClusters(obj.int.allrna.sub, resolution = cluster.res, assay = "integrated", algorithm = 2)

pdf("TS_all_RNA_lognorm_int_postsubset1_umap.pdf")
DimPlot(obj.int.allrna.sub, label = T, repel = T) + NoLegend()
dev.off()
pdf("TS_all_RNA_lognorm_int_postsubset1_vln.pdf", width = 13, height = 10)
VlnPlot(obj.int.allrna.sub, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1, pt.size = -1) + NoLegend()
dev.off()
# low quality cluster: 26

obj.int.allrna.sub <- subset(obj.int.allrna.sub, cells = row.names(obj.int.allrna.sub@meta.data[!obj.int.allrna.sub$integrated_snn_res.0.5 %in% c("26"),]))
obj.int.allrna.sub # 84640
obj.int.allrna.sub <- RunPCA(object = obj.int.allrna.sub, npcs = 50, verbose = F, assay = "integrated")
ElbowPlot(obj.int.allrna.sub, ndims = 50)
obj.int.allrna.sub <- RunUMAP(object = obj.int.allrna.sub, reduction = "pca", assay = "integrated",
                              dims = 1:30, metric = "cosine")
obj.int.allrna.sub <- FindNeighbors(obj.int.allrna.sub, k.param = 50, reduction = "pca", assay = "integrated", dims = 1:30)
cluster.res <- 0.5
obj.int.allrna.sub <- FindClusters(obj.int.allrna.sub, resolution = cluster.res, assay = "integrated", algorithm = 2)

pdf("TS_all_RNA_lognorm_int_postsubset2_umap.pdf")
DimPlot(obj.int.allrna.sub, label = T, repel = T) + NoLegend()
dev.off()
pdf("TS_all_RNA_lognorm_int_postsubset2_vln.pdf", width = 13, height = 10)
VlnPlot(obj.int.allrna.sub, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1, pt.size = -1) + NoLegend()
dev.off()
# looking ok for now

save(obj.int.allrna.sub, file = "SeuratObj/key_output/TS_all_RNA_logNorm_int_20230119.rda")
load("SeuratObj/key_output/TS_all_RNA_logNorm_int_20230119.rda")

# -----------------------------------------------------------------
### Annotation, map all RNA to 4H1/JS fetal, by SCANVI

### load 4H1
load('/media/Scratch_SSD_Voyager/sammi/RNA_editing/WT4H1_YWsubset_integratedObj_20230102.Robj')
WT_4H1 <- obj.integrated
rm(obj.integrated)

DefaultAssay(WT_4H1) <- "RNA"
DefaultAssay(obj.int.allrna.sub) <- "RNA"

common_genes <- intersect(row.names(obj.int.allrna.sub), row.names(WT_4H1));length(common_genes)

so_tmp <- CreateSeuratObject(counts = obj.int.allrna.sub@assays$RNA@counts[common_genes,], assay = "RNA", 
                             project = "TS", meta.data = obj.int.allrna.sub@meta.data)
dim(so_tmp)
SaveH5Seurat(so_tmp, filename = paste0("h5ad/", "TS_all_RNA_int_20230119", ".h5Seurat"))
Convert(paste0("h5ad/", "TS_all_RNA_int_20230119", ".h5Seurat"), dest = "h5ad")

WT_4H1_tmp <- CreateSeuratObject(counts = WT_4H1@assays$RNA@counts[common_genes,], assay = "RNA", project = "4H1",
                                 meta.data = WT_4H1@meta.data)
dim(WT_4H1_tmp)
SaveH5Seurat(WT_4H1_tmp, filename = paste0("h5ad/", "4H1_20230119", ".h5Seurat"), overwrite = T)
Convert(paste0("h5ad/", "4H1_20230119", ".h5Seurat"), dest = "h5ad")


# run Scanvi pipeline in command line python



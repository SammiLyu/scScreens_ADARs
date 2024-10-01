### data from Muotri et. al. cortical organoids, time point 1mo, 3mo, 6mo, 10mo
### cellranger using GRCh38 reference
### to save time, using cell type annotations from Muotri et. al.

# ----------------------------------------------------------------------
### couldn't find metadata - already emailed Alysson, but let's try preprocessing

## load cellranger output - sparse count matrix (human gene + mouse gene)
library(Seurat)
library(swne)
library(cowplot)
library(stringr)
library(DoubletFinder)

CO <- list(Read10X('/media/NAS1/Sammi_NAS1/CellRanger_Output/cortical_organoids_cr/1mo/co_1mo/outs/filtered_feature_bc_matrix/'),
           Read10X('/media/NAS1/Sammi_NAS1/CellRanger_Output/cortical_organoids_cr/3mo/co_3mo/outs/filtered_feature_bc_matrix/'),
           Read10X('/media/NAS1/Sammi_NAS1/CellRanger_Output/cortical_organoids_cr/6mo/co_6mo/outs/filtered_feature_bc_matrix/'),
           Read10X('/media/NAS1/Sammi_NAS1/CellRanger_Output/cortical_organoids_cr/10mo/co_10mo/outs/filtered_feature_bc_matrix/'))
names(CO) <- c("co1mo","co3mo","co6mo","co10mo")

## create initial seuratobj
CO_SO.list <- list()
for (i in 1:length(CO)){
  CO_SO.list[i][[1]] <- CreateSeuratObject(counts = CO[[i]], project = "CO_WTiPSCs",
                                               min.cells = round(0.001*ncol(CO[[i]])), min.features = 50)
}

## gene clean up (mito & ribo)
# count mito gene percent & remove mito gene
CO_SO.list.cleanup <- list()
for (i in 1:length(CO_SO.list)){
  CO_SO.list[i][[1]][["percent.mt"]] <- PercentageFeatureSet(CO_SO.list[i][[1]], pattern = "^MT-")
  
  genes.keep <- rownames(CO_SO.list[i][[1]])[!grepl("^MT-", rownames(CO_SO.list[i][[1]]))]
  CO_SO.list.cleanup[i][[1]] <- subset(CO_SO.list[i][[1]], features = genes.keep)
  
  counts.temp <- GetAssayData(CO_SO.list.cleanup[i][[1]], assay = "RNA")
  CO_SO.list.cleanup[i][[1]] <- CreateSeuratObject(counts.temp,
                                                       meta.data = CO_SO.list.cleanup[i][[1]]@meta.data)
  rm(counts.temp,genes.keep)
}
head(CO_SO.list.cleanup[1][[1]]@meta.data)
row.names(CO_SO.list.cleanup[1][[1]])
grep("^MT-", row.names(CO_SO.list.cleanup[1][[1]]))

CO_SO.list <- CO_SO.list.cleanup
rm(CO_SO.list.cleanup)

names(CO_SO.list) <- c("co1mo","co3mo","co6mo","co10mo")

## remove noncoding RNAs - RP11 is NOT in the gene names rn, decide to skip this removal for now unless ribo genes are driving
which(grepl("^RP11", rownames(CO_SO.list[1][[1]])))


# putting in some metadata
info_age <- c("1mo","3mo","6mo","10mo")

for (i in 1:length(CO_SO.list)) {
  CO_SO.list[i][[1]][['age']] <- info_age[i]
  
  CO_SO.list[i][[1]] <- RenameCells(object = CO_SO.list[i][[1]], 
                                        new.names = paste0(info_age[i], "_",
                                                           colnames(CO_SO.list[i][[1]])))
}

head(CO_SO.list[1][[1]]@meta.data)

rm(info_age)

p <- list(list())
for (i in 1:length(CO_SO.list)){
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  p[i][[1]] <- plot_grid(VlnPlot(CO_SO.list[i][[1]], features = c("nFeature_RNA")) + NoLegend(),
                         VlnPlot(CO_SO.list[i][[1]], features = c("nCount_RNA")) + NoLegend(),
                         VlnPlot(CO_SO.list[i][[1]], features = c("percent.mt")) + NoLegend(),
                         ncol = 3, labels = names(CO_SO.list)[i], scale = 0.9)
}
pdf("Figure_CO_R1/CO_4co_Vln_postcr_prefilter_v1.pdf", width = 11, height = 6)
p
dev.off()

### Filtering for features
p <- list(list())
for (i in 1:length(CO_SO.list)){
  CO_SO.list[i][[1]] <- subset(CO_SO.list[i][[1]], subset = nFeature_RNA > 200 & 
                                     nFeature_RNA < 5000 & percent.mt < 10) #& percent.mm10 < 10)
  p[i][[1]] <- plot_grid(VlnPlot(CO_SO.list[i][[1]], features = c("nFeature_RNA")) + NoLegend(),
                         VlnPlot(CO_SO.list[i][[1]], features = c("nCount_RNA")) + NoLegend(),
                         VlnPlot(CO_SO.list[i][[1]], features = c("percent.mt")) + NoLegend(),
                         ncol = 3, labels = names(CO_SO.list)[i], scale = 0.9)
}
pdf("Figure_CO_R1/CO_4co_Vln_postcr_postfilter_v1.pdf", width = 11, height = 6)
p
dev.off()

# preprocessing
for (i in 1:length(CO_SO.list)){
  CO_SO.list[i][[1]] <- NormalizeData(CO_SO.list[i][[1]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  CO_SO.list[i][[1]] <- FindVariableFeatures(CO_SO.list[i][[1]], selection.method = "vst", nfeatures = 4000, verbose = F)
  CO_SO.list[i][[1]] <- ScaleData(CO_SO.list[i][[1]]) # this step could include vars.to.regress
}
hist(log10(CO_SO.list[1][[1]]@assays$RNA@counts@x)) # raw
hist(CO_SO.list[1][[1]]@assays$RNA@data@x) # normalized data
hist(CO_SO.list[1][[1]]@assays$RNA@scale.data) # scaled

p <- list(list())
for (i in 1:length(CO_SO.list)){
  CO_SO.list[i][[1]] <- RunPCA(CO_SO.list[i][[1]], npcs = 30, verbose = F)
  p[i][[1]] <- ElbowPlot(CO_SO.list[i][[1]], ndims = 30)
}
pdf("Figure_CO_R1/CO_4co_Elbow_v1.pdf", width = 11, height = 6)
p
dev.off()

for (i in 1:length(CO_SO.list)){
  CO_SO.list[i][[1]] <- RunUMAP(CO_SO.list[i][[1]], dims = 1:30, verbose = F)
  CO_SO.list[i][[1]] <- FindNeighbors(CO_SO.list[i][[1]], k.param = 10, reduction = "pca", dims = 1:30, verbose = F)
  CO_SO.list[i][[1]] <- FindClusters(CO_SO.list[i][[1]], resolution = 0.6, 
                                         algorithm = 1, print.output = F)
}

save(CO_SO.list, file ='Figure_CO_R1/CO_4co_preprocessed_for_DF_20230607.Rda')
#write.csv(QC_table, "Figure3_ADARCO_R1/CO_4ter_QC_table_20230201.csv", row.names = TRUE, quote = FALSE)

### doublet finder
ListTest <- CO_SO.list

DF_File_Names <- paste0(names(ListTest), "_DFed_SO.rds")
npcs = 75

# model for linear fit between multiplate rates & recovered nuclei)
#doubleRateDF <- data.frame(multiplate_rate = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6)/100,
#                           recovered_cell = c(0.5e3, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3))
#lmDR <- lm(multiplate_rate~recovered_cell, data = doubleRateDF)
#estimated_doublet_rate = sapply(ListTest,dim)[2,]*lmDR$coefficients[2]+lmDR$coefficients[1]
estimated_doublet_rate = 0.076 # loaded 10k cells for all datasets
#QC_table$estimated_doublet_rate <- rep(estimated_doublet_rate, length(ListTest))

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
  
  saveRDS(seuratTemp_DF, file = paste0("Figure_CO_R1/preprocess_20230607/DF/",DF_File_Names[i]))
}

# Load DFed datasets
DF_File_Names <- paste0(c("co1mo","co3mo","co6mo","co10mo"), "_DFed_SO.rds")

CO <- list()
for (i in 1:length(DF_File_Names)){
  CO[[i]] <- readRDS(paste0("Figure_CO_R1/preprocess_20230607/DF/", DF_File_Names[[i]]))
}
names(CO) <- c("co1mo","co3mo","co6mo","co10mo")

rm(ListTest)
rm(sweep.res.list_seurat,sweep.stats_seurat,seuratTemp,seuratTemp_DF,bcmvn)
rm(npcs,pK,nExp_poi,nExp_poi.adj)
gc()

# brute force modifying the colnames of DF results
for (i in 1:length(CO)) {
  colnames(CO[i][[1]]@meta.data)[ncol(CO[i][[1]]@meta.data)-1] <- "DF_pANN"
  colnames(CO[i][[1]]@meta.data)[ncol(CO[i][[1]]@meta.data)] <- "DF_Classification"
  #QC_table$DF_Singlet[i] <- length(which(CO[i][[1]][["DF_Classification"]] == "Singlet"))
}
#QC_table
#write.csv(QC_table, "Figure3_ADARCO_R1/CO_4ter_QC_table_20230201.csv", row.names = TRUE, quote = FALSE)

save(CO, file = 'Figure_CO_R1/preprocess_20230607/CO_4co_DFed_Seurat_list_20230607.rda')
load('Figure_CO_R1/preprocess_20230607/CO_4co_DFed_Seurat_list_20230607.rda')

### Filtering for Singlet Human cells
for (i in 1:length(CO)) {
  Idents(CO[i][[1]]) <- CO[i][[1]]$DF_Classification
  CO[i][[1]] <- subset(CO[i][[1]], ident = "Singlet")
  
  Idents(CO[i][[1]]) <- CO[i][[1]]$orig.ident
}
sum(sapply(CO, dim)[2,]) # 15747

# Remove noncoding RNAs and mitochondrial genes (RP pseudogenes not functional)
CO.cleanup <- list()
for (i in 1:length(CO)){
  genes.keep <- rownames(CO[i][[1]])[!grepl("^RP11|^MT-", rownames(CO[i][[1]]))]
  CO[i][[1]] <- subset(CO[i][[1]], features = genes.keep)
  
  counts.temp <- GetAssayData(CO[i][[1]], assay = "RNA", slot = "counts")
  CO.cleanup[i][[1]] <- CreateSeuratObject(counts.temp,
                                                         meta.data = CO[i][[1]]@meta.data)
  rm(counts.temp,genes.keep)
}
row.names(CO.cleanup[i][[1]])

CO <- CO.cleanup
rm(CO.cleanup)

names(CO) <- c("co1mo","co3mo","co6mo","co10mo")
CO

head(CO[1][[1]]@meta.data)

grep(pattern = "^mm10-", x = row.names(x = CO[i][[1]]), value = TRUE)
grep(pattern = "^GRCh38-", x = row.names(x = CO[i][[1]]), value = TRUE)
grep(pattern = "^MT-", x = row.names(x = CO[i][[1]]), value = TRUE)

save(CO, file = 'Figure_CO_R1/CO_4co_DFed_Seurat_humanCells_list_20230607.rda')
load('Figure_CO_R1/CO_4co_DFed_Seurat_humanCells_list_20230607.rda')


# ----------------------------------------------------------------------
### integration

for(i in 1:length(CO)) {
  CO[[i]] <- NormalizeData(object = CO[[i]], verbose = FALSE)
  CO[[i]] <- FindVariableFeatures(object = CO[[i]], selection.method = "vst", 
                                              nfeatures = 4000, verbose = F)
}

features <- SelectIntegrationFeatures(object.list = CO, nfeatures = 4000)
# remove RP features during integration
grep(pattern = "^RP[SL][[:digit:]]|^RPSA",
     features,
     value=TRUE, invert=F)
features <- grep(pattern = "^RP[SL][[:digit:]]|^RPSA",features,value=TRUE, invert=T);length(features)
CO <- lapply(X = CO, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE)#, vars.to.regress = "percent.rb")
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = CO, reference = c(1,2,3,4), reduction = "rpca",
                                  dims = 1:30, anchor.features = 4000)
obj.int.co <- IntegrateData(anchorset = anchors, dims = 1:30)

obj.int.co <- ScaleData(obj.int.co, verbose = TRUE, features = rownames(GetAssayData(obj.int.co)))
obj.int.co <- RunPCA(obj.int.co, verbose = FALSE, npcs = 50)
ElbowPlot(obj.int.co, ndims = 50)
obj.int.co <- RunUMAP(obj.int.co, dims = 1:30, reduction = "pca", metric = "cosine")
DimPlot(obj.int.co, group.by = "age") 

obj.int.co <- FindNeighbors(obj.int.co, k.param = 10, reduction = "pca", dims = 1:30, verbose = F)
obj.int.co <- FindClusters(obj.int.co, resolution = 0.6, algorithm = 2, print.output = F)

pdf("Figure_CO_R1/CO_4co_all_RNA_lognorm_int_presubset_umap.pdf", width = 13, height = 10)
plot_grid(
DimPlot(obj.int.co, label = T, repel = T) + NoLegend(),
DimPlot(obj.int.co, group.by = "age") 
)
dev.off()
pdf("Figure_CO_R1/CO_4co_all_RNA_lognorm_int_presubset_vln.pdf", width = 13, height = 10)
VlnPlot(obj.int.co, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1, pt.size = -1) + NoLegend()
dev.off()

# drop cluster 6 and 14 and rerun clustering
obj.int.co.sub <- subset(obj.int.co, cells = row.names(obj.int.co@meta.data[!obj.int.co$integrated_snn_res.0.6 %in% c("6","14"),]))
obj.int.co.sub # 14562
obj.int.co.sub <- RunPCA(object = obj.int.co.sub, npcs = 50, verbose = F, assay = "integrated")
obj.int.co.sub <- RunUMAP(object = obj.int.co.sub, reduction = "pca", assay = "integrated",
                       dims = 1:50, metric = "cosine")
obj.int.co.sub <- FindNeighbors(obj.int.co.sub, k.param = 50, reduction = "pca", assay = "integrated", dims = 1:50)
cluster.res <- 0.5
obj.int.co.sub <- FindClusters(obj.int.co.sub, resolution = cluster.res, assay = "integrated", algorithm = 2)

pdf("Figure_CO_R1/CO_4co_all_RNA_lognorm_int_postsubset1_umap.pdf", width = 13, height = 10)
plot_grid(
DimPlot(obj.int.co.sub, label = T, repel = T) + NoLegend(),
DimPlot(obj.int.co.sub, group.by = "age") 
)
dev.off()
pdf("Figure_CO_R1/CO_4co_all_RNA_lognorm_int_postsubset1_vln.pdf", width = 13, height = 10)
VlnPlot(obj.int.co.sub, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 1, pt.size = -1) + NoLegend()
dev.off()

save(obj.int.co.sub, file = "Figure_CO_R1/CO_4co_RNA_logNorm_int_20230607.rda")
load("Figure_CO_R1/CO_4co_RNA_logNorm_int_20230607.rda")

# ----------------------------------------------------------------------
### cell type annotation
### call cell types based on the manuscript
### cell types: GABAergic, Glutamatergic, Glia, IP, Progenitor, MC, Other
### markers: GABAergic - GAD2, DLX1, DLX5
### markers: Glutamatergic - NEUROD2, MAPT, SNAP25, GRIA2
### markers: Glia - TTYH1, MT2A, SLC1A3
### markers: IP - GADD45G, EOMES, TAC3, NHLH1
### markers: Progenitor - ID3, ID1, PAX6, CA2, OTX2
### markers: Microglia - AIF1, TMEM119

DefaultAssay(obj.int.co.sub) <- "RNA"
manu_markers <- c("GAD2", "DLX1", "DLX5",
                  "NEUROD2","MAPT","SNAP25","GRIA2",
                  "TTYH1","MT2A","SLC1A3",
                  "GADD45G","EOMES","TAC3","NHLH1",
                  "ID3","ID1","PAX6","CA2","OTX2",
                  "AIF1","TMEM119","SALL1","PTPRC","CD68","ITGAM")
pdf("Figure_CO_R1/CO_4co_all_RNA_dotplot_manuscriptmarkers.pdf", width = 10, height = 5)
DotPlot(obj.int.co.sub, features = manu_markers) + RotatedAxis()
dev.off()

pdf("Figure_CO_R1/CO_4co_all_RNA_feature_manuscriptmarkers_gabaergic.pdf", width = 10, height = 8)
FeaturePlot(obj.int.co.sub, features = c("GAD2","DLX1","DLX5"))
dev.off()

obj.int.co.sub <- FindSubCluster(obj.int.co.sub, cluster = c(4), graph.name = "integrated_snn", subcluster.name = "sub4", resolution = 0.8, algorithm = 1)
pdf("Figure_CO_R1/CO_4co_all_RNA_feature_manuscriptmarkers_gabaergic_2.pdf", width = 20, height = 8)
plot_grid(
DimPlot(obj.int.co.sub, label = T, repel = T, group.by = "sub4") + NoLegend(),
DimPlot(obj.int.co.sub, cells.highlight = row.names(obj.int.co.sub@meta.data[obj.int.co.sub$sub4 == c("4_3"),]),label = T, repel = T, group.by = "sub4") + NoLegend(),
FeaturePlot(obj.int.co.sub, features = c("GAD2","DLX1","DLX5")),
DotPlot(obj.int.co.sub, features = c("GAD2","DLX1","DLX5"), group.by = "sub4") + RotatedAxis(), ncol = 4, rel_widths = c(3,3,3,2)
)
dev.off()

cluster.mapping <- c("0"="Glutamatergic",
                     "1"="Glia", 
                     "2"="Progenitor", 
                     "3"="Glutamatergic",
                     "4_0"=NA,
                     "4_1"=NA,
                     "4_2"=NA,
                     "4_3"="GABAergic",
                     "4_4"=NA,
                     "5"="Progenitor", 
                     "6"="IP",
                     "7"="Glutamatergic",
                     "8"="Glia",
                     "9"="Glutamatergic",
                     "10"=NA,
                     "11"="Progenitor",
                     "12"=NA, 
                     "13"=NA, 
                     "14"=NA,
                     "15"="Progenitor",
                     "16"=NA)

obj.int.co.sub$cluster_l1 <- plyr::revalue(obj.int.co.sub$sub4, replace = cluster.mapping)

Idents(obj.int.co.sub) <- obj.int.co.sub$cluster_l1
pdf("Figure_CO_R1/CO_4co_all_RNA_annotated_dotplot_manuscriptmarkers_20230627.pdf", width = 7, height = 5)
DotPlot(obj.int.co.sub[,!is.na(obj.int.co.sub$cluster_l1)], features = manu_markers) + RotatedAxis()
dev.off()

DotPlot(obj.int.co.sub[,!is.na(obj.int.co.sub$cluster_l1)], features = c("ADAR","ADARB1","ADARB2")) + RotatedAxis()
FeaturePlot(obj.int.co.sub[,!is.na(obj.int.co.sub$cluster_l1)], features = c("ADAR","ADARB1","ADARB2")) + RotatedAxis()
DimPlot(obj.int.co.sub, label = T, repel = T, group.by = "cluster_l1") + NoLegend()


# get levels correct
obj.int.co.sub$age <- factor(obj.int.co.sub$age, levels = c('1mo','3mo','6mo','10mo'))
obj.int.co.sub$cluster_l1 <- factor(obj.int.co.sub$cluster_l1, levels = c("Progenitor","IP","Glia","Glutamatergic","GABAergic"))

save(obj.int.co.sub, file = "Figure_CO_R1/CO_4co_RNA_logNorm_int_annot_20230627.rda")
load("Figure_CO_R1/CO_4co_RNA_logNorm_int_annot_20230627.rda")

cluster_levels <- c("Progenitor","IP","Glia","Glutamatergic","GABAergic")
cluster_col <- c("#442288","#6CA2EA","#B5D33D","#A62639","#EB7D5B")

names(cluster_col) <- cluster_levels

pdf('Figure_CO_R1/CO_4co_intUMAP_20230717.pdf', height = 10, width = 10)
DimPlot(obj.int.co.sub[,!is.na(obj.int.co.sub$cluster_l1)], group.by = "cluster_l1", label = T, repel = T, cols = cluster_col) + NoLegend()
dev.off()

DefaultAssay(obj.int.co.sub) <- "RNA"
DotPlot(obj.int.co.sub[,!is.na(obj.int.co.sub$cluster_l1)], features = c("EOMES","GADD45G","NHLH1","SLC17A7","SLC17A6","GAD1","GAD2","RBFOX3","NEUROD2"), group.by = "cluster_l1") + RotatedAxis()


# plot mo10 only
pdf('Figure_CO_R1/CO_10mo_intUMAP_20230723.pdf', height = 10, width = 10)
Idents(obj.int.co.sub) <- obj.int.co.sub$age
obj.int.co.sub.10mo <- subset(obj.int.co.sub, ident = "10mo")
DimPlot(obj.int.co.sub.10mo[,!is.na(obj.int.co.sub.10mo$cluster_l1)], group.by = "cluster_l1", label = T, repel = T, cols = cluster_col) + NoLegend()
dev.off()


table(obj.int.co.sub$age, obj.int.co.sub$cluster_l1)
write.table(table(obj.int.co.sub$age, obj.int.co.sub$cluster_l1), file = "Figure_CO_R1/CO_4co_celltype_count.txt", sep = ',', row.names = T, col.names = T, quote = F)

# ----------------------------------------------------------------------
### save cell barcodes into separate bam files for parsing bams

for (i in names(table(obj.int.co.sub$cluster_l1))){
  for (j in names(table(obj.int.co.sub$age))){
    if (length (row.names(obj.int.co.sub@meta.data[obj.int.co.sub$cluster_l1 == i & obj.int.co.sub$age == j, ])) != 0){
      write.table(paste0("CB:Z:",sapply(strsplit(row.names(obj.int.co.sub@meta.data[obj.int.co.sub$cluster_l1 == i & !is.na(obj.int.co.sub$cluster_l1) & obj.int.co.sub$age == j, ]), split = "_"), "[[",2)), 
                  paste0('cb_listsCO/v2/',i,'_',j,'.txt'), col.names = FALSE, row.names = FALSE, 
                  quote = FALSE, sep = ',')
    }
  }
}

neuro_ct_markers <- c("HMGB2","PCLAF","TOP2A","CENPF","MKI67",
                      "ATP1A2","ANGPT1",
                      "ID2","SOX2","NES","FABP7",
                      "NCAM1","MBP","BCAN","RFX4",
                      "FOXG1","LHX2",
                      "EOMES","NEUROG2",
                      "SLC17A7","NEUROD6","FEZF2",
                      "SLC17A6","KCNJ6","LMX1B",
                      "POU3F2","ASCL1","PRMT8",
                      "ASIC2","TAFA1","FSTL5","GRID1",
                      "DLX6-AS1","GAD2",
                      "SLC6A5","RBFOX1","GAD1","ADARB2","LAMP5",
                      "TFEC","YAP1",
                      "OTX2","MITF","PMEL","BEST1","RLBP1",
                      "SIX6","VSX2","RAX","HES1",
                      "TRPM1",
                      "PRDM13","ONECUT1",
                      "GAP43","POU4F2","PRPH",
                      "RSPO2","SOX2-OT","ZIC5",
                      "LMX1A","FOXJ1","PIFO"
                      )
pdf('Figure_CO_R1/CO_4co_all_RNA_annotated_dotplot_neuroctmarkers_20230627.pdf', width = 15, height = 5)
DotPlot(obj.int.co.sub[,!is.na(obj.int.co.sub$cluster_l1)], features = neuro_ct_markers) + RotatedAxis()
dev.off()


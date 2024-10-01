### script pre-process H1 Teratoma Time Series snRNAseq data
library(Seurat)
library(ggplot2)
library(cowplot)
library(DoubletFinder)
library(SeuratDisk)
library(cellMapper)
library(swne)
library(RColorBrewer)
library(stringr)


# read Scanvi output
TS_allRNA_scanvi.meta <- read.table('anndata_anno_20230121.csv', sep = ",", header = 1, row.names = 1)
head(TS_allRNA_scanvi.meta)

# quick fix on DF_classification
DF_Classification <- c(TS_allRNA_scanvi.meta$DF_Classifications[1:57431],TS_allRNA_scanvi.meta$DF_Classifications[57432:84079])
length(DF_Classification)
dim(TS_allRNA_scanvi.meta)
TS_allRNA_scanvi.meta$DF_Classifications <- DF_Classification
TS_allRNA_scanvi.meta$DF_Classification <- NULL
names(TS_allRNA_scanvi.meta)[13] <- "DF_Classification"

head(TS_allRNA_scanvi.meta)

# put back the annotation with some cleanup
load("SeuratObj/key_output/TS_all_RNA_logNorm_int_20230119.rda")
obj.int.allrna.sub

head(obj.int.allrna.sub@meta.data)
obj.int.allrna.sub@meta.data <- TS_allRNA_scanvi.meta

DimPlot(obj.int.allrna.sub, group.by = "scANVI_predictions", label = T, repel = T) + NoLegend()

obj.int.allrna.sub <- RunUMAP(object = obj.int.allrna.sub, reduction = "pca", assay = "integrated",
                              dims = 1:30, metric = "cosine")
obj.int.allrna.sub <- FindNeighbors(obj.int.allrna.sub, k.param = 50, reduction = "pca", assay = "integrated", dims = 1:30)
cluster.res <- 0.5
obj.int.allrna.sub <- FindClusters(obj.int.allrna.sub, resolution = cluster.res, assay = "integrated", algorithm = 2)

plot_grid(DimPlot(obj.int.allrna.sub, group.by = "scANVI_predictions", label = T, repel = T) + NoLegend(),
          DimPlot(obj.int.allrna.sub, group.by = "integrated_snn_res.0.5", label = T, repel = T) + NoLegend()
)

DimPlot(obj.int.allrna.sub, group.by = "age", label = T, repel = T)

# for each cell type, plot a separate cell highlight plot
p <- list()
for (i in names(table(obj.int.allrna.sub$scANVI_predictions))){
  p[[i]] <- DimPlot(obj.int.allrna.sub, cells.highlight = row.names(obj.int.allrna.sub@meta.data[obj.int.allrna.sub$scANVI_predictions == i,]),
                    group.by = "scANVI_predictions", label = T, repel = T) + NoLegend() + ggtitle(i)

}
pdf("TS_all_RNA_scanvi_cell_type_highlight.pdf")
p
dev.off()


## confirming some specific cell types
DefaultAssay(obj.int.allrna.sub) <- "RNA"

# DimPlot(obj.int.allrna.sub, cells.highlight = row.names(obj.int.allrna.sub@meta.data[obj.int.allrna.sub$scANVI_predictions == "CycProg",]),group.by = "scANVI_predictions", label = T, repel = T) + NoLegend()

# DotPlot(obj.int.allrna.sub,features = YanMarkersPub) + RotatedAxis()

# DotPlot(obj.int.allrna.sub,features = c("SOX17","FOXA2")) + RotatedAxis() # Definitive Endoderm, Ranking programming factors

# DotPlot(obj.int.allrna.sub,features = c("RBPJL","GP2")) + RotatedAxis() # Pancreatic Prog
# # https://pubmed.ncbi.nlm.nih.gov/34638511/#:~:text=RBPJL%20is%20specifically%20expressed%20in,in%20pancreatic%20tumour%20cell%20lines.
# # https://www.nature.com/articles/s41467-017-00561-0#:~:text=Glycoprotein%202%20is%20a%20specific%20cell%20surface%20marker%20of%20human%20pancreatic%20progenitors,-Kathryn%20F.

# markers12 <- FindMarkers(obj.int.allrna.sub, ident.1 = 12, only.pos = T, logfc.threshold = 0.4, min.pct = 0.4)
# markers12 <- markers12[order(markers12$avg_log2FC, decreasing = T),]
# head(markers12,50)

# DotPlot(obj.int.allrna.sub,features = c("RORB","CLVS1","CRB1")) + RotatedAxis() # 

# markers17 <- FindMarkers(obj.int.allrna.sub, ident.1 = 17, only.pos = T, logfc.threshold = 0.4, min.pct = 0.4)
# markers17 <- markers17[order(markers17$avg_log2FC, decreasing = T),]
# head(markers17,50)

# DotPlot(obj.int.allrna.sub,features = c("FMN1","PKHD1","PAX2","PAX8","EYA1")) + RotatedAxis() # also Kidney, but different marker gene
# https://journals.lww.com/jasn/Abstract/2021/08000/Uncovering_Modifier_Genes_of_X_Linked_Alport.18.aspx
# https://pubmed.ncbi.nlm.nih.gov/15108277/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3876955/

# FeaturePlot(obj.int.allrna.sub, features = c("FOXJ1","CDHR3"))
obj.int.allrna.sub <- FindSubCluster(obj.int.allrna.sub, cluster = "19", graph.name = "integrated_snn", subcluster.name = "sub.cluster.19", algorithm = 2)
# DimPlot(obj.int.allrna.sub, group.by = "sub.cluster.19", label = T, repel = T)
# DotPlot(obj.int.allrna.sub, features = c("FOXJ1","CDHR3"), group.by = "sub.cluster.19")
# DimPlot(obj.int.allrna.sub, cells.highlight = row.names(obj.int.allrna.sub@meta.data[obj.int.allrna.sub$sub.cluster.19 == "19_4",]),group.by = "sub.cluster.19", label = T, repel = T) + NoLegend()

# Amir wonders if we see liver prog
# FeaturePlot(obj.int.allrna.sub, features = c("HNF4A","NR5A2")) # hapatic progenitor
# DotPlot(obj.int.allrna.sub, features = c("HNF4A","NR5A2","PTPRC","LMO2"), group.by = "sub.cluster.19") # CD45 = PTPRC
#https://www.nature.com/articles/s41598-017-15973-7

# quick clean-up
obj.int.allrna.sub$annot_v1 <- NA
# keep chondrogenic MSC/Fib Annotations from scANVI
obj.int.allrna.sub@meta.data[obj.int.allrna.sub$scANVI_predictions %in% c("Chondrogenic MSC/Fib"), "annot_v1"] = "Chondrogenic MSC/Fib"
obj.int.allrna.sub@meta.data[obj.int.allrna.sub$integrated_snn_res.0.5 %in% c("6","8","18") & is.na(obj.int.allrna.sub$annot_v1), "annot_v1"] = 
  obj.int.allrna.sub@meta.data[obj.int.allrna.sub$integrated_snn_res.0.5 %in% c("6","8","18") & is.na(obj.int.allrna.sub$annot_v1), "scANVI_predictions"]

obj.int.allrna.sub@meta.data[obj.int.allrna.sub$sub.cluster.19 == "19_4" , "annot_v1"] = "Airway Epi"
obj.int.allrna.sub@meta.data[obj.int.allrna.sub$integrated_snn_res.0.5 %in% c("19") & obj.int.allrna.sub$sub.cluster.19 != "19_4" , "annot_v1"] = "Retinal Epi"

cluster.mapping <- c("0"="Adipogenic MSC/Fib",
                     "1"="Early Neurons", 
                     "2"="Radial Glia", 
                     "3"="MSC/Fib",
                     "4"="Foregut Epi",
                     "5"="Mid/Hindgut Epi", 
                     "7"="Radial Glia",
                     "9"="MyoFib",
                     "10"="Cycling MSC/Fib",
                     "11"="Retinal Neurons",
                     "12"="Radial Glia", 
                     "13"="Schwann Cells", 
                     "14"="Early Neurons",
                     "15"="Pericytes",
                     "16"="CycProg", 
                     "17"="Kidney Prog", 
                     "20"="Smooth Muscle", 
                     "21"="Pancreatic Prog", 
                     "22"="Early Neurons",
                     "23"="Definitive Endoderm", 
                     "24"="Immune", 
                     "25"="Mid/Hindgut Epi", 
                     "26"="Early Neurons",
                     "27"="Kidney Prog")
for (i in 1:ncol(obj.int.allrna.sub)){
  if (obj.int.allrna.sub@meta.data[i,"integrated_snn_res.0.5"] %in% names(cluster.mapping)){
    if (is.na(obj.int.allrna.sub@meta.data[i,"annot_v1"])){
      obj.int.allrna.sub@meta.data[i,"annot_v1"] = cluster.mapping[which(names(cluster.mapping) == obj.int.allrna.sub@meta.data[i,"integrated_snn_res.0.5"])][[1]]
    }
  }
}

### check plot again
plot_grid(DimPlot(obj.int.allrna.sub, group.by = "scANVI_predictions", label = T, repel = T) + NoLegend(),
          DimPlot(obj.int.allrna.sub, group.by = "annot_v1", label = T, repel = T, shuffle = T) + NoLegend()
)
dim(obj.int.allrna.sub@meta.data[is.na(obj.int.allrna.sub$annot_v1),]) # sanity check

### make sure some elements in the seurat object metadata is factored
obj.int.allrna.sub$age <- factor(obj.int.allrna.sub$age, levels = c("4","6","8","10"))
obj.int.allrna.sub$mouse <- factor(obj.int.allrna.sub$mouse, levels = c("M1","M14","M15","M9","M8","M10","A","B"))
obj.int.allrna.sub$annot_v1 <- factor(obj.int.allrna.sub$annot_v1, levels = c("Radial Glia","CycProg","Early Neurons","Retinal Neurons","Retinal Epi","Schwann Cells",
                                                                              "Definitive Endoderm","Foregut Epi","Airway Epi","Pancreatic Prog","Mid/Hindgut Epi","Immune","Adipogenic MSC/Fib",
                                                                              "Chondrogenic MSC/Fib","MSC/Fib","Cycling MSC/Fib","MyoFib",
                                                                              "Muscle Prog","Cardiac/Skeletal Muscle","Pericytes","Smooth Muscle","Kidney Prog"))


# save seurat object & metadata
save(obj.int.allrna.sub, file = "SeuratObj/key_output/TS_all_RNA_logNorm_int_annot_v1_20230124.rda")
write.csv(obj.int.allrna.sub@meta.data, file = "SeuratObj/key_output/TS_all_RNA_logNorm_int_annot_v1_metadata_20230124.csv", row.names = TRUE, quote = FALSE)
read.csv("SeuratObj/key_output/TS_all_RNA_logNorm_int_annot_v1_metadata_20230124.csv", header = T, row.names = 1, sep = ",")

load("/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/SeuratObj/key_output/TS_all_RNA_logNorm_int_annot_v1_20230124.rda")
## save some plots

# pre-annotation plot
pdf("TS_all_RNA_scanvi_cell_type_annot_20230124.pdf", height = 10, width = 15)
plot_grid(DimPlot(obj.int.allrna.sub, group.by = "scANVI_predictions", label = T, repel = T) + NoLegend(),
          DimPlot(obj.int.allrna.sub, group.by = "integrated_snn_res.0.5", label = T, repel = T) + NoLegend()
)
dev.off()

# scanvi_hist
pdf("TS_all_RNA_scanvi_cell_type_scanvi_hist_20230124.pdf")
hist(obj.int.allrna.sub$scANVI_prediction_scores)
dev.off()

# dot plot of marker gene exp
pdf("TS_all_RNA_scanvi_cell_type_annot_dotplot_20230124.pdf", height = 10, width = 15)
DotPlot(obj.int.allrna.sub, features = c(YanMarkersPub,"SOX17","FOXA2",
                                        "RBPJL","GP2",
                                        "RORB","CLVS1","CRB1",
                                        "FMN1","PKHD1","PAX2","PAX8","EYA1")) + RotatedAxis()
dev.off()

# heatmap of marker gene exp
YanMarkersSub <- intersect(c("VIM","SOX2","HES5","HMGB2","RORB","CLVS1","CRB1","DCX","MAP2",
                             "OTX2","NRL","VSX2","ASCL1","ATOH7","NEUROD1","MITF","TTR",
                             "MPZ","SOX10","MLANA",
                             "SOX17","FOXA2","ELF3","PAX9","KRT4","FOXJ1","CDHR3","CDX2","RBPJL","GP2",
                             "CD74","CD34","HHEX","THY1","ITM2A","SHOX2","COL2A1","SOX9","COL15A1",
                             "MYOD1","PAX7","TNNT1","TNNI2","FOXC1","CYP1B1","ACTA2","RGS5","WT1","FMN1","PKHD1","PAX2","PAX8","EYA1"
                             ), VariableFeatures(obj.int.allrna.sub, assay = "integrated"))
heat.mat <- t(apply(GetAssayData(obj.int.allrna.sub, assay = "integrated", slot = "scale.data")[YanMarkersSub,], 1, function(x) {
  tapply(x, obj.int.allrna.sub$annot_v1, mean)
}))
heat.mat[heat.mat > 5] <- 5
ggHeat(heat.mat, clustering = "none", x.lab.size = 11, y.lab.size = 11) + ggtitle('Marker Gene Expression for Teratoma Cell Types')


YanMarkersSub <- c("VIM","SOX2","HES5","HMGB2","RORB","CLVS1","CRB1","DCX","MAP2",
                   "OTX2","NRL","VSX2","ASCL1","ATOH7","NEUROD1","MITF","TTR",
                   "MPZ","SOX10","MLANA",
                             "SOX17","FOXA2","ELF3","PAX9","KRT4","FOXJ1","CDHR3","CDX2","RBPJL","GP2",
                             "CD74","CD34","HHEX","COL14A1","THY1","PDGFRA","ITM2A","SHOX2","COL2A1","SOX9","COL15A1",
                             "MYOD1","PAX7","TNNT1","TNNI2","FOXC1","CYP1B1","ACTA2","RGS5","WT1","FMN1","PKHD1","PAX2","PAX8","EYA1"
)
#obj.int.allrna.sub <- ScaleData(obj.int.allrna.sub)
heat.mat <- t(apply(GetAssayData(obj.int.allrna.sub, assay = "RNA", slot = "scale.data")[YanMarkersSub,], 1, function(x) {
  tapply(x, obj.int.allrna.sub$annot_v1, mean)
}))
heat.mat[heat.mat > 5] <- 5
ggHeat(heat.mat, clustering = "none", x.lab.size = 11, y.lab.size = 11) + ggtitle('Marker Gene Expression for Teratoma Cell Types')


pdf("TS_allRNA_MarkerExpression_20230124.pdf", width = 6, height = 8)
ggHeat(heat.mat, clustering = "none", x.lab.size = 11, y.lab.size = 11) + ggtitle('Marker Gene Expression for Teratoma Cell Types')
dev.off()

pdf("TS_all_RNA_scanvi_cell_type_annotated_v1_20230124.pdf")
DimPlot(obj.int.allrna.sub, group.by = "annot_v1", label = T, repel = T, shuffle = T, cols = cluster_col) + NoLegend()
dev.off()



cluster_col <- c("#4795F5","#0B62CD","#7474FF","#DADAFF","#F2B77C","#9EEA00","#f07cab",
                 "#002700","#00A700","#d30b94","#005A00","#04409B",
                 "#0000FF","#FF0000","#AD07E3","#00C000","#000000",
                 "#B35A00","#FF9933","#A00000","#94641F","#610051")
names(cluster_col) <- levels(obj.int.allrna.sub$annot_v1)

pdf("/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/TS_all_RNA_scanvi_cell_type_annotated_v1_20230210.pdf", height = 8, width = 12)
DimPlot(obj.int.allrna.sub, group.by = "annot_v1", label = F, repel = T, shuffle = T, cols = cluster_col)
dev.off()

obj.int.allrna.sub_ds <- SplitObject(obj.int.allrna.sub, split.by = "age")
sapply(obj.int.allrna.sub_ds, dim)
for (i in 1:length(obj.int.allrna.sub_ds)){
  obj.int.allrna.sub_ds[[i]] <- obj.int.allrna.sub_ds[[i]][, sample(colnames(obj.int.allrna.sub_ds[[i]]), size = min(sapply(obj.int.allrna.sub_ds, dim)[2,]), replace = F )]
}
obj.int.allrna.sub_ds <- merge(obj.int.allrna.sub_ds[[1]], y = c(obj.int.allrna.sub_ds[[2]], obj.int.allrna.sub_ds[[3]], obj.int.allrna.sub_ds[[4]]), 
                               merge.data = T, merge.dr = "umap")
table(obj.int.allrna.sub_ds$age)

# downsample to keep even # of cells in plot for age
pdf("TS_all_RNA_scanvi_cell_type_annotated_v1_groupbyage_ds_20230130.pdf", height = 12, width = 15)
DimPlot(obj.int.allrna.sub_ds, group.by = "annot_v1", label = F, repel = T, shuffle = T, split.by = "age", ncol = 2, cols = cluster_col)
dev.off()

obj.int.allrna.sub_ds$annot_v1 <- factor(obj.int.allrna.sub_ds$annot_v1, levels = levels(obj.int.allrna.sub$annot_v1))
obj.int.allrna.sub_ds$age <- factor(obj.int.allrna.sub_ds$age, levels = c("4","6","8","10"))

pdf("TS_all_RNA_scanvi_cell_type_annotated_v1_groupbyage_ds_20230201.pdf", height = 12, width = 15)
DimPlot(obj.int.allrna.sub_ds, group.by = "annot_v1", label = F, repel = T, shuffle = T, split.by = "age", ncol = 2, cols = cluster_col)
dev.off()



# gene exp validation by Jay Shendure's fetal cell atlas
load("/media/Home_Raid1_Voyager/sammi/Teratoma_Analysis/TimeSeries/JSHF_SO_down100k.Robj")
JSHF_ct_collapse <- read.table("/media/Scratch_SSD_Voyager/sammi/20220406_Teratoma_TS/JSHF_Collapsing.csv", header = TRUE, sep = ",")

levels(factor(JSHF_SO_down$Main_cluster_name))
list_cluster <- factor(JSHF_SO_down$Main_cluster_name)
cl.replace <- JSHF_ct_collapse$New.Ident
names(cl.replace) <- JSHF_ct_collapse$Orig.Ident
list_cluster <- plyr::revalue(list_cluster, replace = cl.replace)
JSHF_SO_down$Collapsed_CT <- list_cluster
table(JSHF_SO_down$Collapsed_CT)
JSHF_ST_co_order <- c("Astrocytes","Neurons","Ex.Neurons","In.Neurons","In.Interneurons",
                      "Retinal.Prog","Retinal.Neurons","Retinal.Interneurons",
                      "Ganglion.Cells","Horizontal.Cells","Retinal.Pig",
                      "Retinal.Epi","Retinal.Fiber","Oligodendrocytes",
                      "Neuroendocrine.Cells",
                      "Schwann.Cells","Glia",
                      
                      "Foregut.Epi","Lung.Epi","Gut.Epi","Heart.Epi","Heart.Cells",
                      "Lymphatic.Endo","Lymphatic.Epi","Endometrial.Epi","Vascular.Endo",
                      
                      "Microglia","Immune.Cells","Blood.Cells",
                      "MSC","Stromal.Cells","Skeletal.Muscle.Cells", "Cardiac.Muscle.Cells",
                      "Mesothelial.Cells","Smooth.Muscle.Cells",
                      "Kidney",
                      
                      "Squamous.Epithelial.Cells","Hepatoblasts","Sympathoblasts","Trophoblasts")
JSHF_SO_down$Collapsed_CT <- factor(JSHF_SO_down$Collapsed_CT, levels = JSHF_ST_co_order)
length(levels(factor(JSHF_SO_down$Collapsed_CT)))

hist(JSHF_SO_down@assays$RNA@counts@x)
JSHF_SO_down <- NormalizeData(JSHF_SO_down)
JSHF_SO_down <- FindVariableFeatures(JSHF_SO_down, nfeatures = 3000)
JSHF_SO_down <- ScaleData(JSHF_SO_down, features = VariableFeatures(JSHF_SO_down))
genes.use <- unique(c(as.character(YanMarkersPub), VariableFeatures(JSHF_SO_down), VariableFeatures(obj.int.allrna.sub)))
genes.use <- Reduce(intersect, list(genes.use, rownames(GetAssayData(JSHF_SO_down, assay = "RNA", slot = "scale.data")), 
                                    rownames(GetAssayData(obj.int.allrna.sub, assay = "integrated"))))
length(genes.use) # 970
Tera.r.assay <- GetAssayData(obj.int.allrna.sub, assay = "integrated")

cor.result <- MapClustersCor(query.data = Tera.r.assay, query.clusters = obj.int.allrna.sub$annot_v1, 
                             train.data = GetAssayData(JSHF_SO_down, assay = "RNA", slot = "scale.data"),
                             train.clusters = JSHF_SO_down$Collapsed_CT,
                             genes.use = genes.use)
ggHeat(t(cor.result$cluster.correlations[complete.cases(cor.result$cluster.correlations),]), clustering = "both", x.lab.size = 11, y.lab.size = 11)
# when plotting, remove neuroendocrine cells, blood cells, heart epi, heart cells, lymphatic endo, lymphatic epi, endometrial epi
# mesothelial cells, squamous epithelial cells, hepatoblasts, sympathoblasts, trophoblasts
ggHeat(t(cor.result$cluster.correlations[complete.cases(cor.result$cluster.correlations),][JSHF_ST_co_order[-c(15,21,22,23,24,25,29,34,37,38,39,40)],]), x.lab.size = 11, y.lab.size = 11)

pdf("TS_all_RNA_correlateJSHF_20230124.pdf", height = 7, width = 10)
ggHeat(t(cor.result$cluster.correlations[complete.cases(cor.result$cluster.correlations),][JSHF_ST_co_order[-c(15,21,22,23,24,25,29,34,37,38,39,40)],]), x.lab.size = 11, y.lab.size = 11) + 
  labs(title='Correlation between Cell Types from Teratoma and Fetal Cell Atlas', x='Fetal Cell Atlas Cell Types', y='Teratoma Time Series Cell Types')
dev.off()

## all Jay Shendure's cell types
# cor.result <- MapClustersCor(query.data = Tera.r.assay, query.clusters = obj.int.allrna.sub$annot_v1, 
#                              train.data = GetAssayData(JSHF_SO_down, assay = "RNA", slot = "scale.data"),
#                              train.clusters = JSHF_SO_down$Main_cluster_name,
#                              genes.use = genes.use)
# ggHeat(t(cor.result$cluster.correlations[complete.cases(cor.result$cluster.correlations),]), clustering = "both", x.lab.size = 11, y.lab.size = 11)

rm(JSHF_SO_down)



### some supplementary plots

# stacked bar plots of cell type change over development time
obj.int.allrna.sub$age_mouse_replicate <- paste0("wk",obj.int.allrna.sub$age,"_",obj.int.allrna.sub$mouse,"_",obj.int.allrna.sub$replicate)
cc_per <- table(obj.int.allrna.sub$annot_v1, obj.int.allrna.sub$age_mouse_replicate)
cc_per <- apply(cc_per, 2, function(x) x/sum(x))
colSums(cc_per)

cc_per_order <- colnames(cc_per)[c(6:14,1:5)]
cc_per <- cc_per[,cc_per_order]

df <- data.frame(celltype=rep(row.names(cc_per),ncol(cc_per)),
                 age=rep(colnames(cc_per),each = nrow(cc_per)),
                 per=unlist(data.frame(cc_per)))
df
#n <- nrow(cc_per)
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

cluster_col <- c("#4795F5","#0B62CD","#7474FF","#DADAFF","#F2B77C","#9EEA00","#f07cab",
                 "#002700","#00A700","#d30b94","#005A00","#04409B",
                 "#0000FF","#FF0000","#AD07E3","#00C000","#000000",
                 "#B35A00","#FF9933","#A00000","#94641F","#610051")
names(cluster_col) <- c("Radial Glia","CycProg","Early Neurons","Retinal Neurons","Retinal Epi","Schwann Cells",
                        "Definitive Endoderm","Foregut Epi","Airway Epi","Pancreatic Prog","Mid/Hindgut Epi","Immune","Adipogenic MSC/Fib",
                        "Chondrogenic MSC/Fib","MSC/Fib","Cycling MSC/Fib","MyoFib",
                        "Muscle Prog","Cardiac/Skeletal Muscle","Pericytes","Smooth Muscle","Kidney Prog")

pdf("TS_all_RNA_stackedbar_2023201.pdf", height = 8)
ggplot(df, aes(x = factor(age, levels = cc_per_order) , y = per, 
               fill = factor(celltype, levels = c("Radial Glia","CycProg","Early Neurons","Retinal Neurons","Retinal Epi","Schwann Cells",
                                                  "Definitive Endoderm","Foregut Epi","Airway Epi","Pancreatic Prog","Mid/Hindgut Epi","Immune","Adipogenic MSC/Fib",
                                                  "Chondrogenic MSC/Fib","MSC/Fib","Cycling MSC/Fib","MyoFib",
                                                  "Muscle Prog","Cardiac/Skeletal Muscle","Pericytes","Smooth Muscle","Kidney Prog")))) +
  geom_col() + scale_fill_manual(values = cluster_col) +
  xlab("age") + ylab("percentage of cells") + theme_classic() + RotatedAxis() + guides(fill=guide_legend(title="Cell Types"))
dev.off()

# instead of stacked bar plots, let's try a line plot
head(obj.int.allrna.sub)
cc_per <- DistMatPlot(obj.int.allrna.sub, "annot_v1", "age_mouse_replicate")
cc_per <- apply(cc_per, 2, function(x) x/sum(x))
colSums(cc_per)

cc_per_order <- colnames(cc_per)[c(6:14,1:5)]
cc_per <- cc_per[,cc_per_order]

df <- data.frame(celltype=rep(row.names(cc_per),ncol(cc_per)),
                 exp=rep(colnames(cc_per),each = nrow(cc_per)),
                 per=unlist(data.frame(cc_per)))
df$age <- sapply(strsplit(df$exp, split = "_"),"[[",1)
df[df$celltype %in% c("Radial Glia","CycProg","Early Neurons","Retinal Neurons","Retinal Epi","Schwann Cells"), "lineage"] <- "ecto"
df[df$celltype %in% c("Definitive Endoderm","Foregut Epi","Airway Epi","Pancreatic Prog","Mid/Hindgut Epi"), "lineage"] <- "endo"
df[df$celltype %in% c("Immune","Adipogenic MSC/Fib",
                      "Chondrogenic MSC/Fib","MSC/Fib","Cycling MSC/Fib","MyoFib",
                      "Muscle Prog","Cardiac/Skeletal Muscle","Pericytes","Smooth Muscle","Kidney Prog"), "lineage"] <- "meso"

ggplot(data=df, aes(x=factor(age, levels = c("wk4","wk6","wk8","wk10")), y=per, group=celltype)) +
  geom_line(aes(color=celltype))+
  geom_point()

# try summary plot

#df <- data_summary(df, varname="per", 
#                    groupnames=c("age", "celltype"))
#head(df)

# plot a bit crowded, let's plot for each lineage?

df_ecto <- data_summary(df[df$lineage == "ecto",], varname="per", 
                    groupnames=c("age", "celltype"))
df_meso <- data_summary(df[df$lineage == "meso",], varname="per", 
                        groupnames=c("age", "celltype"))
df_endo <- data_summary(df[df$lineage == "endo",], varname="per", 
                        groupnames=c("age", "celltype"))

pdf("TS_all_RNA_celltype_line_bylineage_ecto_20230124.pdf")
ggplot(df_ecto, aes(x=factor(age, levels = c("wk4","wk6","wk8","wk10")), y=per, group=celltype, color=celltype)) + 
  geom_errorbar(aes(ymin=per-sd, ymax=per+sd), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point() + theme_classic()
dev.off()

pdf("TS_all_RNA_celltype_line_bylineage_endo_20230124.pdf")
ggplot(df_endo, aes(x=factor(age, levels = c("wk4","wk6","wk8","wk10")), y=per, group=celltype, color=celltype)) + 
  geom_errorbar(aes(ymin=per-sd, ymax=per+sd), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point() + theme_classic()
dev.off()

pdf("TS_all_RNA_celltype_line_bylineage_meso_20230124.pdf")
ggplot(df_meso, aes(x=factor(age, levels = c("wk4","wk6","wk8","wk10")), y=per, group=celltype, color=celltype)) + 
  geom_errorbar(aes(ymin=per-sd, ymax=per+sd), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point() + theme_classic()
dev.off()

# Use position_dodge to move overlapped errorbars horizontally
#pdf("TS_all_RNA_celltype_line_20230124.pdf")
#ggplot(df, aes(x=factor(age, levels = c("wk4","wk6","wk8","wk10")), y=per, group=celltype, color=celltype)) + 
#  geom_errorbar(aes(ymin=per-sd, ymax=per+sd), width=.1, 
#                position=position_dodge(0.05)) +
#  geom_line() + geom_point() + theme_classic()
#dev.off()

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


### separately save the cell barcode for each cell into a file (cell barcodes per cell type per mouse per batch)

# may need to collapse some cell types
obj.int.allrna.sub$collapsed_annot_v1 <- factor(obj.int.allrna.sub$annot_v1, levels = unique(c(levels(obj.int.allrna.sub$annot_v1), "Neural-Progenitors","Neurons","Retinal Epi","SCP",
                                                                                           "Definitive Endoderm","Gut-Epi","Pancreatic Prog","Hematopoietic","Adipogenic MSC/Fib",
                                                                                           "Chondrogenic MSC/Fib","MSC/Fib","Cycling MSC/Fib","MyoFib",
                                                                                           "Muscle","Pericytes","Smooth Muscle","Kidney Prog" )))

#Neural Progenitors: RG + CycProg
obj.int.allrna.sub$collapsed_annot_v1[obj.int.allrna.sub$annot_v1 %in% c("Radial Glia","CycProg")] <- "Neural-Progenitors"
#Neurons: Early Neurons + Retinal Neurons
obj.int.allrna.sub$collapsed_annot_v1[obj.int.allrna.sub$annot_v1 %in% c("Early Neurons","Retinal Neurons")] <- "Neurons"
#SCP: Schwann + Melanobalsts
obj.int.allrna.sub$collapsed_annot_v1[obj.int.allrna.sub$annot_v1 %in% c("Schwann Cells","Melanoblasts")] <- "SCP"

#Gut Epithelial: Foregut Epi + Mid/Hindgut Epi + Airway Epi
obj.int.allrna.sub$collapsed_annot_v1[obj.int.allrna.sub$annot_v1 %in% c("Airway Epi","Foregut Epi","Mid/Hindgut Epi")] <- "Gut-Epi"

#Hematopoietic: HSC + immune + Erythrocyte
obj.int.allrna.sub$collapsed_annot_v1[obj.int.allrna.sub$annot_v1 %in% c("HSC","Immune","Erythrocyte")] <- "Hematopoietic"
#Muscle: Muscle Prog + Cardiac/Skeletal Muscle
obj.int.allrna.sub$collapsed_annot_v1[obj.int.allrna.sub$annot_v1 %in% c("Muscle Prog","Cardiac/Skeletal Muscle")] <- "Muscle"


# fix cell type names
obj.int.allrna.sub$collapsed_annot_v1 <- str_replace_all(string = obj.int.allrna.sub$collapsed_annot_v1, pattern = "/", replacement = "-")
obj.int.allrna.sub$collapsed_annot_v1 <- str_replace_all(string = obj.int.allrna.sub$collapsed_annot_v1, pattern = fixed(" "), replacement = "-")
obj.int.allrna.sub$collapsed_annot_v1 <- factor(obj.int.allrna.sub$collapsed_annot_v1, levels = c("Neural-Progenitors","Neurons","Retinal-Epi","SCP",
                                                                                                  "Definitive-Endoderm","Gut-Epi","Pancreatic-Prog","Hematopoietic","Adipogenic-MSC-Fib",
                                                                                                  "Chondrogenic-MSC-Fib","MSC-Fib","Cycling-MSC-Fib","MyoFib",
                                                                                                  "Muscle","Pericytes","Smooth-Muscle","Kidney-Prog"))

table(obj.int.allrna.sub$collapsed_annot_v1,obj.int.allrna.sub$age_mouse_replicate )[,cc_per_order]

# save the seurat object to RNA editing folder, will switch to RNA editing project
save(obj.int.allrna.sub, file = "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_logNorm_int_annot_v1_col_v1_20230124.rda")

for (i in names(table(obj.int.allrna.sub$collapsed_annot_v1))){
  for (j in names(table(obj.int.allrna.sub$age_mouse_replicate))){
    if (length (row.names(obj.int.allrna.sub@meta.data[obj.int.allrna.sub$collapsed_annot_v1 == i & obj.int.allrna.sub$age_mouse_replicate == j, ])) != 0){
      write.table(paste0("CB:Z:",sapply(strsplit(row.names(obj.int.allrna.sub@meta.data[obj.int.allrna.sub$collapsed_annot_v1 == i & obj.int.allrna.sub$age_mouse_replicate == j, ]), split = "_"), "[[",5)), 
                  paste0('/media/Scratch_SSD_Voyager/sammi/RNA_editing/cb_listsTS/v1/',i,'_',j,'.txt'), col.names = FALSE, row.names = FALSE, 
                  quote = FALSE, sep = ',')
    }
  }
}

### view transcriptomic similarities between consecutive time points

load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_logNorm_int_annot_v1_col_v1_20230124.rda")
obj.int.allrna.sub

obj.int.allrna.sub <- FindVariableFeatures(obj.int.allrna.sub)
genes.use <- unique(c(as.character(YanMarkersPub), VariableFeatures(obj.int.allrna.sub), VariableFeatures(obj.int.allrna.sub)))
genes.use <- Reduce(intersect, list(genes.use, rownames(GetAssayData(obj.int.allrna.sub, assay = "integrated", slot = "scale.data")), 
                                    rownames(GetAssayData(obj.int.allrna.sub, assay = "integrated"))))
length(genes.use) # 1718
Tera.r.assay <- GetAssayData(obj.int.allrna.sub, assay = "integrated", slot = "scale.data")
cor.result <- MapClustersCor(query.data = Tera.r.assay, query.clusters = obj.int.allrna.sub$age_mouse_replicate, 
                             train.data = GetAssayData(obj.int.allrna.sub, assay = "integrated", slot = "scale.data"),
                             train.clusters = obj.int.allrna.sub$age_mouse_replicate,
                             genes.use = genes.use, metric = "pearson")
amr_order <- c("wk4_M1_1","wk4_M1_2","wk6_M15_1","wk6_M15_2","wk6_M9_1","wk8_M8_1","wk8_M8_2","wk8_M10_1","wk10_A_2","wk10_A_3","wk10_A_4","wk10_B_1","wk10_B_3")
ggHeat(t(cor.result$cluster.correlations[amr_order,amr_order]), x.lab.size = 11, y.lab.size = 11)

apply(cor.result$cluster.correlations[amr_order,amr_order], 2, max)

pdf("/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/TS_all_RNA_tpsamples_self_correlation_20230206.pdf", height = 7, width = 8)
ggHeat(t(cor.result$cluster.correlations[amr_order,amr_order]), x.lab.size = 11, y.lab.size = 11)
dev.off()


# check some more markers for MSC/Fib
load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_logNorm_int_annot_v1_col_v1_20230124.rda")
obj.int.allrna.sub

DefaultAssay(obj.int.allrna.sub) <- "RNA"

Stromal_Markers_info <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/Stromal_Annotation_Marker_Table.txt', sep = "\t", quote = "", header = 1)
Stromal_Markers <- noquote(Stromal_Markers_info$Additional.Markers)
gsub("\"", ",",paste(Stromal_Markers, collapse = ''))
Stromal_Markers <- strsplit(gsub("\"", ",",paste(Stromal_Markers, collapse = '')), split = ", |,,|,|,  ")[[1]]
Stromal_Markers <- Stromal_Markers[2:113]
length(Stromal_Markers)
length(unique(Stromal_Markers))
DotPlot(obj.int.allrna.sub, features = unique(Stromal_Markers))

pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_MSC_markers_20240102.pdf', width = 18, height = 12)
DotPlot(obj.int.allrna.sub, features = unique(Stromal_Markers)) + RotatedAxis()
dev.off()

pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_MSC_additional_markers_annotated_20240102.pdf', width = 18, height = 12)
DotPlot(obj.int.allrna.sub, features = unique(Stromal_Markers), group.by = "annot_v1") + RotatedAxis()
dev.off()

Fib_Markers_info <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/Human_Kidney_Adv_Fibroblast_Top_Markers.txt', sep = "\t", quote = "", header = 1)
pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_aFIBRSPO3_markers_annotated_20240102.pdf', width = 18, height = 12)
DotPlot(obj.int.allrna.sub, features = unique(Fib_Markers_info[["aFIB.RSPO3"]]), group.by = "annot_v1") + RotatedAxis()
dev.off()

pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_aFIBCD34_markers_annotated_20240102.pdf', width = 18, height = 12)
DotPlot(obj.int.allrna.sub, features = unique(Fib_Markers_info[["aFIB.CD34"]]), group.by = "annot_v1") + RotatedAxis()
dev.off()

pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_CMMYOF_markers_annotated_20240102.pdf', width = 18, height = 12)
DotPlot(obj.int.allrna.sub, features = unique(Fib_Markers_info[["C.M.MYOF"]]), group.by = "annot_v1") + RotatedAxis()
dev.off()

# check some more markers for MSC/Fib
load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_logNorm_int_annot_v1_col_v1_20230124.rda")
obj.int.allrna.sub

DefaultAssay(obj.int.allrna.sub) <- "RNA"

pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_markers_annotated_1_20240111.pdf', width = 18, height = 12)
DotPlot(obj.int.allrna.sub, features = c("GFAP", "MAP2", "PCDHA1", "GAD1"), group.by = "annot_v1") + RotatedAxis()
dev.off()

SaveH5Seurat(obj.int.allrna.sub, filename = "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_logNorm_int_annot_v1_col_v1_20230124.h5Seurat")
Convert("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_logNorm_int_annot_v1_col_v1_20230124.h5Seurat", dest = "h5ad")

# check some more markers for neural cells
load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_logNorm_int_annot_v1_col_v1_20230124.rda")
obj.int.allrna.sub

DefaultAssay(obj.int.allrna.sub) <- "RNA"

Idents(obj.int.allrna.sub) <- obj.int.allrna.sub$annot_v1
pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_markers_annotated_1_20240813.pdf', width = 18, height = 12)
cowplot:::plot_grid(DotPlot(obj.int.allrna.sub, features = c("OTX2", "GBX2"), group.by = "annot_v1") + RotatedAxis(),
FeaturePlot(obj.int.allrna.sub, features = c("OTX2", "GBX2"), label = T), rel_widths = c(1,2))
dev.off()

obj.int.allrna.sub@meta.data[WhichCells(object=obj.int.allrna.sub, expression = OTX2 > 1 & GBX2 > 1),]
obj.int.allrna.sub@meta.data[WhichCells(object=obj.int.allrna.sub, expression = OTX2 > 1),]
obj.int.allrna.sub@meta.data[WhichCells(object=obj.int.allrna.sub, expression = GBX2 > 1),]

# check some more markers for neural cells
load('/media/NAS1/Sammi_NAS1/SO/TS_neu_only_20230215.rda')
obj.int.allrna.sub_neu_sub

DefaultAssay(obj.int.allrna.sub_neu_sub) <- "RNA"

Idents(obj.int.allrna.sub_neu_sub) <- obj.int.allrna.sub_neu_sub$annot_l2_v1
pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_neuroectoderm_RNA_markers_annotated_2_20240813.pdf', width = 18, height = 12)
cowplot:::plot_grid(DotPlot(obj.int.allrna.sub_neu_sub, features = c("OTX2", "GBX2"), group.by = "annot_l2_v1") + RotatedAxis(),
FeaturePlot(obj.int.allrna.sub_neu_sub, features = c("OTX2", "GBX2"), label = T), rel_widths = c(1,2))
dev.off()

obj.int.allrna.sub_neu_sub@meta.data[WhichCells(object=obj.int.allrna.sub_neu_sub, expression = MAFB > 1 & HNF1B > 1),]


### scripts for correlating teratoma time series data with Human Fetal Atlas data


library(Seurat)
library(dplyr)


### transcriptomic similarity

# fetal data from Shendure et. al
HF_Cbrum <- readRDS('/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/ref/Fetal_Atlas/Cerebrum_gene_count.RDS')
dim(HF_Cbrum)

# switch gene name
HF_gene_name <- readRDS('/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/ref/Fetal_Atlas/df_gene.RDS')
length(unique(HF_gene_name$gene_id)) ; length(unique(HF_gene_name$gene_short_name))
HF_gene_name <- HF_gene_name[match(HF_gene_name$gene_id, row.names(HF_Cbrum)),] # make sure order matches
row.names(HF_Cbrum) <- HF_gene_name$gene_short_name
rm(HF_gene_name)

# check meta data
HF_cell_meta <- readRDS('/media/Scratch_SSD_Voyager/sammi/Teratoma_TS/ref/Fetal_Atlas/df_cell.RDS')
dim(HF_cell_meta)
dim(HF_cell_meta[complete.cases(HF_cell_meta),])
HF_cell_meta <- HF_cell_meta[complete.cases(HF_cell_meta) & HF_cell_meta$sample %in% colnames(HF_Cbrum),]
dim(HF_cell_meta)
row.names(HF_cell_meta) <- HF_cell_meta$sample

HF_Cbrum <- HF_Cbrum[, colnames(HF_Cbrum) %in% row.names(HF_cell_meta)]
dim(HF_Cbrum)
HF_neu <- CreateSeuratObject(counts = HF_Cbrum, meta.data = HF_cell_meta, assay = "RNA", project = "HF", 
                             min.cells = round(0.001*ncol(HF_Cbrum)), min.features = 200)
head(HF_neu@meta.data)

JSHF_ct_collapse <- read.table("/media/NAS1/Sammi_NAS1/TS_Archive/20220406_Teratoma_TS/JSHF_Collapsing.csv", header = TRUE, sep = ",")
cl.replace <- JSHF_ct_collapse$New.Ident
names(cl.replace) <- JSHF_ct_collapse$Orig.Ident
list_cluster <- factor(HF_neu$Main_cluster_name)
list_cluster <- plyr::revalue(list_cluster, replace = cl.replace)
HF_neu$Collapsed_CT <- list_cluster
table(HF_neu$Collapsed_CT)

Idents(HF_neu) <- HF_neu$Collapsed_CT
HF_neu <- subset(HF_neu, ident = c("Astrocytes","Ex.Neurons","In.Neurons","Neurons"))
HF_neu$Collapsed_CT <- droplevels(HF_neu$Collapsed_CT)

# down sample HF_neu to pcd 89, 115, 122
Idents(HF_neu) <- HF_neu$Development_day
HF_neu <- subset(HF_neu, ident = c(89, 115,122))
Idents(HF_neu) <- HF_neu$Collapsed_CT
HF_neu$Collapsed_CT <- droplevels(HF_neu$Collapsed_CT)

# down sample HF_neu to equal # of cells for each time point
HF_neu <- SplitObject(HF_neu, split.by = "Development_day")
sapply(HF_neu, dim)
for (i in 1:length(HF_neu)){
  HF_neu[[i]] <- HF_neu[[i]][, sample(colnames(HF_neu[[i]]), size = min(sapply(HF_neu, dim)[2,]), replace = F )]
}
HF_neu <- merge(HF_neu[[1]], y = c(HF_neu[[2]], HF_neu[[3]], HF_neu[[4]], HF_neu[[5]], HF_neu[[6]]), 
                               merge.data = T, merge.dr = "umap")
table(HF_neu$Development_day,HF_neu$Collapsed_CT)


HF_neu <- NormalizeData(HF_neu)
HF_neu <- FindVariableFeatures(HF_neu, nfeatures = 3000)
HF_neu <- ScaleData(HF_neu, features = VariableFeatures(HF_neu))

load('/media/NAS1/Sammi_NAS1/SO/TS_neu_only_20230215.rda')
table(obj.int.allrna.sub_neu_sub$annot_l2_v1, obj.int.allrna.sub_neu_sub$teratoma)
table(obj.int.allrna.sub_neu_sub$age)
Idents(obj.int.allrna.sub_neu_sub) <- obj.int.allrna.sub_neu_sub$annot_l2_v1
obj.int.allrna.sub_neu_sub <- subset(obj.int.allrna.sub_neu_sub, ident = c("CycProg","RG/Astro","Radial Glia","pre-OPC",
"Neural Precursor Cells","Intermediate Progenitors","Early Ex Neu","Dopaminergic Neu","Early Differentiating Neu","Differentiating Ex Neu",
"Differentiating Inh Neu","Inhib Neu"))
obj.int.allrna.sub_neu_sub$annot_l2_v1 <- droplevels(obj.int.allrna.sub_neu_sub$annot_l2_v1)
obj.int.allrna.sub_neu_sub

obj.int.allrna.sub_neu_sub <- FindVariableFeatures(obj.int.allrna.sub_neu_sub, nfeatures = 3000)
obj.int.allrna.sub_neu_sub <- ScaleData(obj.int.allrna.sub_neu_sub, features = VariableFeatures(obj.int.allrna.sub_neu_sub))

YanMarkersPub_neu <- c("VIM","SOX2","HES5","HMGB2","DCX","MAP2","SLC17A6","PAX6","NEUROD6","SLC17A7")
genes.use <- unique(c(as.character(YanMarkersPub_neu), VariableFeatures(HF_neu), VariableFeatures(obj.int.allrna.sub_neu_sub)))
genes.use <- Reduce(intersect, list(genes.use, rownames(GetAssayData(HF_neu, assay = "RNA", slot = "scale.data")), 
                                    rownames(GetAssayData(obj.int.allrna.sub_neu_sub, assay = "RNA", slot = "scale.data"))))
length(genes.use) # 1102
Tera.r.assay <- GetAssayData(obj.int.allrna.sub_neu_sub, assay = "RNA", slot = "scale.data")

cor.result <- cellMapper:::MapClustersCor(query.data = Tera.r.assay, query.clusters = obj.int.allrna.sub_neu_sub$annot_l2_v1, 
                             train.data = GetAssayData(HF_neu, assay = "RNA", slot = "scale.data"),
                             train.clusters = HF_neu$Collapsed_CT,
                             genes.use = genes.use)

pdf("Figure2_TimeSeries_R1/TS_neuroectoderm_correlating_shendurecerebellum_celltypesonly_notp_ds_v2_20230906.pdf", height = 7, width = 5)
ggHeat(t(cor.result$cluster.correlations[complete.cases(cor.result$cluster.correlations),]), x.lab.size = 11, y.lab.size = 11) + 
  labs(title='Correlation between Cell Types from Teratoma Neuroectoderm and Fetal Cortex', x='Fetal Cortex Cell Types', y='Teratoma Time Series Neuroectoderm Cell Types')
dev.off()

obj.int.allrna.sub_neu_sub$age_ct <- paste0(obj.int.allrna.sub_neu_sub$age, "_", obj.int.allrna.sub_neu_sub$annot_l2_v1)
HF_neu$age_ct <- paste0(HF_neu$Development_day, "_", HF_neu$Collapsed_CT)

cor.result <- cellMapper:::MapClustersCor(query.data = Tera.r.assay, query.clusters = obj.int.allrna.sub_neu_sub$age_ct, 
                             train.data = GetAssayData(HF_neu, assay = "RNA", slot = "scale.data"),
                             train.clusters = HF_neu$age_ct,
                             genes.use = genes.use)

pdf("Figure2_TimeSeries_R1/TS_neuroectoderm_correlating_shendurecerebellum_celltypesonly_notp_ds_v2_20230906.pdf", height = 7, width = 5)
ggHeat(t(cor.result$cluster.correlations[complete.cases(cor.result$cluster.correlations),]), x.lab.size = 11, y.lab.size = 11, cluster = "both") + 
  labs(title='Correlation between Cell Types from Teratoma Neuroectoderm and Fetal Cortex', x='Fetal Cortex Cell Types', y='Teratoma Time Series Neuroectoderm Cell Types')
dev.off()



### all cell type AEI correlation

aei_JS <- read.table('Figure_HFA_R1/aei_v2_all.txt', sep = ",", quote = "", header = T)
head(aei_JS)

agg_JS_tbl <- aei_JS %>% group_by(CellType) %>% 
  summarise(mean_AEI=mean(A2GEditingIndex),
            .groups = 'drop')

aei_ter <- read.table('Figure2_TimeSeries_R1/TS_allRNA_AEI_byCT_byExp_20230510.txt', header = T, sep =',')
head(aei_ter)
aei_ter$Sample <- row.names(aei_ter)
aei_ter$Subject <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 2)
#aei_ter$CellType <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 1)
#aei_ter$Time <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 3)
aei_ter <- aei_ter[,c("Sample","A2GEditingIndex","Subject","CellType","Age")]
aei_ter <- aei_ter[aei_ter$Age == "wk10",]

agg_ter_tbl <- aei_ter %>% group_by(CellType) %>% 
  summarise(mean_AEI=mean(A2GEditingIndex),
            .groups = 'drop')

# define cell type matching relationships
df_matchct <- data.frame("ter_ct" = c("Neural-Progenitors","Neurons","Neurons","Neurons","Neurons",
#"Neural-Progenitors","Neurons","Neurons","Neurons","Neurons","Retinal-Epi","Retinal-Epi",
"SCP","SCP",
                    #"Definitive-Endoderm",
                    "Gut-Epi",
                    #"Pancreatic-Prog",
                    "Hematopoietic","Hematopoietic","Hematopoietic","Hematopoietic",
                    "Adipogenic-MSC-Fib","Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Pericytes",
                    "Muscle","Muscle","Smooth-Muscle",
                    "Kidney-Prog"),
                    "ter_AEI" = NA,
                    "js_ct" = c("Astrocytes","Neurons","Ex.Neurons","In.Neurons","In.Interneurons",
                      #"Retinal.Prog","Retinal.Neurons","Retinal.Interneurons",
                      #"Ganglion.Cells","Horizontal.Cells","Retinal.Pig",
                      #"Retinal.Epi",
                      "Schwann.Cells","Glia",
                      
                      #"Foregut.Epi",
                      "Gut.Epi","Vascular.Endo",
                      "Microglia","Immune.Cells","Blood.Cells",
                      #"MSC","MSC","MSC","MSC","MSC","MSC",
                      "Stromal.Cells","Stromal.Cells","Stromal.Cells","Stromal.Cells","Stromal.Cells","Stromal.Cells",
                      "Skeletal.Muscle.Cells", "Cardiac.Muscle.Cells",
                      "Smooth.Muscle.Cells",
                      "Kidney"),
                      "js_AEI" = NA)


df_matchct$js_ct <- stringr:::str_replace_all(df_matchct$js_ct, "[.]", "-")

which(!df_matchct$ter_ct %in% agg_ter_tbl$CellType)
which(!df_matchct$js_ct %in% agg_JS_tbl$CellType)

for (i in 1:nrow(df_matchct)){
  df_matchct[i,"ter_AEI"] <- agg_ter_tbl[agg_ter_tbl$CellType == df_matchct[i,"ter_ct"], "mean_AEI"]
  df_matchct[i,"js_AEI"] <- agg_JS_tbl[agg_JS_tbl$CellType == df_matchct[i,"js_ct"], "mean_AEI"]

}

write.table(df_matchct, file = "Figure2_TimeSeries_R1/TS_allRNA_wk10_AEI_corr_JS_AEI.txt", sep = ",", quote = F, col.names = T)

### neuroectoderm AEI correlation

aei_JS <- read.table('Figure_HFA_R1/aei_v2_all.txt', sep = ",", quote = "", header = T)
head(aei_JS)

agg_JS_tbl <- aei_JS %>% group_by(CellType,Development_day) %>% 
  summarise(mean_AEI=mean(A2GEditingIndex),
            .groups = 'drop')

aei_ter <- read.table('AEI_TS_v2neuroectoderm/output_index_20230605/AEI_TS_v2neuroectoderm_outputindex/EditingIndex.csv', header = T, sep =',')
head(aei_ter)
row.names(aei_ter) <- aei_ter$Sample
aei_ter$CellType <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 1)
aei_ter$Teratoma <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 3)
aei_ter$Age <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 2)
aei_ter <- aei_ter[,c("Sample","A2GEditingIndex","Age","CellType","Teratoma")]
# drop M9
aei_ter <- aei_ter[aei_ter$Teratoma != "M9",]

agg_ter_tbl <- aei_ter %>% group_by(CellType,Age) %>% 
  summarise(mean_AEI=mean(A2GEditingIndex),
            .groups = 'drop')

# define cell type matching relationships
df_matchct <- data.frame("ter_ct" = paste0(rep(c("wk4","wk6","wk8","wk10"), each = 15),"_",rep(c("CycProg","RG-Astro","Radial-Glia","pre-OPC","Neural-Precursor-Cells",
"Intermediate-Progenitors","Early-Ex-Neu","Dopaminergic-Neu","Differentiating-Ex-Neu","Differentiating-Inh-Neu",
"Inhib-Neu","Dopaminergic-Neu","Differentiating-Ex-Neu","Differentiating-Inh-Neu",
"Inhib-Neu"), n = 4)),
                    "ter_AEI" = NA,
                    "js_ct" = paste0(rep(c(89,89,115,122), each = 15), "_", rep(c("Astrocytes","Astrocytes","Astrocytes","Astrocytes","Ex-Neurons",
                    "Ex-Neurons","Ex-Neurons","In-Neurons","In-Neurons","In-Neurons","In-Neurons",
                    "Neurons","Neurons","Neurons","Neurons"), n = 4)),
                      "js_AEI" = NA)


df_matchct$js_ct <- stringr:::str_replace_all(df_matchct$js_ct, "[.]", "-")

#which(!df_matchct$ter_ct %in% agg_ter_tbl$CellType)
#which(!df_matchct$js_ct %in% agg_JS_tbl$CellType)

for (i in 1:nrow(df_matchct)){
    if (dim(agg_ter_tbl[agg_ter_tbl$CellType == strsplit(df_matchct[i,"ter_ct"], split = "_")[[1]][2] & agg_ter_tbl$Age == strsplit(df_matchct[i,"ter_ct"], split = "_")[[1]][1], "mean_AEI"])[1] > 0) {
        df_matchct[i,"ter_AEI"] <- agg_ter_tbl[agg_ter_tbl$CellType == strsplit(df_matchct[i,"ter_ct"], split = "_")[[1]][2] & agg_ter_tbl$Age == strsplit(df_matchct[i,"ter_ct"], split = "_")[[1]][1], "mean_AEI"]
    }
    if (dim(agg_JS_tbl[agg_JS_tbl$CellType == strsplit(df_matchct[i,"js_ct"], split = "_")[[1]][2] & agg_JS_tbl$Development_day == strsplit(df_matchct[i,"js_ct"], split = "_")[[1]][1], "mean_AEI"])[1] > 0) {
        df_matchct[i,"js_AEI"] <- agg_JS_tbl[agg_JS_tbl$CellType == strsplit(df_matchct[i,"js_ct"], split = "_")[[1]][2] & agg_JS_tbl$Development_day == strsplit(df_matchct[i,"js_ct"], split = "_")[[1]][1], "mean_AEI"]
    }

}




aei_JS <- read.table('Figure_HFA_R1/aei_v2_all.txt', sep = ",", quote = "", header = T)
head(aei_JS)

agg_JS_tbl <- aei_JS %>% group_by(Development_day) %>% 
  summarise(mean_AEI=mean(A2GEditingIndex),
            .groups = 'drop')
agg_JS_tbl[agg_JS_tbl$Development_day %in% c(89,115,122),]

aei_ter <- read.table('AEI_TS_v2neuroectoderm/output_index_20230605/AEI_TS_v2neuroectoderm_outputindex/EditingIndex.csv', header = T, sep =',')
head(aei_ter)
row.names(aei_ter) <- aei_ter$Sample
aei_ter$CellType <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 1)
aei_ter$Teratoma <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 3)
aei_ter$Age <- sapply(strsplit(row.names(aei_ter), split = "_"), "[[", 2)
aei_ter <- aei_ter[,c("Sample","A2GEditingIndex","Age","CellType","Teratoma")]
# drop M9
aei_ter <- aei_ter[aei_ter$Teratoma != "M9",]
aei_ter

# remove low cell count samples
cc_ter <- read.table('Figure2_TimeSeries_R1/TS_all_RNA_neuroectoderm_celltype_dist_20230428.txt', header = T, sep = ",", quote = "")
row.names(cc_ter) <- stringr:::str_replace_all(row.names(cc_ter), " ", "-")
row.names(cc_ter) <- stringr:::str_replace_all(row.names(cc_ter), "/", "-")
colnames(cc_ter) <- paste0("wk",sapply(strsplit(colnames(cc_ter), split = "X"), "[[", 2))

rows_to_rm <- c()
for (i in row.names(aei_ter)) {
  ct_tmp <- strsplit(i, split = "_")[[1]][1]
  sample_tmp <- paste0(strsplit(i, split = "_")[[1]][2], "_", strsplit(i, split = "_")[[1]][3], "_", strsplit(i, split = "_")[[1]][4])
  if (cc_ter[ct_tmp, sample_tmp] < 10){
    rows_to_rm <- c(rows_to_rm, i)
  }
}
aei_ter <- aei_ter[!(row.names(aei_ter) %in% rows_to_rm),]

agg_ter_tbl <- aei_ter %>% group_by(Age) %>% 
  summarise(mean_AEI=mean(A2GEditingIndex),
            .groups = 'drop')
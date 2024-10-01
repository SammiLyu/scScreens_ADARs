### generate plots for Figure 3 for both CO & H1

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)


### for Cortical Organoid

aei_co <- read.table("AEI_CO_v2/output_index_20230627/AEI_CO_v2_outputindex/EditingIndex.csv", header = 1, sep = ",")

cluster_levels <- c("Progenitor","IP","Glia","Glutamatergic","GABAergic")
cluster_col <- c("#442288","#6CA2EA","#B5D33D","#A62639","#EB7D5B")

names(cluster_col) <- cluster_levels

aei_co$CellType <- sapply(strsplit(aei_co$Sample, split = "_"),"[[",2)
aei_co$CellType <- factor(aei_co$CellType, levels = cluster_levels)
aei_co$Age <- sapply(strsplit(aei_co$Sample, split = "_"),"[[",1)
aei_co$Age <- factor(aei_co$Age, levels = c('1mo','3mo','6mo','10mo'))


co_count_table <- read.table("Figure_CO_R1/CO_4co_celltype_count.txt", quote = "", sep = ",", header = 1)
# remove entry 1mo_GABAergic, 6mo_GABAergic, 10mo_Progenitor for cell count <10?
#aei_co <- aei_co[!aei_co$Sample %in% c("1mo_GABAergic","6mo_GABAergic","10mo_Progenitor")]

pdf('Figure_CO_R1/CO_4co_AEI_20230717.pdf', height = 8, width = 8)
ggplot(data=aei_co, aes(x=Age, y=A2GEditingIndex, colour=CellType)) +
  geom_point(stat="identity", position=position_dodge(width = 0.4), size = 3) + RotatedAxis() +
  scale_color_manual(values=cluster_col) + theme_classic()
dev.off()

pdf('Figure_CO_R1/CO_4co_AEI_10mo_20230717.pdf', height = 8, width = 8)
ggplot(data=aei_co[aei_co$Age == "10mo",], aes(x=CellType, y=A2GEditingIndex, fill=CellType)) +
  geom_bar(stat="identity") + #RotatedAxis() +
  scale_fill_manual(values=cluster_col) + theme_classic()
dev.off()

adar_tpm <- read.table("fc_CO_v1/adar_exp_tpm.txt", header = 1, sep = "\t")
adar_tpm$CellType <- sapply(strsplit(adar_tpm$Sample, split = "_"),"[[",2)
adar_tpm$Age <- sapply(strsplit(adar_tpm$Sample, split = "_"),"[[",3)
adar_tpm$Sample <- sapply(strsplit(adar_tpm$Sample, split = "20230627_"),"[[",2)
adar_tpm$Sample <- paste0(sapply(strsplit(adar_tpm$Sample, split = "_"),"[[",2),"_", sapply(strsplit(adar_tpm$Sample, split = "_"),"[[",1))

ggplot(data=adar_tpm, aes(x=Age, y=TPM, colour=CellType)) +
  geom_point(stat="identity", position=position_dodge(width = 0.4), size = 3) + RotatedAxis() +
  scale_color_manual(values=cluster_col) + theme_classic()

df_aei_by_tpm <- data.frame(matrix(NA, nrow = nrow(aei_co), ncol = 6, dimnames = list(aei_co$Sample, c("AEI_A2G","ADAR_TPM","ADARB1_TPM","ADARB2_TPM","cell_type","age"))))
for (i in row.names(df_aei_by_tpm)){
  df_aei_by_tpm[i,"AEI_A2G"] <- aei_co[aei_co$Sample == i,"A2GEditingIndex"]
  df_aei_by_tpm[i,"ADAR_TPM"] <- adar_tpm[adar_tpm$Sample == i & adar_tpm$Gene == 'ADAR',"TPM"]
  df_aei_by_tpm[i,"ADARB1_TPM"] <- adar_tpm[adar_tpm$Sample == i & adar_tpm$Gene == 'ADARB1',"TPM"]
  df_aei_by_tpm[i,"ADARB2_TPM"] <- adar_tpm[adar_tpm$Sample == i & adar_tpm$Gene == 'ADARB2',"TPM"]
}

df_aei_by_tpm$cell_type <- sapply(strsplit(row.names(df_aei_by_tpm), split = "_"),"[[",2)
df_aei_by_tpm$age <- sapply(strsplit(row.names(df_aei_by_tpm), split = "_"),"[[",1)

pdf("Figure_CO_R1/CO_4co_AEI_by_TPM_20230717.pdf", width = 6, height = 6)
ggplot(df_aei_by_tpm, aes(x = log2(ADAR_TPM), y = AEI_A2G)) + geom_jitter(size=3, aes(color = cell_type)) + scale_color_manual(values = cluster_col) + 
  #geom_smooth(method='lm', fullrange=TRUE) + 
  #geom_text(x = 6, y = 0.65, label = eq(log2(df_aei_by_tpm$ADAR_TPM+0.01),df_aei_by_tpm$AEI_A2G), parse = TRUE) +
  ylim(0.2, 0.7) + xlim(2,5.5) + theme_classic()
dev.off()

p <- list()
for (i in names(table(df_aei_by_tpm$age))){
    p[[i]] <- ggplot(df_aei_by_tpm[df_aei_by_tpm$age == i,], aes(x = log2(ADAR_TPM), y = AEI_A2G)) + geom_jitter(size=3, aes(color = cell_type)) + scale_color_manual(values = cluster_col) + 
  #geom_smooth(method='lm', fullrange=TRUE) + 
  #geom_text(x = 6, y = 0.65, label = eq(log2(df_aei_by_tpm$ADAR_TPM+0.01),df_aei_by_tpm$AEI_A2G), parse = TRUE) +
  ylim(0.2, 0.7) + xlim(2,5.5) + theme_classic() + ggtitle(i)
}

pdf("Figure_CO_R1/CO_4co_AEI_by_TPM_bymo_20230717.pdf", width = 6, height = 6)
p
dev.off()


### for H1

aei_H1 <- read.table("AEI_H1_v1/output_index_20230518/AEI_H1_v1_outputindex/EditingIndex.csv", header = 1, sep = ",")

cluster_levels <- c("Radial Glia","CycProg","Intermediate Progenitor","Glutamatergic","GABAergic","Early Neurons","Retinal Neurons","Retinal Epi","Schwann Cells","Melanoblasts",
                    "Foregut Epi","Mid/Hindgut Epi","Airway Epi","HSC","Immune","Erythrocyte",
                    "Adipogenic MSC/Fib","Chondrogenic MSC/Fib","Cycling MSC/Fib","MSC/Fib","MyoFib",
                    "Muscle Prog","Cardiac/Skeletal Muscle","Pericytes","Smooth Muscle",
                    "Kidney Prog")

cluster_col <- c("#4795F5","#0B62CD","#6CA2EA","#A62639","#EB7D5B","#7474FF","#DADAFF","#F2B77C","#9EEA00","#4FF4A2",
                 "#002700","#005A00","#00A700","#000032","#04409B","#4545FF",
                 "#0000FF","#FF0000","#00C000","#AD07E3","#000000",
                 "#B35A00","#FF9933","#A00000","#94641F","#610051")

names(cluster_col) <- cluster_levels

aei_H1$CellType <- sapply(strsplit(aei_H1$Sample, split = "_"),"[[",2)
aei_H1$CellType <- factor(aei_H1$CellType, levels = cluster_levels)

pdf('Figure_WT_R1/WT_4H1_AEI_20230717.pdf', height = 8, width = 8)
ggplot(data=aei_H1, aes(x=CellType, y=A2GEditingIndex, fill=CellType)) +
  geomgeom_bar_point(stat="identity") + #RotatedAxis() +
  scale_color_manual(values=cluster_col) + theme_classic()
dev.off()

adar_tpm <- read.table("fc_H1_v1/adar_exp_tpm.txt", header = 1, sep = "\t")
adar_tpm$CellType <- sapply(strsplit(adar_tpm$Sample, split = "_"),"[[",2)
adar_tpm$Sample <- sapply(strsplit(adar_tpm$Sample, split = "20230627_"),"[[",2)
adar_tpm$Sample <- paste0(sapply(strsplit(adar_tpm$Sample, split = "_"),"[[",2),"_", sapply(strsplit(adar_tpm$Sample, split = "_"),"[[",1))

df_aei_by_tpm <- data.frame(matrix(NA, nrow = nrow(aei_H1), ncol = 6, dimnames = list(aei_co$Sample, c("AEI_A2G","ADAR_TPM","ADARB1_TPM","ADARB2_TPM","cell_type","age"))))
for (i in row.names(df_aei_by_tpm)){
  df_aei_by_tpm[i,"AEI_A2G"] <- aei_H1[aei_H1$Sample == i,"A2GEditingIndex"]
  df_aei_by_tpm[i,"ADAR_TPM"] <- adar_tpm[adar_tpm$Sample == i & adar_tpm$Gene == 'ADAR',"TPM"]
  df_aei_by_tpm[i,"ADARB1_TPM"] <- adar_tpm[adar_tpm$Sample == i & adar_tpm$Gene == 'ADARB1',"TPM"]
  df_aei_by_tpm[i,"ADARB2_TPM"] <- adar_tpm[adar_tpm$Sample == i & adar_tpm$Gene == 'ADARB2',"TPM"]
}

df_aei_by_tpm$cell_type <- sapply(strsplit(row.names(df_aei_by_tpm), split = "_"),"[[",2)

pdf("Figure_WT_R1/WT_4H1_AEI_by_TPM_20230717.pdf", width = 6, height = 6)
ggplot(df_aei_by_tpm, aes(x = log2(ADAR_TPM), y = AEI_A2G)) + geom_jitter(size=3, aes(color = cell_type)) + scale_color_manual(values = cluster_col) + 
  #geom_smooth(method='lm', fullrange=TRUE) + 
  #geom_text(x = 6, y = 0.65, label = eq(log2(df_aei_by_tpm$ADAR_TPM+0.01),df_aei_by_tpm$AEI_A2G), parse = TRUE) +
  ylim(0.2, 0.7) + xlim(2,5.5) + theme_classic()
dev.off()


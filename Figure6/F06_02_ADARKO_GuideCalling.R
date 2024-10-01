### determine which guide the cells in ADAR-KO teratoma carries

library(dplyr)
library(stringr)

# working path
setwd('/media/Scratch_SSD_Voyager/sammi/RNA_editing/')
input_dir = '/media/NAS1/Sammi_NAS1/Cellranger_Output/ADAR_20220614/dual_ref/'
ter_list <- c("ter1","ter2","ter3","ter4")

#Load dataset
ter_KO <- list(Read10X(data.dir = paste0(input_dir, "/ter1/outs/filtered_feature_bc_matrix")),
               Read10X(data.dir = paste0(input_dir, "/ter2/outs/filtered_feature_bc_matrix")),
               Read10X(data.dir = paste0(input_dir, "/ter3/outs/filtered_feature_bc_matrix")),
               Read10X(data.dir = paste0(input_dir, "/ter4/outs/filtered_feature_bc_matrix")))
names(ter_KO) <- ter_list

### several visualizing
# # of total guide UMI for each cell
apply(ter_KO[[1]][[2]], 2, sum)
# fraction of guide UMI of top guide UMI
dim(apply(ter_KO[[1]][[2]], 2, function(x) x/sum(x))) # fraction for all guides, 20 by cells
colSums(apply(ter_KO[[1]][[2]], 2, function(x) x/sum(x))) # each column should sum up to 1 (or NaN for cells with no guide UMI at all)
apply(apply(ter_KO[[1]][[2]], 2, function(x) x/sum(x)), 2, max) # fraction of top guide

### distinguish the human and the mouse cells using cellranger output
ter_gem_cla <- list()
for (i in 1:length(ter_list)){
  ter_gem_cla[[i]] <- read.table(paste0(input_dir, ter_list[[i]], "/outs/analysis/gem_classification.csv"), sep = ",", header = TRUE)
  names(ter_gem_cla)[[i]] <- ter_list[[i]]
  row.names(ter_gem_cla[[i]]) <- paste0(names(ter_gem_cla)[[i]], "_", ter_gem_cla[[i]]$barcode)
}
length(ter_gem_cla)
head(ter_gem_cla[[1]],6)

rm(input_dir)

### select only singlet human cells
guide_count_matrix_cu <- list()
for (i in 1:length(ter_KO)){
  singlet_human <- which(paste0(names(ter_KO)[i], "_", colnames(ter_KO[[i]][[2]])) %in% row.names(ter_gem_cla[[i]][ter_gem_cla[[i]]$call == "GRCh38",]))
  print(length(singlet_human))
  guide_count_matrix_cu[[i]] <- ter_KO[[i]][[2]][, singlet_human]
  names(guide_count_matrix_cu)[i] <- names(ter_KO)[i]
}

sapply(guide_count_matrix_cu, dim)

sc_df <- list()
for (min_guide in seq(0,10,1)){
  for (i in 1:length(guide_count_matrix_cu)){
    results_temp <- gUMI_vis(guide_count_matrix_cu[[i]], names(guide_count_matrix_cu)[i], min_guide = min_guide)
    row.names(results_temp[[1]]) <- paste0(names(guide_count_matrix_cu)[i], "_", row.names(results_temp[[1]]))
    sc_df[[paste0(names(guide_count_matrix_cu)[i], "_", min_guide)]] <- results_temp
    rm(results_temp)
  }
  ter_guide <- rbind(sc_df[[paste0(names(guide_count_matrix_cu)[1], "_", min_guide)]][[1]], 
                     sc_df[[paste0(names(guide_count_matrix_cu)[2], "_", min_guide)]][[1]], 
                     sc_df[[paste0(names(guide_count_matrix_cu)[3], "_", min_guide)]][[1]], 
                     sc_df[[paste0(names(guide_count_matrix_cu)[4], "_", min_guide)]][[1]])
  #print(dim(ter_guide)) # 50339 cells total
  #head(ter_guide)
  
  ter_guide$teratoma <- sapply(strsplit(row.names(ter_guide), split = "_"),"[[",1)
  ter_guide$guide_info <- sapply(strsplit(ter_guide$g_id, split = c("[-_]+")),"[[",1)
  #write.table(ter_guide, paste0("varying_cutoff_max/guide_info_4tera_min", min_guide, "_v2.csv"),col.names = TRUE, row.names = TRUE, quote = FALSE, sep = ',')
}


pdf("Figure3_ADARKO_R1/guide_calling_scatter_plot_v1.pdf", width = 10, height = 6)
for (i in 1:length(sc_df)){
  print(sc_df[[i]][[2]])
}
dev.off()

pdf("Figure3_ADARKO_R1/scatter_and_hist_guide_calling_min_5_20230201.pdf", width = 19, height = 12)
plot_grid(plotlist = list(sc_df[["ter1_5"]][[2]], sc_df[["ter2_5"]][[2]], sc_df[["ter3_5"]][[2]], sc_df[["ter4_5"]][[2]]))
dev.off()

# a function that identifies the guide identity of a cell based on max gUMI identity,
# if a cell has more than 1 "max gUMI", note as multiplet
gUMI_vis <- function(guide_count_matrix, sample_name, min_guide){
  sc_max_temp <- apply(guide_count_matrix, 2, max)
  sc_df_temp <- data.frame(tot_UMI = apply(guide_count_matrix, 2, sum),
                           top_UMI = sc_max_temp)
  
  for (it in 1:length(sc_max_temp)){
    if (sc_df_temp[it, "tot_UMI"] >= min_guide) {
      g_id_it <- names(which(guide_count_matrix[,it] == sc_max_temp[it]))
      if (length(g_id_it) == 1){
        sc_df_temp[names(sc_max_temp)[it],"g_id"] <- g_id_it
      }
      else{
        sc_df_temp[names(sc_max_temp)[it],"g_id"] <- 'multiplet'
      }
    }
  }
  
  highlight_df <- sc_df_temp %>% 
    dplyr::filter(tot_UMI>min_guide)
  
  p1 <- ggplot(sc_df_temp, aes(x=log10(tot_UMI), y=top_UMI)) + geom_point(size = 0.4) + geom_point(data=highlight_df, 
                                                                                                   aes(x=log10(tot_UMI),y=top_UMI), 
                                                                                                   color='red', size = 0.4)+ theme_classic() +
    labs(title= paste0("scatter plot of g_UMI of ", sample_name, ' (# of guide UMI: >=', c(min_guide), ') (', c(dim(highlight_df)[1]), ' cells selected)'),
         x="log10(total # of guide UMI)", y = "top guide UMI count") + xlim(0,4)
  
  p2 <- ggplot(sc_df_temp, aes(x=top_UMI)) + geom_histogram(binwidth = 0.05, boundary = 0)  + theme_classic() + 
    labs(title = paste0("Histogram of top g_UMI count"), x = "top guide UMI count", y = "# of cells") + 
    coord_flip()
  
  sc_df_temp %>%
    count(cut_width(tot_UMI,30))
  
  p3 <- ggplot(sc_df_temp, aes(x=log10(tot_UMI))) + geom_histogram(binwidth = 0.05, boundary = 0)  + theme_classic() + 
    labs(title = paste0("Histogram of total guide UMI"), x = "log10(total # of guide UMI)", y = "# of cells") + xlim(0,4) + 
    geom_vline(xintercept = log10(min_guide), linetype="dotted", color = "red", size=1.5)
  p4 <- plot_grid(p3,NULL,p1,p2, nrow = 2, ncol = 2, rel_widths = c(2:1), rel_heights = c(1:2))
  return(list(sc_df_temp,p4))
}
write.table(rbind(sc_df[["ter1_5"]][[1]], sc_df[["ter2_5"]][[1]], sc_df[["ter3_5"]][[1]], sc_df[["ter4_5"]][[1]]), 
            file = 'Figure3_ADARKO_R1/KO_4ter_guidecalling_min5_20230201.txt', sep = ",", quote = F)


### Put guide calling information into Seurat object
load("Figure3_ADARKO_R1/KO_4ter_project_SO_annot_20230201.rda")
KO_gc <- read.table('Figure3_ADARKO_R1/KO_4ter_guidecalling_min5_20230201.txt', sep = ",", header = T)

ko_gid <- KO_gc$g_id
names(ko_gid) <- row.names(KO_gc)
KO <- AddMetaData(KO, metadata = ko_gid, col.name = "guide_id" )
# subset to only cells with an assigned guide ID
length(which(!is.na(KO$guide_id) & KO$guide_id != "multiplet"))
KO <- subset(KO, cells = names(KO$guide_id[!is.na(KO$guide_id) & KO$guide_id != "multiplet"]))
KO

table(KO$guide_id)

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

pdf("Figure3_ADARKO_R1/KO_4ter_project_refUMAP_postGC_20230201.pdf", height = 8, width = 12)
plot_grid(
  DimPlot(KO, group.by = "predicted.id", label = T, repel = T, cols = cluster_col) + NoLegend(),
  DimPlot(KO, group.by = "guide_id", label = T, repel = T) + NoLegend()
)
dev.off()



## collapse cell types similar to wild type
KO$annot_v1 <- KO$predicted.id
#Neural Progenitors: RG + CycProg
KO$annot_v1[KO$predicted.id %in% c("Radial Glia","CycProg")] <- "Neural-Progenitors"
#Neurons: Early Neurons + Retinal Neurons
KO$annot_v1[KO$predicted.id %in% c("Early Neurons","Retinal Neurons")] <- "Neurons"
#SCP: Schwann + Melanobalsts
KO$annot_v1[KO$predicted.id %in% c("Schwann Cells","Melanoblasts")] <- "SCP"

#Gut Epithelial: Foregut Epi + Mid/Hindgut Epi + Airway Epi
KO$annot_v1[KO$predicted.id %in% c("Airway Epi","Foregut Epi","Mid/Hindgut Epi")] <- "Gut-Epi"

#Hematopoietic: HSC + immune + Erythrocyte
KO$annot_v1[KO$predicted.id %in% c("HSC","Immune","Erythrocyte")] <- "Hematopoietic"
#Muscle: Muscle Prog + Cardiac/Skeletal Muscle
KO$annot_v1[KO$predicted.id %in% c("Muscle Prog","Cardiac/Skeletal Muscle")] <- "Muscle"

KO$annot_v1 <- str_replace_all(string = KO$annot_v1, pattern = "/", replacement = "-")
KO$annot_v1 <- str_replace_all(string = KO$annot_v1, pattern = fixed(" "), replacement = "-")
table(KO$annot_v1)

save(KO, file = "Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_20230201.rda")

KO$guide <- sapply(strsplit(KO$guide_id, split = "-|_"),"[[",1)

# ADARB1_1 removed due to inadequate cutting efficiency
KO_sub <- subset(KO, guide_id != "ADARB1_1")
table(KO$guide_id)

save(KO_sub, file = "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")
KO <- KO_sub

### save cell barcode list
for (i in names(table(KO$annot_v1))){
  for (j in names(table(KO$exp))){
    for(k in names(table(KO$guide))){
      if (length (row.names(KO@meta.data[KO$annot_v1 == i & KO$exp == j & KO$guide == k, ])) != 0){
        write.table(paste0("CB:Z:",sapply(strsplit(row.names(KO@meta.data[KO$annot_v1 == i & KO$exp == j & KO$guide == k, ]), split = "_"), "[[", 2)), 
                    paste0('cb_listsKO/v2/',i,'_',j,'_',k,'.txt'), col.names = FALSE, row.names = FALSE, 
                    quote = FALSE, sep = ',')
      }
    }
  }
}

for (j in names(table(KO_sub$exp))){
    for(k in names(table(KO_sub$guide))){
      if (length (row.names(KO_sub@meta.data[KO_sub$exp == j & KO_sub$guide == k, ])) != 0){
        write.table(paste0("CB:Z:",sapply(strsplit(row.names(KO_sub@meta.data[KO_sub$exp == j & KO_sub$guide == k, ]), split = "_"), "[[", 2)), 
                    paste0('cb_listsKO/v7/',j,'_',k,'.txt'), col.names = FALSE, row.names = FALSE, 
                    quote = FALSE, sep = ',')
      }
    }
  }
library(DESeq2)
library(Seurat)
library(ggplot2)
library(ggrepel)

load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")
KO_sub

Idents(KO_sub) <- KO_sub$guide

comps <- data.frame(compA = c("AAVS1","ADAR1","ADAR1","ADARB1","ADARB1","ADARB2","ADARB2"), compB = c("NTC","NTC","AAVS1","NTC","AAVS1","NTC","AAVS1"))
mycolors <- c("blue", "black", "red")
names(mycolors) <- c("DOWN", "NO", "UP")

### repeat the same thing for each germ layer
for (k in names(table(KO_sub$germlayer))){
  so_ct_tmp <- subset(KO_sub, cells = row.names(KO_sub@meta.data[KO_sub$germlayer == k,]))
  for (j in 1:nrow(comps)){
    if (nrow(so_ct_tmp@meta.data[so_ct_tmp$guide == comps[j,"compA"],]) >= 3 & nrow(so_ct_tmp@meta.data[so_ct_tmp$guide == comps[j,"compB"],]) >= 3){
      markers_tmp <- FindMarkers(so_ct_tmp, ident.1 = comps[j,"compA"], ident.2 = comps[j,"compB"], logfc.threshold = 0, test.use = "wilcox",
                                 min.pct = 0.1, min.diff.pct = -Inf, slot = "data", pseudocount.use = 1/10000)
      markers_tmp$p_val_adj_fdr = p.adjust(markers_tmp$p_val, method='fdr')
      markers_tmp[abs(markers_tmp$avg_log2FC) > 0.25 , "differential"] <- "Differential"
      markers_tmp[abs(markers_tmp$avg_log2FC) <= 0.25, "differential"] <- "Not Differential"
      markers_tmp[markers_tmp$p_val_adj_fdr < 0.05, "significance"] <- "Significant"
      markers_tmp[markers_tmp$p_val_adj_fdr >= 0.05, "significance"] <- "Not Significant"
      markers_tmp$average_exp.1 <- AverageExpression(KO_sub, assays = "RNA", slot = "data", features = row.names(markers_tmp), group.by = "guide")$RNA[,comps[j,"compA"]]
      markers_tmp$average_exp.2 <- AverageExpression(KO_sub, assays = "RNA", slot = "data", features = row.names(markers_tmp), group.by = "guide")$RNA[,comps[j,"compB"]]
      markers_tmp$log2FC_from_avExp <- log2(markers_tmp$average_exp.1 / markers_tmp$average_exp.2)
      
      markers_tmp$log2_average_exp.1 <- log2(markers_tmp$average_exp.1)
      markers_tmp$log2_average_exp.2 <- log2(markers_tmp$average_exp.2)
      markers_tmp[markers_tmp$log2_average_exp.1 > 1, "significant_exp.1"] <- "Significant"
      markers_tmp[markers_tmp$log2_average_exp.1 <= 1, "significant_exp.1"] <- "Not Significant"
      markers_tmp[markers_tmp$log2_average_exp.2 > 1, "significant_exp.2"] <- "Significant"
      markers_tmp[markers_tmp$log2_average_exp.2 <= 1, "significant_exp.2"] <- "Not Significant"
      
      markers_tmp$gene_symbol <- row.names(markers_tmp)
      
      write.table(markers_tmp, file = paste0('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/DifferentialGene/SeuratFindMarkers_wilcox/ByGL_minpct10/KO_4ter_diffGene_',k,'_',comps[j,"compA"],'vs',comps[j,"compB"],'_20230505.txt'), quote = F, col.names = T, row.names = T, sep = ",")
    }
  }
}

## for plotting for each germ layer
p <- list()
for (k in names(table(KO_sub$germlayer))){
  for (j in 1:nrow(comps)){
    markers_tmp <- read.table(paste0('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/DifferentialGene/SeuratFindMarkers_wilcox/ByGL/KO_4ter_diffGene_',k,'_',comps[j,"compA"],'vs',comps[j,"compB"],'_20230411.txt'), quote = "", header = T, row.names = 1, sep = ",")
    
    # add a column of NAs, this will be implemented after we decide on some logfc threshold
    markers_tmp$diffexpressed <- "NO"
    # if log2Foldchange > 0.25 and pvalue_fdr < 0.05, set as "UP" 
    markers_tmp$diffexpressed[markers_tmp$avg_log2FC > 0.25 & markers_tmp$p_val_adj_fdr < 0.05] <- "UP"
    # if log2Foldchange < -0.25 and pvalue_fdr < 0.05, set as "DOWN"
    markers_tmp$diffexpressed[markers_tmp$avg_log2FC < -0.25 & markers_tmp$p_val_adj_fdr < 0.05] <- "DOWN"
    markers_tmp$delabel <- NA
    markers_tmp$delabel[markers_tmp$diffexpressed != "NO"] <- markers_tmp$gene_symbol[markers_tmp$diffexpressed != "NO"]
    
    # save plot code for future
    p[[paste0(k, "_", comps[j,"compA"],'vs',comps[j,"compB"])]] <- ggplot(data=markers_tmp, aes(x=avg_log2FC, y=-log10(p_val_adj_fdr), col=diffexpressed, label=delabel)) + scale_color_manual(values=mycolors) +
      geom_point() + geom_vline(xintercept=c(-0.25, 0.25), col="red") +
      geom_hline(yintercept=-log10(0.05), col="red") +
      theme_classic() +
      geom_text_repel() + ggtitle(paste0(k, "_",comps[j,"compA"],'vs',comps[j,"compB"]))
  }
}
pdf('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/DifferentialGene/SeuratFindMarkers_wilcox/KO_diffGene_byGL_20230505.pdf', width = 12, height = 8)
p
dev.off()
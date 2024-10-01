load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")

Idents(KO_sub) <- KO_sub$guide
comps <- data.frame(compA = c("ADAR1","ADARB1","ADARB2"), compB = rep("AAVS1", 3))
JNKgenes <- c("TRAF3","TRAF4","GADD45A","GADD45B","CCN2","GADD45G","TNFRSF19","BECN1")
mycolors <- c("DOWN"="blue", "NO"="grey", "UP"="red")
## for plotting for each germ layer
p <- list()
for (k in names(table(KO_sub$germlayer))){
  for (j in 1:nrow(comps)){
    markers_tmp <- read.table(paste0(input_dir,'/KO_4ter_diffGene_',k,'_',comps[j,"compA"],'vs',comps[j,"compB"],'_20230505.txt'), quote = "", header = T, row.names = 1, sep = ",")
    print(markers_tmp$gene_symbol[markers_tmp$gene_symbol %in% JNKgenes])
    # add a column of NAs, this will be implemented after we decide on some logfc threshold
    markers_tmp$diffexpressed <- "NO"
    # if log2Foldchange > 0.25 and pvalue_fdr < 0.05, set as "UP" 
    markers_tmp$diffexpressed[markers_tmp$avg_log2FC > 0.4 & markers_tmp$p_val_adj_fdr < 0.05] <- "UP"
    # if log2Foldchange < -0.25 and pvalue_fdr < 0.05, set as "DOWN"
    markers_tmp$diffexpressed[markers_tmp$avg_log2FC < -0.4 & markers_tmp$p_val_adj_fdr < 0.05] <- "DOWN"
    markers_tmp$delabel <- NA
    markers_tmp$delabel[markers_tmp$gene_symbol %in% JNKgenes] <- markers_tmp$gene_symbol[markers_tmp$gene_symbol %in% JNKgenes]
    #markers_tmp$delabel[markers_tmp$diffexpressed != "NO"] <- markers_tmp$gene_symbol[markers_tmp$diffexpressed != "NO"]
    
    highlight_df <- markers_tmp[markers_tmp$gene_symbol %in% JNKgenes,]
    # save plot code for future
    p[[paste0(k, "_", comps[j,"compA"],'vs',comps[j,"compB"])]] <- ggplot(data=markers_tmp, aes(x=avg_log2FC, y=-log10(p_val_adj_fdr), col=diffexpressed, label=delabel)) + scale_color_manual(values=mycolors) +
      geom_point(size = 0.5) + geom_vline(xintercept=c(-0.4, 0.4), col="red") +
      geom_point(data=highlight_df, aes(x=avg_log2FC,y=-log10(p_val_adj_fdr)), color='black', size=2) + xlim(-2,2) + ylim(0,5) +
      geom_hline(yintercept=-log10(0.05), col="red") +
      theme_classic() +
      geom_text_repel(colour = "black") + ggtitle(paste0(k, "_",comps[j,"compA"],'vs',comps[j,"compB"]))
  }
}
pdf('Figure3_ADARKO_R1/DiffGene_VolcanoPlot_test.pdf', width = 12, height = 8)
p
dev.off()

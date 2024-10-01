library(Seurat)
library(ggrepel)

load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")
KO_sub

comps <- data.frame(compA = c("ADAR1","AAVS1"), compB = c("AAVS1","NTC"))

adipgenes <- c("FABP7","SULT1E1","DLK1","FRZB","KLF4","KLF6","CD24")
stressgenes <- c("JUN","SOD3","TRAF4","GADD45B","CCN2","XBP1")
glycolysisgenes <- c("ENO1","LDHA","IRX3","IRX5")

combgenes <- c(adipgenes,stressgenes,glycolysisgenes)

mycolors <- c("DOWN"="blue", "NO"="grey", "UP"="red")

so_ct_tmp <- subset(KO_sub, cells = row.names(KO_sub@meta.data[KO_sub$annot_v1 == "Adipogenic-MSC-Fib",]))
Idents(so_ct_tmp) <- so_ct_tmp$guide
markers_tmp <- FindMarkers(so_ct_tmp, ident.1 = "ADAR1", ident.2 = "AAVS1", logfc.threshold = 0, test.use = "wilcox", assay = "RNA",
                                 min.pct = 0.1, min.diff.pct = -Inf, slot = "data", pseudocount.use = 1/10000)
markers_tmp$p_val_adj_fdr = p.adjust(markers_tmp$p_val, method='fdr')
markers_tmp$diffexpressed <- "NO"
# if log2Foldchange > 0.25 and pvalue_fdr < 0.05, set as "UP" 
markers_tmp$diffexpressed[markers_tmp$avg_log2FC > 0.25 & markers_tmp$p_val_adj_fdr < 0.05] <- "UP"
# if log2Foldchange < -0.25 and pvalue_fdr < 0.05, set as "DOWN"
markers_tmp$diffexpressed[markers_tmp$avg_log2FC < -0.25 & markers_tmp$p_val_adj_fdr < 0.05] <- "DOWN"
markers_tmp$delabel <- NA
markers_tmp$delabel[markers_tmp$gene_symbol %in% combgenes] <- markers_tmp$gene_symbol[markers_tmp$gene_symbol %in% combgenes]
#markers_tmp$delabel[markers_tmp$diffexpressed != "NO"] <- markers_tmp$gene_symbol[markers_tmp$diffexpressed != "NO"]

markers_tmp$gene_symbol <- row.names(markers_tmp)
highlight_df <- markers_tmp[markers_tmp$gene_symbol %in% combgenes,]
# save plot code for future
pdf('Figure3_ADARKO_R1/DiffGene_VolcanoPlot_adipmscfib_ADAR1vsAAVS1_20230814.pdf', width = 12, height = 8)
ggplot(data=markers_tmp, aes(x=avg_log2FC, y=-log10(p_val_adj_fdr), col=diffexpressed, label=delabel)) + scale_color_manual(values=mycolors) +
      geom_point(size = 0.5) + geom_vline(xintercept=c(-0.25, 0.25), col="red") +
      geom_point(data=highlight_df, aes(x=avg_log2FC,y=-log10(p_val_adj_fdr)), color='black', size=2) + xlim(-2,2) + ylim(0,5) +
      geom_hline(yintercept=-log10(0.05), col="red") +
      theme_classic() +
      geom_text_repel(data = highlight_df, aes(label = gene_symbol), colour = "black") + ggtitle("ADAR1 vs. AAVS1 for Adipogenic MSC Fib")
dev.off()

write.table(markers_tmp, "Figure3_ADARKO_R1/DiffGene_VolcanoPlot_adipmscfib_ADAR1vsAAVS1_20230814.txt", quote = F, col.names = T, sep = ",")

so_ct_tmp <- subset(KO_sub, cells = row.names(KO_sub@meta.data[KO_sub$annot_v1 == "Adipogenic-MSC-Fib",]))
Idents(so_ct_tmp) <- so_ct_tmp$guide
markers_tmp <- FindMarkers(so_ct_tmp, ident.1 = "AAVS1", ident.2 = "NTC", logfc.threshold = 0, test.use = "wilcox", assay = "RNA",
                                 min.pct = 0.1, min.diff.pct = -Inf, slot = "data", pseudocount.use = 1/10000)
markers_tmp$p_val_adj_fdr = p.adjust(markers_tmp$p_val, method='fdr')
markers_tmp$diffexpressed <- "NO"
# if log2Foldchange > 0.25 and pvalue_fdr < 0.05, set as "UP" 
markers_tmp$diffexpressed[markers_tmp$avg_log2FC > 0.25 & markers_tmp$p_val_adj_fdr < 0.05] <- "UP"
# if log2Foldchange < -0.25 and pvalue_fdr < 0.05, set as "DOWN"
markers_tmp$diffexpressed[markers_tmp$avg_log2FC < -0.25 & markers_tmp$p_val_adj_fdr < 0.05] <- "DOWN"
markers_tmp$delabel <- NA
markers_tmp$delabel[markers_tmp$gene_symbol %in% combgenes] <- markers_tmp$gene_symbol[markers_tmp$gene_symbol %in% combgenes]
#markers_tmp$delabel[markers_tmp$diffexpressed != "NO"] <- markers_tmp$gene_symbol[markers_tmp$diffexpressed != "NO"]

markers_tmp$gene_symbol <- row.names(markers_tmp)
highlight_df <- markers_tmp[markers_tmp$gene_symbol %in% combgenes,]
# save plot code for future
pdf('Figure3_ADARKO_R1/DiffGene_VolcanoPlot_adipmscfib_AAVS1vsNTC_20230814.pdf', width = 12, height = 8)
ggplot(data=markers_tmp, aes(x=avg_log2FC, y=-log10(p_val_adj_fdr), col=diffexpressed, label=delabel)) + scale_color_manual(values=mycolors) +
      geom_point(size = 0.5) + geom_vline(xintercept=c(-0.25, 0.25), col="red") +
      geom_point(data=highlight_df, aes(x=avg_log2FC,y=-log10(p_val_adj_fdr)), color='black', size=2) + xlim(-2,2) + ylim(0,5) +
      geom_hline(yintercept=-log10(0.05), col="red") +
      theme_classic() +
      geom_text_repel(data = highlight_df, aes(label = gene_symbol), colour = "black") + ggtitle("AAVS1 vs. NTC for Adipogenic MSC Fib")
dev.off()

write.table(markers_tmp, "Figure3_ADARKO_R1/DiffGene_VolcanoPlot_adipmscfib_AAVS1vsNTC_20230814.txt", quote = F, col.names = T, sep = ",")

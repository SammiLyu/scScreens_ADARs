library(Seurat)
library(swne)

load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")
KO_sub

adipgenes <- c("AXIN1","AXIN2","CEBPA","CEBPB","CEBPG","FABP7","KLF4","KLF6","SULT1E1","DLK1","FRZB")
stressgenes <- c("CASP6","CCN2","GADD45A","GADD45B","GADD45G","JUN","TRAF3","TRAF4","SOD3")
glycolysisgenes <- c("IRX3","IRX5","ENO1","ENO2","GAPDH","HK1","HK2","LDHA")

KO_sub$annot_guide <- paste0(KO_sub$annot_v1,"_",KO_sub$guide)
so_ct_tmp <- subset(KO_sub, cells = row.names(KO_sub@meta.data[KO_sub$annot_v1 == "Adipogenic-MSC-Fib",]))
so_ct_tmp <- ScaleData(so_ct_tmp, features = row.names(so_ct_tmp))
#combgenes[which(!combgenes %in% row.names(GetAssayData(so_ct_tmp, assay = "RNA", slot = "scale.data")))]

heat.mat <- t(apply(GetAssayData(so_ct_tmp, assay = "RNA", slot = "scale.data")[rev(c(adipgenes, stressgenes, glycolysisgenes)),], 1, function(x) {
  tapply(x, so_ct_tmp$annot_guide, mean)
}))
heat.mat[heat.mat > 2.5] <- 2.5

heat.mat.sub <- heat.mat[,sapply(strsplit(colnames(heat.mat), split = "_"),"[[",1) == "Adipogenic-MSC-Fib"]

#pdf("Figure3_ADARKO_R1/DiffGene_Heatmap_adipmscfib_allprtb_20230814.pdf", width = 6, height = 7)
#swne:::ggHeat(heat.mat.sub, clustering = "none", x.lab.size = 11, y.lab.size = 11) + ggtitle('Selected Gene Expression for Adipogenic MSC Fib')
#dev.off()

heat.mat.sub <- heat.mat[,sapply(strsplit(colnames(heat.mat), split = "_"),"[[",2) %in% c("ADAR1","AAVS1","NTC")]

pdf("Figure3_ADARKO_R1/DiffGene_Heatmap_adipmscfib_adar1only_glycolysisgenes_mergingexp_20230815.pdf", width = 3, height = 6)
swne:::ggHeat(heat.mat.sub, clustering = "none", x.lab.size = 11, y.lab.size = 11) + ggtitle('Selected Gene Expression for Adipogenic MSC Fib')
dev.off()

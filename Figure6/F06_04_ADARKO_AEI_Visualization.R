library(stats)

aei_table_df <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_AEI_byCT_byExp_20230510.txt', header = T, sep =',')
guide_order <- c("NTC","AAVS1","ADAR1","ADARB1","ADARB2")

# quick plot
df_aei_by_prtb <- data.frame(matrix(NA, nrow = nrow(aei_table_df), ncol = 3, dimnames = list(row.names(aei_table_df), c("AEI_A2G","cell_type","prtb"))))
df_aei_by_prtb$AEI_A2G <- aei_table_df$A2GEditingIndex
df_aei_by_prtb$cell_type <- aei_table_df$CellType
df_aei_by_prtb$prtb <- aei_table_df$Guide
df_aei_by_prtb$ter <- sapply(strsplit(row.names(df_aei_by_prtb), split = "_"), "[[", 2)

cluster_col <- c("#0000FF","#FF0000","#00C000","#F07CAB","#000080",
                 "#AD07E3","#FF8000","#000000","#90BFF9","#C0C0FF",
                 "#D30B94","#A00000","#F2B77C","#00FF00","#94641F")
names(cluster_col) <- names(table(aei_table_df$CellType))

pdf("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_AEI_20230204.pdf", width = 18, height = 8)
plot_grid(
  ggplot(df_aei_by_prtb, aes(x = factor(prtb, levels = guide_order), y = AEI_A2G)) + geom_jitter(size=3, aes(color = cell_type)) + scale_color_manual(values = cluster_col) + 
    geom_smooth(method='lm', fullrange=TRUE) +
    #geom_text(x = 5, y = 0.65, label = eq(log2(df_aei_by_prtb$ADAR_TPM+0.01)+log2(df_aei_by_prtb$ADARB1_TPM+0.01),df_aei_by_prtb$AEI_A2G), parse = TRUE) +
    theme_classic(),
  
  ggplot(df_aei_by_prtb, aes(x = factor(prtb, levels = guide_order), y = AEI_A2G)) + geom_jitter(size=3, aes(color = factor(ter, levels = c("ter1","ter2","ter3","ter4")))) + geom_smooth(method='lm', fullrange=TRUE) +
    #geom_text(x = 5, y = 0.65, label = eq(log2(df_aei_by_prtb$ADAR_TPM+0.01)+log2(df_aei_by_prtb$ADARB1_TPM+0.01),df_aei_by_prtb$AEI_A2G), parse = TRUE) +
    theme_classic()
)
dev.off()

### anova for ADARKO AEI

res.aov <- aov(A2GEditingIndex ~ Guide, data = aei_table_df)
summary(res.aov)
TukeyHSD(res.aov)

# script for running cohens'd analysis on AEI and adar exp for organ development + visulization

aei_table_df <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_OrganDev_v2/output_index_20230607/OrganDev_v1_outputindex/EditingIndex.csv', header = T, sep =',')
adar_tpm <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/fc_OD_v1/adar_exp_tpm.txt', header = T, sep ='\t')

# formulate "meta data"
metadata_df <- data.frame(matrix(NA, nrow = dim(aei_table_df)[1], ncol = 7))
row.names(metadata_df) <- aei_table_df$Sample
colnames(metadata_df) <- c("sex","age","organ","ADAR1exp","ADAR2exp","ADAR3exp","AEI_A2GEditing")
metadata_df$sex <- sapply(strsplit(row.names(metadata_df), split = "[.]"), "[[", 5)
metadata_df$age <- sapply(strsplit(row.names(metadata_df), split = "[.]"), "[[", 4)
metadata_df$organ <- sapply(strsplit(row.names(metadata_df), split = "[.]"), "[[", 3)

OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable.csv', sep = ',', header = T)
metadata_df$library.ID <- sapply(strsplit(row.names(metadata_df), split = 'sTS'), "[[", 1)
OrganDev_meta <- OrganDev_meta[OrganDev_meta$library.ID %in% metadata_df$library.ID, ]
metadata_df <- metadata_df[match(OrganDev_meta$library.ID, metadata_df$library.ID),]

metadata_df$DevelopmentStage <- OrganDev_meta[['Developmental.stage']]
metadata_df$Organ <- OrganDev_meta[['Organ']]

row.names(aei_table_df) <- aei_table_df$Sample

adar_tpm$Sample <- sapply(strsplit(adar_tpm$Sample, split = "v1_nomultiseg_"), "[[",2)
adar_tpm$Sample <- sapply(strsplit(adar_tpm$Sample, split = "_nomultiseg_bam"), "[[",1)
adar_tpm$Sample <- stringr:::str_replace_all(adar_tpm$Sample, "_", ".")

for (i in 1:nrow(metadata_df)){
    if (length(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADAR", "TPM"]) != 0 &
        adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADAR", "TPM"] != 0){
        metadata_df$ADAR1exp[i] <- log2(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADAR", "TPM"])
    }
    if (length(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB1", "TPM"]) != 0 &
        adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB1", "TPM"] != 0){
        metadata_df$ADAR2exp[i] <- log2(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB1", "TPM"])
    }
    if (length(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB2", "TPM"]) != 0 &
    adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB2", "TPM"] != 0){
        metadata_df$ADAR3exp[i] <- log2(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB2", "TPM"])
    }

    if (length(aei_table_df[row.names(metadata_df)[i], "A2GEditingIndex"]) != 0 ){
        metadata_df$AEI_A2GEditing[i] <- aei_table_df[row.names(metadata_df)[i], "A2GEditingIndex"]
    }
}
dim(metadata_df)

metadata_df$sex <- factor(metadata_df$sex, levels = c('Female','Male'))
metadata_df$Organ <- factor(metadata_df$Organ, levels = c('Forebrain/Cerebrum','Hindbrain/Cerebellum','Heart','Liver','Kidney','Ovary','Testis'))
metadata_df$DevelopmentStage <- factor(metadata_df$DevelopmentStage, levels = c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc','11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc','newborn','infant','toddler','school','youngTeenager','teenager','oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'))

# ensure continuous x scale
lat.order.levels <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc','11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc','newborn','infant','toddler','school','youngTeenager','teenager','oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')
DevelopmentStage.data <- data.frame(
  DevelopmentStage = lat.order.levels,
  DevelopmentStage_num = 1:length(lat.order.levels)
)
df <- dplyr::left_join(metadata_df, DevelopmentStage.data, by= "DevelopmentStage")

min_exp <- floor(min(c(df$ADAR1exp, df$ADAR2exp)))
max_exp <- ceiling(max(c(df$ADAR1exp, df$ADAR2exp)))

min_aei <- floor(min(df$AEI_A2GEditing))
max_aei <- ceiling(max(df$AEI_A2GEditing))

coeff <- (max_aei - min_aei) / (max_exp - min_exp) ; coeff

df <- tidyr::complete(df, DevelopmentStage_num, organ)

p <- list()
for (i in names(table(df$Organ))){
    if (dim(df[df$Organ == i,])[2] != 0){
    
    metadata_df_sub <- reshape2:::melt(df[df$Organ == i,c("ADAR1exp","ADAR2exp","AEI_A2GEditing","DevelopmentStage_num")], id.vars="DevelopmentStage_num")
    metadata_df_sub[metadata_df_sub$variable == "AEI_A2GEditing","value"] <- metadata_df_sub[metadata_df_sub$variable == "AEI_A2GEditing","value"] * coeff + min_exp

    p[[i]] <- ggplot(metadata_df_sub, aes(x = DevelopmentStage_num, y = value)) + geom_point(size=7.5, aes(color = variable, shape = variable)) + scale_color_manual(values = c('#619CFF','#00BA38','#F8766D')) + scale_shape_manual(values = c(21,21,16)) +
    geom_smooth(aes(color = variable, linetype = variable), linewidth=1.5, se = F) +  scale_linetype_manual(values = c(2,2,1)) +
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + scale_y_continuous(
    # Features of the first axis
    name = "Gene Expression",
    limits = c(min_exp, max_exp),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.-3, name="AEI")
  ) + theme_classic() + ggtitle(i) + coord_cartesian(xlim = c(1,25)) + 
theme(axis.text=element_text(size=30,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)
    }
}
#p[[1]]
pdf('Figure_OrganDevelopment_R1/OD_allOrgan_ADAR1ADARB1AEI_splitbyOrgan_20230821.pdf', height = 10, width = 22)
p
dev.off()

## ADAR3 expression not plotting but code is here
#min_adar3exp <- floor(min(df$ADAR3exp, na.rm = T))
#max_adar3exp <- ceiling(max(df$ADAR3exp, na.rm = T))

#p <- list()
#for (i in names(table(df$Organ))){
#    if (dim(df[df$Organ == i,])[2] != 0){
#     metadata_df_sub <- df[df$Organ == i,]

#     p[[i]] <- ggplot(metadata_df_sub, aes(x = DevelopmentStage_num, y = ADAR3exp)) + geom_point(size=7.5, aes(color = sex)) + geom_smooth(color = "black", linewidth=1.5) +
#   scale_x_continuous(
#     breaks=1:length(lat.order.levels),
#     labels=lat.order.levels
#   ) + theme_classic() + ggtitle(i) + geom_hline(yintercept = log2(5), col = "red", linetype="dotdash") + 
# theme(axis.text=element_text(size=30,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"),
# legend.title = element_text(size=18), #change legend title font size
#         legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") + ylim(min_adar3exp, max_adar3exp)
#     }
# }
# pdf('Figure_OrganDevelopment_R1/OD_allOrgan_ADARB2_splitbyOrgan_20230814.pdf', height = 10, width = 22)
# p
# dev.off()

#metadata_df <- metadata_df[complete.cases(metadata_df),]
#dim(metadata_df)

dim(metadata_df)
for (i in row.names(metadata_df)){
    if (!is.na(metadata_df[i,"ADAR3exp"]) & metadata_df[i,"ADAR3exp"] < log2(5)){
        metadata_df[i, "ADAR3exp"] <- NA
    }
}
head(metadata_df)

table(metadata_df$DevelopmentStage)
table(metadata_df$Organ)

metadata_df_subembryo <- metadata_df[metadata_df$DevelopmentStage %in% c("4wpc","5wpc","6wpc","7wpc","8wpc","9wpc","10wpc","11wpc","12wpc","13wpc","16wpc","18wpc","19wpc","20wpc"),]
metadata_df_subembryo

earlier_time <- c("4wpc","5wpc","6wpc","7wpc","8wpc","9wpc","10wpc")
later_time <- c("11wpc","12wpc","13wpc","16wpc","18wpc","19wpc","20wpc")

cohensd_res <- data.frame(matrix(NA, nrow = length(table(metadata_df_subembryo$Organ)), ncol = 16))
row.names(cohensd_res) <- names(table(metadata_df_subembryo$Organ))
colnames(cohensd_res) <- paste0(rep(c("AEI","ADAR1","ADAR2","ADAR3"), each = 4), "_", rep(c("d_estimate","lower_ci","upper_ci","pvalue")))

for (ct in names(table(metadata_df_subembryo$Organ))){
    if (dim(metadata_df_subembryo[metadata_df_subembryo$Organ == ct,])[1] !=0){
        print(ct)

        group_early <- metadata_df_subembryo[metadata_df_subembryo$Organ == ct & metadata_df_subembryo$DevelopmentStage %in% earlier_time,  ]
        group_late <- metadata_df_subembryo[metadata_df_subembryo$Organ == ct & metadata_df_subembryo$DevelopmentStage %in% later_time,  ]
        
        cohensd_res[ct, "AEI_d_estimate"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "AEI_lower_ci"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "AEI_upper_ci"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR1_d_estimate"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR1_lower_ci"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR1_upper_ci"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR2_d_estimate"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR2_lower_ci"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR2_upper_ci"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR3_d_estimate"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR3_lower_ci"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR3_upper_ci"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$conf.int["upper"]

        if (length(group_early$AEI_A2GEditing[!is.na(group_early$AEI_A2GEditing)]) > 1 & length(group_late$AEI_A2GEditing[!is.na(group_late$AEI_A2GEditing)]) > 1 ) {
            cohensd_res[ct, "AEI_pvalue"] <- t.test(group_early$AEI_A2GEditing, group_late$AEI_A2GEditing, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR1exp[!is.na(group_early$ADAR1exp)]) > 1 & length(group_late$ADAR1exp[!is.na(group_late$ADAR1exp)]) > 1 ) {
            cohensd_res[ct, "ADAR1_pvalue"] <- t.test(group_early$ADAR1exp, group_late$ADAR1exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR2exp[!is.na(group_early$ADAR2exp)]) > 1 & length(group_late$ADAR2exp[!is.na(group_late$ADAR2exp)]) > 1 ) {
            cohensd_res[ct, "ADAR2_pvalue"] <- t.test(group_early$ADAR2exp, group_late$ADAR2exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR3exp[!is.na(group_early$ADAR3exp)]) > 1 & length(group_late$ADAR3exp[!is.na(group_late$ADAR3exp)]) > 1 ) {
            cohensd_res[ct, "ADAR3_pvalue"] <- t.test(group_early$ADAR3exp, group_late$ADAR3exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
    }
}

write.table(cohensd_res, "Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_4wpc_20wpc_202300815.txt", sep = ",", row.names = T, quote = F)

csd_organ_order <- c("Forebrain/Cerebrum","Hindbrain/Cerebellum","Heart","Kidney","Liver","Testis")
organ_fill_col <- c("Forebrain/Cerebrum"="#990000","Hindbrain/Cerebellum"="#FF0000","Heart"="#FF9933",
"Kidney"="#66FF33","Liver"="#66FFFF","Testis"="#9933FF")
cohensd_res <- cohensd_res[row.names(cohensd_res) != "Ovary",]
bound_limit <- ceiling(max(abs(cohensd_res[!is.na(cohensd_res)])))

pdf("Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_4wpc_20wpc_20230815.pdf", height = 12, width = 25)
cowplot:::plot_grid(
ggplot(cohensd_res[,c(1:4)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = AEI_lower_ci, 
    upper = AEI_upper_ci, 
    ymin = AEI_lower_ci, 
    ymax = AEI_upper_ci, 
    middle = AEI_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit,bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('AEI'),
ggplot(cohensd_res[,c(5:8)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR1_lower_ci, 
    upper = ADAR1_upper_ci, 
    ymin = ADAR1_lower_ci, 
    ymax = ADAR1_upper_ci, 
    middle = ADAR1_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit,bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADAR1 Expression in TPM'),
ggplot(cohensd_res[,c(9:12)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR2_lower_ci, 
    upper = ADAR2_upper_ci, 
    ymin = ADAR2_lower_ci, 
    ymax = ADAR2_upper_ci, 
    middle = ADAR2_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit,bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADARB1 Expression in TPM'),
ggplot(cohensd_res[,c(13:16)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR3_lower_ci, 
    upper = ADAR3_upper_ci, 
    ymin = ADAR3_lower_ci, 
    ymax = ADAR3_upper_ci, 
    middle = ADAR3_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit,bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADARB2 Expression in TPM'),
ncol = 4
)
dev.off()

### second comp

metadata_df_subbirth <- metadata_df[metadata_df$DevelopmentStage %in% c("11wpc","12wpc","13wpc","16wpc","18wpc","19wpc","20wpc","newborn","infant","toddler","school","teenager","youngTeenager","oldTeenager"),]
metadata_df_subbirth

earlier_time <- c("11wpc","12wpc","13wpc","16wpc","18wpc","19wpc","20wpc")
later_time <- c("newborn","infant","toddler","school","teenager","youngTeenager","oldTeenager")

cohensd_res <- data.frame(matrix(NA, nrow = length(table(metadata_df_subbirth$Organ)), ncol = 16))
row.names(cohensd_res) <- names(table(metadata_df_subbirth$Organ))
colnames(cohensd_res) <- paste0(rep(c("AEI","ADAR1","ADAR2","ADAR3"), each = 4), "_", rep(c("d_estimate","lower_ci","upper_ci","pvalue")))

for (ct in names(table(metadata_df_subbirth$Organ))){
    if (dim(metadata_df_subbirth[metadata_df_subbirth$Organ == ct,])[1] !=0){
        print(ct)

        group_early <- metadata_df_subbirth[metadata_df_subbirth$Organ == ct & metadata_df_subbirth$DevelopmentStage %in% earlier_time,  ]
        group_late <- metadata_df_subbirth[metadata_df_subbirth$Organ == ct & metadata_df_subbirth$DevelopmentStage %in% later_time,  ]
        
        cohensd_res[ct, "AEI_d_estimate"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "AEI_lower_ci"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "AEI_upper_ci"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR1_d_estimate"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR1_lower_ci"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR1_upper_ci"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR2_d_estimate"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR2_lower_ci"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR2_upper_ci"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR3_d_estimate"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR3_lower_ci"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR3_upper_ci"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$conf.int["upper"]

        if (length(group_early$AEI_A2GEditing[!is.na(group_early$AEI_A2GEditing)]) > 1 & length(group_late$AEI_A2GEditing[!is.na(group_late$AEI_A2GEditing)]) > 1 ) {
            cohensd_res[ct, "AEI_pvalue"] <- t.test(group_early$AEI_A2GEditing, group_late$AEI_A2GEditing, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR1exp[!is.na(group_early$ADAR1exp)]) > 1 & length(group_late$ADAR1exp[!is.na(group_late$ADAR1exp)]) > 1 ) {
            cohensd_res[ct, "ADAR1_pvalue"] <- t.test(group_early$ADAR1exp, group_late$ADAR1exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR2exp[!is.na(group_early$ADAR2exp)]) > 1 & length(group_late$ADAR2exp[!is.na(group_late$ADAR2exp)]) > 1 ) {
            cohensd_res[ct, "ADAR2_pvalue"] <- t.test(group_early$ADAR2exp, group_late$ADAR2exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR3exp[!is.na(group_early$ADAR3exp)]) > 1 & length(group_late$ADAR3exp[!is.na(group_late$ADAR3exp)]) > 1 ) {
            cohensd_res[ct, "ADAR3_pvalue"] <- t.test(group_early$ADAR3exp, group_late$ADAR3exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
    }
}

write.table(cohensd_res, "Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_20wpc_olderteenager_20230815.txt", sep = ",", row.names = T, quote = F)

csd_organ_order <- c("Forebrain/Cerebrum","Hindbrain/Cerebellum","Heart","Kidney","Liver","Testis")
organ_fill_col <- c("Forebrain/Cerebrum"="#990000","Hindbrain/Cerebellum"="#FF0000","Heart"="#FF9933",
"Kidney"="#66FF33","Liver"="#66FFFF","Testis"="#9933FF")
cohensd_res <- cohensd_res[row.names(cohensd_res) != "Ovary",]
bound_limit <- ceiling(max(abs(cohensd_res[!is.na(cohensd_res)])))

pdf("Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_20wpc_olderteenager_20230815.pdf", height = 12, width = 25)
cowplot:::plot_grid(
ggplot(cohensd_res[,c(1:4)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = AEI_lower_ci, 
    upper = AEI_upper_ci, 
    ymin = AEI_lower_ci, 
    ymax = AEI_upper_ci, 
    middle = AEI_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit, bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('AEI'),
ggplot(cohensd_res[,c(5:8)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR1_lower_ci, 
    upper = ADAR1_upper_ci, 
    ymin = ADAR1_lower_ci, 
    ymax = ADAR1_upper_ci, 
    middle = ADAR1_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit, bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADAR1 Expression in TPM'),
ggplot(cohensd_res[,c(9:12)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR2_lower_ci, 
    upper = ADAR2_upper_ci, 
    ymin = ADAR2_lower_ci, 
    ymax = ADAR2_upper_ci, 
    middle = ADAR2_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit, bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADARB1 Expression in TPM'),
ggplot(cohensd_res[,c(13:16)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR3_lower_ci, 
    upper = ADAR3_upper_ci, 
    ymin = ADAR3_lower_ci, 
    ymax = ADAR3_upper_ci, 
    middle = ADAR3_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit, bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADARB2 Expression in TPM'),
ncol = 4
)
dev.off()

### third comp

metadata_df_afterbirth <- metadata_df[metadata_df$DevelopmentStage %in% c("newborn","infant","toddler","school","teenager","youngTeenager","oldTeenager","youngAdult","youngMidAge","olderMidAge","senior"),]
metadata_df_afterbirth

earlier_time <- c("newborn","infant","toddler","school","teenager","youngTeenager")
later_time <- c("oldTeenager","youngAdult","youngMidAge","olderMidAge","senior")

cohensd_res <- data.frame(matrix(NA, nrow = length(table(metadata_df_afterbirth$Organ)), ncol = 16))
row.names(cohensd_res) <- names(table(metadata_df_afterbirth$Organ))
colnames(cohensd_res) <- paste0(rep(c("AEI","ADAR1","ADAR2","ADAR3"), each = 4), "_", rep(c("d_estimate","lower_ci","upper_ci","pvalue")))

for (ct in names(table(metadata_df_afterbirth$Organ))){
    if (dim(metadata_df_afterbirth[metadata_df_afterbirth$Organ == ct,])[1] !=0){
        print(ct)

        group_early <- metadata_df_afterbirth[metadata_df_afterbirth$Organ == ct & metadata_df_afterbirth$DevelopmentStage %in% earlier_time,  ]
        group_late <- metadata_df_afterbirth[metadata_df_afterbirth$Organ == ct & metadata_df_afterbirth$DevelopmentStage %in% later_time,  ]
        
        cohensd_res[ct, "AEI_d_estimate"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "AEI_lower_ci"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "AEI_upper_ci"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR1_d_estimate"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR1_lower_ci"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR1_upper_ci"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR2_d_estimate"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR2_lower_ci"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR2_upper_ci"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR3_d_estimate"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR3_lower_ci"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR3_upper_ci"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$conf.int["upper"]

        if (length(group_early$AEI_A2GEditing[!is.na(group_early$AEI_A2GEditing)]) > 1 & length(group_late$AEI_A2GEditing[!is.na(group_late$AEI_A2GEditing)]) > 1 ) {
            cohensd_res[ct, "AEI_pvalue"] <- t.test(group_early$AEI_A2GEditing, group_late$AEI_A2GEditing, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR1exp[!is.na(group_early$ADAR1exp)]) > 1 & length(group_late$ADAR1exp[!is.na(group_late$ADAR1exp)]) > 1 ) {
            cohensd_res[ct, "ADAR1_pvalue"] <- t.test(group_early$ADAR1exp, group_late$ADAR1exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR2exp[!is.na(group_early$ADAR2exp)]) > 1 & length(group_late$ADAR2exp[!is.na(group_late$ADAR2exp)]) > 1 ) {
            cohensd_res[ct, "ADAR2_pvalue"] <- t.test(group_early$ADAR2exp, group_late$ADAR2exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR3exp[!is.na(group_early$ADAR3exp)]) > 1 & length(group_late$ADAR3exp[!is.na(group_late$ADAR3exp)]) > 1 ) {
            cohensd_res[ct, "ADAR3_pvalue"] <- t.test(group_early$ADAR3exp, group_late$ADAR3exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
    }
}

write.table(cohensd_res, "Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_newborn_senior_20230820.txt", sep = ",", row.names = T, quote = F)
cohensd_res <- read.table("Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_newborn_senior_20230820.txt", sep = ",", header = T, quote = "")

csd_organ_order <- c("Forebrain/Cerebrum","Hindbrain/Cerebellum","Heart","Liver","Testis")
organ_fill_col <- c("Forebrain/Cerebrum"="#990000","Hindbrain/Cerebellum"="#FF0000","Heart"="#FF9933",
"Liver"="#66FFFF","Testis"="#9933FF")
cohensd_res <- cohensd_res[!row.names(cohensd_res) %in% c("Ovary","Kidney"),]
bound_limit <- ceiling(max(abs(cohensd_res[!is.na(cohensd_res)])))

pdf("Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_newborn_senior_20230820.pdf", height = 12, width = 25)
cowplot:::plot_grid(
ggplot(cohensd_res[,c(1:4)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = AEI_lower_ci, 
    upper = AEI_upper_ci, 
    ymin = AEI_lower_ci, 
    ymax = AEI_upper_ci, 
    middle = AEI_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit, bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('AEI'),
ggplot(cohensd_res[,c(5:8)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR1_lower_ci, 
    upper = ADAR1_upper_ci, 
    ymin = ADAR1_lower_ci, 
    ymax = ADAR1_upper_ci, 
    middle = ADAR1_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit, bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADAR1 Expression in TPM'),
ggplot(cohensd_res[,c(9:12)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR2_lower_ci, 
    upper = ADAR2_upper_ci, 
    ymin = ADAR2_lower_ci, 
    ymax = ADAR2_upper_ci, 
    middle = ADAR2_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit, bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADARB1 Expression in TPM'),
ggplot(cohensd_res[,c(13:16)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR3_lower_ci, 
    upper = ADAR3_upper_ci, 
    ymin = ADAR3_lower_ci, 
    ymax = ADAR3_upper_ci, 
    middle = ADAR3_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit, bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADARB2 Expression in TPM'),
ncol = 4
)
dev.off()


### wanted to check to see if for heart, liver and testis, ADAR1 might drive a drop in AEI in early gestation
metadata_df_hlt <- metadata_df[metadata_df$DevelopmentStage %in% c("4wpc","5wpc","6wpc","7wpc","8wpc","9wpc","10wpc","11wpc"),]
metadata_df_hlt

earlier_time <- c("4wpc","5wpc","6wpc","7wpc")
later_time <- c("8wpc","9wpc","10wpc","11wpc")

cohensd_res <- data.frame(matrix(NA, nrow = length(table(metadata_df_hlt$Organ)), ncol = 16))
row.names(cohensd_res) <- names(table(metadata_df_hlt$Organ))
colnames(cohensd_res) <- paste0(rep(c("AEI","ADAR1","ADAR2","ADAR3"), each = 4), "_", rep(c("d_estimate","lower_ci","upper_ci","pvalue")))

for (ct in names(table(metadata_df_hlt$Organ))){
    if (dim(metadata_df_hlt[metadata_df_hlt$Organ == ct,])[1] !=0){
        print(ct)

        group_early <- metadata_df_hlt[metadata_df_hlt$Organ == ct & metadata_df_hlt$DevelopmentStage %in% earlier_time,  ]
        group_late <- metadata_df_hlt[metadata_df_hlt$Organ == ct & metadata_df_hlt$DevelopmentStage %in% later_time,  ]
        
        cohensd_res[ct, "AEI_d_estimate"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "AEI_lower_ci"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "AEI_upper_ci"] <- effsize:::cohen.d(group_late$AEI_A2GEditing,group_early$AEI_A2GEditing, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR1_d_estimate"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR1_lower_ci"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR1_upper_ci"] <- effsize:::cohen.d(group_late$ADAR1exp,group_early$ADAR1exp, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR2_d_estimate"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR2_lower_ci"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR2_upper_ci"] <- effsize:::cohen.d(group_late$ADAR2exp,group_early$ADAR2exp, pooled = T, na.rm = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR3_d_estimate"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$estimate
        cohensd_res[ct, "ADAR3_lower_ci"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR3_upper_ci"] <- effsize:::cohen.d(group_late$ADAR3exp,group_early$ADAR3exp, pooled = T, na.rm = T)$conf.int["upper"]

        if (length(group_early$AEI_A2GEditing[!is.na(group_early$AEI_A2GEditing)]) > 1 & length(group_late$AEI_A2GEditing[!is.na(group_late$AEI_A2GEditing)]) > 1 ) {
            cohensd_res[ct, "AEI_pvalue"] <- t.test(group_early$AEI_A2GEditing, group_late$AEI_A2GEditing, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR1exp[!is.na(group_early$ADAR1exp)]) > 1 & length(group_late$ADAR1exp[!is.na(group_late$ADAR1exp)]) > 1 ) {
            cohensd_res[ct, "ADAR1_pvalue"] <- t.test(group_early$ADAR1exp, group_late$ADAR1exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR2exp[!is.na(group_early$ADAR2exp)]) > 1 & length(group_late$ADAR2exp[!is.na(group_late$ADAR2exp)]) > 1 ) {
            cohensd_res[ct, "ADAR2_pvalue"] <- t.test(group_early$ADAR2exp, group_late$ADAR2exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
        if (length(group_early$ADAR3exp[!is.na(group_early$ADAR3exp)]) > 1 & length(group_late$ADAR3exp[!is.na(group_late$ADAR3exp)]) > 1 ) {
            cohensd_res[ct, "ADAR3_pvalue"] <- t.test(group_early$ADAR3exp, group_late$ADAR3exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)$p.value
        }
    }
}

write.table(cohensd_res, "Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_4wpc_11wpc_20231020.txt", sep = ",", row.names = T, quote = F)

csd_organ_order <- c("Forebrain/Cerebrum","Hindbrain/Cerebellum","Heart","Kidney","Liver","Testis")
organ_fill_col <- c("Forebrain/Cerebrum"="#990000","Hindbrain/Cerebellum"="#FF0000","Heart"="#FF9933",
"Kidney"="#66FF33","Liver"="#66FFFF","Testis"="#9933FF")
cohensd_res <- cohensd_res[row.names(cohensd_res) != "Ovary",]
bound_limit <- ceiling(max(abs(cohensd_res[!is.na(cohensd_res)])))

pdf("Figure_OrganDevelopment_R1/OD_cohensd_logTPM_effsizebtwTime_plot_4wpc_11wpc_20231012.pdf", height = 12, width = 25)
cowplot:::plot_grid(
ggplot(cohensd_res[,c(1:4)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = AEI_lower_ci, 
    upper = AEI_upper_ci, 
    ymin = AEI_lower_ci, 
    ymax = AEI_upper_ci, 
    middle = AEI_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit,bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('AEI'),
ggplot(cohensd_res[,c(5:8)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR1_lower_ci, 
    upper = ADAR1_upper_ci, 
    ymin = ADAR1_lower_ci, 
    ymax = ADAR1_upper_ci, 
    middle = ADAR1_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit,bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADAR1 Expression in TPM'),
ggplot(cohensd_res[,c(9:12)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR2_lower_ci, 
    upper = ADAR2_upper_ci, 
    ymin = ADAR2_lower_ci, 
    ymax = ADAR2_upper_ci, 
    middle = ADAR2_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit,bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADARB1 Expression in TPM'),
ggplot(cohensd_res[,c(13:16)], aes(x = factor(row.names(cohensd_res), levels = csd_organ_order), fill = organ_fill_col)) +
  geom_boxplot(aes(lower = ADAR3_lower_ci, 
    upper = ADAR3_upper_ci, 
    ymin = ADAR3_lower_ci, 
    ymax = ADAR3_upper_ci, 
    middle = ADAR3_d_estimate), stat = "identity", width = 0.5, show.legend = FALSE) + theme_classic() + coord_flip() + ylim(-bound_limit,bound_limit) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Organ') +
    ggtitle('ADARB2 Expression in TPM'),
ncol = 4
)
dev.off()

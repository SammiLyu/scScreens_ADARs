
load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_logNorm_int_annot_v1_col_v1_20230124.rda")
aei_table_df <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_allRNA_AEI_byCT_byExp_20230510.txt', header = T, sep =',')
adar_tpm <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_allRNA_adarTPM_bySample_byTer_20230128.txt', header = T, sep =',')

obj.int.allrna.sub
table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)
sam_tot <- paste0(rep(colnames(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)), 
                        each = nrow(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1))),
                    "_",
                  rep(row.names(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)), 
                        ncol(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1))))


# ------------------------------------------------------------------
# option 2 - logTPM

# formulate "meta data"
metadata_df <- data.frame(matrix(NA, nrow = length(sam_tot), ncol = 7))
row.names(metadata_df) <- sam_tot
colnames(metadata_df) <- c("ter","age","celltype","ADAR1exp","ADAR2exp","ADAR3exp","AEI_A2GEditing")
metadata_df$ter <- sapply(strsplit(row.names(metadata_df), split = "_"), "[[", 3)
metadata_df$age <- sapply(strsplit(row.names(metadata_df), split = "_"), "[[", 2)
metadata_df$celltype <- sapply(strsplit(row.names(metadata_df), split = "_"), "[[", 1)

for (i in 1:nrow(metadata_df)){
    if (length(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADAR", "TPM"]) != 0) {
        if (adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADAR", "TPM"] != 0){
            metadata_df$ADAR1exp[i] <- log2(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADAR", "TPM"])
        }
    }
    if (length(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB1", "TPM"]) != 0){
        if (adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB1", "TPM"] != 0){
            metadata_df$ADAR2exp[i] <- log2(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB1", "TPM"])
        }    
    }
    if (length(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB2", "TPM"]) != 0){
        if (adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB2", "TPM"] != 0){
            metadata_df$ADAR3exp[i] <- log2(adar_tpm[adar_tpm$Sample == row.names(metadata_df)[i] & adar_tpm$Gene == "ADARB2", "TPM"])
        }    
    }

    if (length(aei_table_df[row.names(metadata_df)[i], "A2GEditingIndex"]) != 0 ){
        metadata_df$AEI_A2GEditing[i] <- aei_table_df[row.names(metadata_df)[i], "A2GEditingIndex"]
    }
}
dim(metadata_df)

#metadata_df <- metadata_df[complete.cases(metadata_df),]
#dim(metadata_df)

# remove outlier M9
metadata_df <- metadata_df[metadata_df$ter != "M9",]

earlier_time <- c("wk4","wk6")
later_time <- c("wk8","wk10")

cohensd_res <- data.frame(matrix(NA, nrow = length(table(metadata_df$celltype)), ncol = 16))
row.names(cohensd_res) <- names(table(metadata_df$celltype))
colnames(cohensd_res) <- paste0(rep(c("AEI","ADAR1","ADAR2","ADAR3"), each = 4), "_", rep(c("d_estimate","lower_ci","upper_ci","pvalue")))

for (ct in names(table(metadata_df$celltype))){
    if (dim(metadata_df[metadata_df$celltype == ct,])[1] !=0){
        print(ct)

        group_early <- metadata_df[metadata_df$celltype == ct & metadata_df$age %in% earlier_time,  ]
        group_late <- metadata_df[metadata_df$celltype == ct & metadata_df$age %in% later_time,  ]
        
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
            cohensd_res[ct, "AEI_pvalue"] <- t.test(group_early$AEI_A2GEditing, group_late$AEI_A2GEditing, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95, na.rm = T)$p.value
        }
        if (length(group_early$ADAR1exp[!is.na(group_early$ADAR1exp)]) > 1 & length(group_late$ADAR1exp[!is.na(group_late$ADAR1exp)]) > 1 ) {
            cohensd_res[ct, "ADAR1_pvalue"] <- t.test(group_early$ADAR1exp, group_late$ADAR1exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95, na.rm = T)$p.value
        }
        if (length(group_early$ADAR2exp[!is.na(group_early$ADAR2exp)]) > 1 & length(group_late$ADAR2exp[!is.na(group_late$ADAR2exp)]) > 1 ) {
            cohensd_res[ct, "ADAR2_pvalue"] <- t.test(group_early$ADAR2exp, group_late$ADAR2exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95, na.rm = T)$p.value
        }
        if (length(group_early$ADAR3exp[!is.na(group_early$ADAR3exp)]) > 1 & length(group_late$ADAR3exp[!is.na(group_late$ADAR3exp)]) > 1 ) {
            cohensd_res[ct, "ADAR3_pvalue"] <- t.test(group_early$ADAR3exp, group_late$ADAR3exp, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95, na.rm = T)$p.value
        }
    }
}


bound <- ceiling(max(abs(min(cohensd_res, na.rm = T)), max(cohensd_res, na.rm = T)))
pdf("Figure2_TimeSeries_R1/TS_allRNA_cohensd_logTPM_plot_noM9_20230731.pdf", height = 8, width = 20)
plot_grid(
ggplot(cohensd_res[,c(1:4)], aes(x = row.names(cohensd_res))) +
  geom_boxplot(aes(lower = AEI_lower_ci, 
    upper = AEI_upper_ci, 
    ymin = AEI_lower_ci, 
    ymax = AEI_upper_ci, 
    middle = AEI_d_estimate), stat = "identity", width = 0.5) + theme_classic() + coord_flip() + ylim(-bound,bound) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Cell Type') +
    ggtitle('AEI'),
ggplot(cohensd_res[,c(5:8)], aes(x = row.names(cohensd_res))) +
  geom_boxplot(aes(lower = ADAR1_lower_ci, 
    upper = ADAR1_upper_ci, 
    ymin = ADAR1_lower_ci, 
    ymax = ADAR1_upper_ci, 
    middle = ADAR1_d_estimate), stat = "identity", width = 0.5) + theme_classic() + coord_flip() + ylim(-bound,bound) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Cell Type') +
    ggtitle('ADAR1 Expression in TPM'),
ggplot(cohensd_res[,c(9:12)], aes(x = row.names(cohensd_res))) +
  geom_boxplot(aes(lower = ADAR2_lower_ci, 
    upper = ADAR2_upper_ci, 
    ymin = ADAR2_lower_ci, 
    ymax = ADAR2_upper_ci, 
    middle = ADAR2_d_estimate), stat = "identity", width = 0.5) + theme_classic() + coord_flip() + ylim(-bound,bound) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Cell Type') +
    ggtitle('ADARB1 Expression in TPM'),
ggplot(cohensd_res[,c(13:16)], aes(x = row.names(cohensd_res))) +
  geom_boxplot(aes(lower = ADAR3_lower_ci, 
    upper = ADAR3_upper_ci, 
    ymin = ADAR3_lower_ci, 
    ymax = ADAR3_upper_ci, 
    middle = ADAR3_d_estimate), stat = "identity", width = 0.5) + theme_classic() + coord_flip() + ylim(-bound,bound) + 
    theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1), axis.title=element_text(size=18,face="bold")) + xlab('Cell Type') +
    ggtitle('ADARB2 Expression in TPM'), ncol = 4
)
dev.off()

write.table(cohensd_res, file = "Figure2_TimeSeries_R1/TS_allRNA_cohensd_AEI_vs_avgexp_noM9_20230731.txt", quote = F, sep = ",", col.names = T)







## adding in a normalization
metadata_df[,c(4:7)] <- apply(metadata_df[,c(4:7)],2,function(x){x/max(x)})


cohensd_res <- data.frame(matrix(NA, nrow = ncol(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)), ncol = 9))
row.names(cohensd_res) <- colnames(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1))
colnames(cohensd_res) <- paste0(rep(c("ADAR1","ADAR2","ADAR3"), each = 3), "_", rep(c("d_estimate","lower_ci","upper_ci")))

for (ct in names(table(metadata_df$celltype))){
    if (dim(metadata_df[metadata_df$celltype == ct,])[1] !=0){
        print(ct)

        cohensd_res[ct, "ADAR1_d_estimate"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR1exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$estimate
        cohensd_res[ct, "ADAR1_lower_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR1exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR1_upper_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR1exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR2_d_estimate"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR2exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$estimate
        cohensd_res[ct, "ADAR2_lower_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR2exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR2_upper_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR2exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR3_d_estimate"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR3exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$estimate
        cohensd_res[ct, "ADAR3_lower_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR3exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR3_upper_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR3exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["upper"]
    }
}

write.table(cohensd_res, file = "Figure2_TimeSeries_R1/TS_allRNA_cohensd_AEI_vs_logTPM_aftercleanup_20230623.txt", quote = F, sep = ",", col.names = T)

# ------------------------------------------------------------------
# option 3 - average expression from Seurat

avg_exp <- Seurat::AverageExpression(obj.int.allrna.sub, group.by = c("collapsed_annot_v1", "age_mouse_replicate"))$RNA[c("ADAR","ADARB1","ADARB2"),]
avg_exp <- t(avg_exp)

# formulate "meta data"
metadata_df <- data.frame(matrix(NA, nrow = length(sam_tot), ncol = 7))
row.names(metadata_df) <- sam_tot
colnames(metadata_df) <- c("ter","age","celltype","ADAR1exp","ADAR2exp","ADAR3exp","AEI_A2GEditing")
metadata_df$ter <- sapply(strsplit(row.names(metadata_df), split = "_"), "[[", 3)
metadata_df$age <- sapply(strsplit(row.names(metadata_df), split = "_"), "[[", 2)
metadata_df$celltype <- sapply(strsplit(row.names(metadata_df), split = "_"), "[[", 1)

for (i in 1:nrow(metadata_df)){

    if (row.names(metadata_df)[i] %in% row.names(avg_exp)){
    metadata_df[i,"ADAR1exp"] <- avg_exp[row.names(metadata_df)[i], "ADAR"]
    metadata_df[i,"ADAR2exp"] <- avg_exp[row.names(metadata_df)[i], "ADARB1"]
    metadata_df[i,"ADAR3exp"] <- avg_exp[row.names(metadata_df)[i], "ADARB2"]
    }

    if (length(aei_table_df[row.names(metadata_df)[i], "A2GEditingIndex"]) != 0 ){
        metadata_df$AEI_A2GEditing[i] <- aei_table_df[row.names(metadata_df)[i], "A2GEditingIndex"]
    }
}
dim(metadata_df)

metadata_df <- metadata_df[complete.cases(metadata_df),]
dim(metadata_df)

cohensd_res <- data.frame(matrix(NA, nrow = ncol(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)), ncol = 9))
row.names(cohensd_res) <- colnames(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1))
colnames(cohensd_res) <- paste0(rep(c("ADAR1","ADAR2","ADAR3"), each = 3), "_", rep(c("d_estimate","lower_ci","upper_ci")))

for (ct in names(table(metadata_df$celltype))){
    if (dim(metadata_df[metadata_df$celltype == ct,])[1] !=0){
        print(ct)

        cohensd_res[ct, "ADAR1_d_estimate"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR1exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$estimate
        cohensd_res[ct, "ADAR1_lower_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR1exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR1_upper_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR1exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR2_d_estimate"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR2exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$estimate
        cohensd_res[ct, "ADAR2_lower_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR2exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR2_upper_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR2exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["upper"]

        cohensd_res[ct, "ADAR3_d_estimate"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR3exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$estimate
        cohensd_res[ct, "ADAR3_lower_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR3exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["lower"]
        cohensd_res[ct, "ADAR3_upper_ci"] <- cohen.d(metadata_df[metadata_df$celltype == ct,"ADAR3exp"],
        metadata_df[metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T)$conf.int["upper"]
    }
}

write.table(cohensd_res, file = "Figure2_TimeSeries_R1/TS_allRNA_cohensd_AEI_vs_avgexp_20230623.txt", quote = F, sep = ",", col.names = T)


# ------------------------------------------------------------------
# currently I think the most sensible case would be using logTPM

cohensd_res <- read.table("Figure2_TimeSeries_R1/TS_allRNA_cohensd_AEI_vs_logTPM_aftercleanup_20230623.txt", quote = "", sep = ",", row.names = 1)


pdf("Figure2_TimeSeries_R1/RS_allRNA_cohensd_AEI_vs_logTPM_aftercleanup_afternorm_plot_20230623.pdf", height = 8, width = 20)
plot_grid(
ggplot(cohensd_res, aes(x=row.names(cohensd_res), y=ADAR1_d_estimate, fill = row.names(cohensd_res))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-5,5),
   ggplot(cohensd_res, aes(x=row.names(cohensd_res), y=ADAR2_d_estimate, fill = row.names(cohensd_res))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-5,5),
   ggplot(cohensd_res, aes(x=row.names(cohensd_res), y=ADAR3_d_estimate, fill = row.names(cohensd_res))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-5,5), ncol = 3
)
dev.off()


quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

ggplot(cohensd_res, aes(x=row.names(cohensd_res), y=ADAR1_d_estimate, fill = row.names(cohensd_res))) + 
   coord_flip() +
    stat_summary(fun.data = quantiles_95, geom="boxplot")

barplot(cohensd_res$ADAR1_d_estimate)
barplot(cohensd_res$ADAR2_d_estimate)
barplot(cohensd_res$ADAR3_d_estimate)



# ------------------------------------------------------------------
# test cases
group1 <- c(8, 9, 11, 11, 12, 14, 15, 16, 16, 18, 20, 21)
group2 <- c(7, 9, 10, 10, 11, 11, 12, 14, 14, 16, 20, 23)

library(effsize)
stat <- cohen.d(group1, group2, pooled = T)

library(lsr)
cohensD(group1,group2)

cohens_d <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1
    md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
    csd <- lx * var(x) + ly * var(y)
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation

    cd  <- md/csd                        ## cohen's d
    return(cd)
}
cohens_d(group1,group2)





for (ct in names(table(metadata_df$celltype))){
    for (t in names(table(metadata_df$age))){
        if (dim(metadata_df[metadata_df$age == t & metadata_df$celltype == ct,])[1] !=0){
            print(ct)
            print(t)

            print("adar1")
            print(cohen.d(metadata_df[metadata_df$age == t & metadata_df$celltype == ct,"ADAR1exp"],
            metadata_df[metadata_df$age == t & metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T))

            print("adar2")
            print(cohen.d(metadata_df[metadata_df$age == t & metadata_df$celltype == ct,"ADAR2exp"],
            metadata_df[metadata_df$age == t & metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T))

            print("adar3")
            print(cohen.d(metadata_df[metadata_df$age == t & metadata_df$celltype == ct,"ADAR3exp"],
            metadata_df[metadata_df$age == t & metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = T))
        }
    }
}






# -----------------------------------------------------------------
# another idea - time series cohen's d?

cohensd_res <- data.frame(matrix(NA, nrow = 4*3*ncol(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)), ncol = 4))
colnames(cohensd_res) <- c("cohensd","t","celltype","variable")
cohensd_res$t <- rep(rep(c("wk6_wk4","wk8_wk6","wk10_wk8"), each = 4),ncol(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)))
cohensd_res$celltype <- rep(colnames(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)), each = 4*3)
cohensd_res$variable <- rep(c("ADAR1","ADAR2","ADAR3","AEI"), 3*ncol(table(obj.int.allrna.sub$age_mouse_replicate, obj.int.allrna.sub$collapsed_annot_v1)))

for (ct in names(table(metadata_df$celltype))){
    for (t in c("wk6_wk4","wk8_wk6","wk10_wk8")){
        t1 <- strsplit(t, split = "_")[[1]][1]
        t2 <- strsplit(t, split = "_")[[1]][2]
        if (dim(metadata_df[metadata_df$age == t1 & metadata_df$celltype == ct,])[1] !=0 &
            dim(metadata_df[metadata_df$age == t2 & metadata_df$celltype == ct,])[1] !=0){
            print(ct)
            print(t)

            cohensd_res[cohensd_res$t == t & cohensd_res$celltype == ct & cohensd_res$variable == "ADAR1", "cohensd"] <- effsize:::cohen.d(metadata_df[metadata_df$age == t1 & metadata_df$celltype == ct,"ADAR1exp"],
            metadata_df[metadata_df$age == t2 & metadata_df$celltype == ct,"ADAR1exp"], pooled = T, paired = F, hedges.correction=T)$estimate
            cohensd_res[cohensd_res$t == t & cohensd_res$celltype == ct & cohensd_res$variable == "ADAR2", "cohensd"] <- effsize:::cohen.d(metadata_df[metadata_df$age == t1 & metadata_df$celltype == ct,"ADAR2exp"],
            metadata_df[metadata_df$age == t2 & metadata_df$celltype == ct,"ADAR2exp"], pooled = T, paired = F, hedges.correction=T)$estimate
            cohensd_res[cohensd_res$t == t & cohensd_res$celltype == ct & cohensd_res$variable == "ADAR3", "cohensd"] <- effsize:::cohen.d(metadata_df[metadata_df$age == t1 & metadata_df$celltype == ct,"ADAR3exp"],
            metadata_df[metadata_df$age == t2 & metadata_df$celltype == ct,"ADAR3exp"], pooled = T, paired = F, hedges.correction=T)$estimate
            cohensd_res[cohensd_res$t == t & cohensd_res$celltype == ct & cohensd_res$variable == "AEI", "cohensd"] <- effsize:::cohen.d(metadata_df[metadata_df$age == t1 & metadata_df$celltype == ct,"AEI_A2GEditing"],
            metadata_df[metadata_df$age == t2 & metadata_df$celltype == ct,"AEI_A2GEditing"], pooled = T, paired = F, hedges.correction=T)$estimate
        }
    }
}

cohensd_res$celltype <- factor(cohensd_res$celltype,levels(obj.int.allrna.sub$collapsed_annot_v1))

pdf("Figure2_TimeSeries_R1/RS_allRNA_cohensd_AEI_vs_logTPM_effsizebtwTime_plot_20230626.pdf", height = 8, width = 20)
plot_grid(
ggplot(cohensd_res[complete.cases(cohensd_res[,c(1:3)]) & cohensd_res$variable == "ADAR1",], aes(x=factor(celltype,levels(obj.int.allrna.sub$collapsed_annot_v1)), y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5),
ggplot(cohensd_res[complete.cases(cohensd_res[,c(1:3)]) & cohensd_res$variable == "ADAR2",], aes(x=celltype, y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5),
ggplot(cohensd_res[complete.cases(cohensd_res[,c(1:3)]) & cohensd_res$variable == "ADAR3",], aes(x=celltype, y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5), ncol = 3)
dev.off()


pdf("Figure2_TimeSeries_R1/RS_allRNA_cohensd_AEI_vs_logTPM_effsizebtwTime_plot_20230628.pdf", height = 12, width = 25)
plot_grid(
ggplot(cohensd_res[cohensd_res$variable == "ADAR1",], aes(x=factor(celltype,levels(obj.int.allrna.sub$collapsed_annot_v1)), y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5) + xlab('Cell Type') + ggtitle('ADAR1 Expression in TPM') + geom_hline(yintercept = c(-0.5, 0.5), linetype="dotdash") +
theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)),
ggplot(cohensd_res[cohensd_res$variable == "ADAR2",], aes(x=factor(celltype,levels(obj.int.allrna.sub$collapsed_annot_v1)), y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5) + xlab('Cell Type') + ggtitle('ADARB1 Expression in TPM') + geom_hline(yintercept = c(-0.5, 0.5), linetype="dotdash") +
theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)),
ggplot(cohensd_res[cohensd_res$variable == "ADAR3",], aes(x=factor(celltype,levels(obj.int.allrna.sub$collapsed_annot_v1)), y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5) + xlab('Cell Type') + ggtitle('ADARB2 Expression in TPM') + geom_hline(yintercept = c(-0.5, 0.5), linetype="dotdash") +
theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)), ncol = 3)
dev.off()

cohensd_res$celltype <- factor(cohensd_res$celltype, levels = sort(names(table(cohensd_res$celltype))))

pdf("Figure2_TimeSeries_R1/RS_allRNA_cohensd_AEI_vs_logTPM_effsizebtwTime_plot_alphabetical_20230628.pdf", height = 12, width = 25)
plot_grid(
ggplot(cohensd_res[cohensd_res$variable == "ADAR1",], aes(x=celltype, y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5) + xlab('Cell Type') + ggtitle('ADAR1 Expression in TPM') + geom_hline(yintercept = c(-0.5, 0.5), linetype="dotdash") +
theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)),
ggplot(cohensd_res[cohensd_res$variable == "ADAR2",], aes(x=celltype, y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5) + xlab('Cell Type') + ggtitle('ADARB1 Expression in TPM') + geom_hline(yintercept = c(-0.5, 0.5), linetype="dotdash") +
theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)),
ggplot(cohensd_res[cohensd_res$variable == "ADAR3",], aes(x=celltype, y=cohensd, fill = row.names(celltype))) + 
   coord_flip() + geom_boxplot() + theme_classic() + ylim(-1.5,1.5) + xlab('Cell Type') + ggtitle('ADARB2 Expression in TPM') + geom_hline(yintercept = c(-0.5, 0.5), linetype="dotdash") +
theme(axis.text=element_text(size=18,angle = 0, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)), ncol = 3)
dev.off()

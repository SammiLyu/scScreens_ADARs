library(cowplot)
library(ggplot2)
library(Seurat)
library(stringr)
library(dplyr)
library(cowplot)

library(limma)
library(edgeR)
library(stats)


### per cite RNA editing levels

### run detecting RNA editing pipeline with rediportal editing site table for reference GRCh38

load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")
KO_sub


### read in editing output and plot editing level per site with "significant" editing for each cell type
input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_KO_v7/"
edit_table <- list()
#table(KO$annot_v1, KO$guide)
#celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Hematopoietic",
#                    "Adipogenic-MSC-Fib","Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle",
#                    "Kidney-Prog")
exp_order <- names(table(KO_sub$exp))
guide_order <- names(table(KO_sub$guide))

  for (j in exp_order){
    for (k in guide_order){
      file_tmp <- paste0(input_dir, j,"_",k, "_dedup.txt")
      if (file.exists(file_tmp)) {
        print("file is there")
        edit_table_tmp <- read.table(file_tmp, header = 1)
        if (dim(edit_table_tmp)[1] == 0){
          print("empty file ignored")
        }
        else{
          edit_table_tmp$Exp <- j
          edit_table_tmp$Guide <- k
          edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
          edit_table[[paste0( j,"_",k)]] <- edit_table_tmp
        }
      }
    }
  }

length(edit_table)
QC_table <- data.frame(matrix(NA, nrow = length(edit_table), ncol = 3, dimnames = list(names(edit_table),c("prefilter","qcfiltered","snpfiltered"))))
QC_table$prefilter <- sapply(edit_table, dim)[1,]

# dbsnp filter
REDIref <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_GRCh38.txt", header = T, sep = "\t")
REDIref$Chrom_Pos <- paste0(REDIref$Region,"_",REDIref$Position)
head(REDIref)
dim(REDIref) # [1] 15681871       27
dim(REDIref[REDIref$dbsnp == "-",]) # [1] 13569494       27
SNP_sites <- REDIref[REDIref$dbsnp != "-","Chrom_Pos"]
head(SNP_sites)

#REDIref[REDIref$Chrom_Pos == "GRCh38_chr4_157336723", ] #GRIA2 editing site
#REDIref[REDIref$Chrom_Pos == "GRCh38_chrX_269035", ] #ADAR1 specific
#REDIref[REDIref$Chrom_Pos == "GRCh38_chrX_269037", ] #ADAR1 specific

#
#for (i in names(edit_table)){
#  print(i)
#  print(edit_table[[i]][edit_table[[i]]$Chrom_Pos == "GRCh38_chr4_157336723",])
#}

edit_table_cb <- do.call(rbind, edit_table)
edit_table_cb$Guide <- factor(edit_table_cb$Guide, levels = c("NTC","AAVS1","ADAR1","ADARB1","ADARB2"))

length(unique(edit_table_cb$Chrom_Pos))
edit_table_cb_avgter <- edit_table_cb %>% group_by(Chrom_Pos, Guide) %>%
    reframe(Chrom = Chrom, Position = Position, Gene = Gene, Conversion = Conversion, Alu = Alu, Region = Region, 
              TotalEditedReads = sum(TotalEditedReads), NonRefEdited = sum(NonRefEdited), 
              EditLevel = mean(EditLevel), ColEditLevel = sum(NonRefEdited)/sum(TotalEditedReads),Guide = Guide) %>%
              #CollapsedTotalEditedReads = toString(TotalEditedReads), CollapsedNonRefEdited = toString(NonRefEdited), 
              #MeanEditLevel = toString(EditLevel)) %>%
    #ungroup()
    distinct(Chrom_Pos, Guide, .keep_all = T) #%>% distinct(Guide, .keep_all = T)
dim(edit_table_cb_avgter)

edit_table_qc <- edit_table_cb_avgter
# remove editing sites with fewer than 3 edited reads
# remove all 0 editing level sites
edit_table_qc <- edit_table_qc[edit_table_qc$NonRefEdited >= 5,]
# filter for # of reads per site (at least 5 (Cellular))
edit_table_qc <- edit_table_qc[edit_table_qc$TotalEditedReads >= 10,]

dim(edit_table_qc)
table(edit_table_qc$Guide)

edit_table_qc <- edit_table_qc[!edit_table_qc$Chrom_Pos %in% SNP_sites,]
table(edit_table_qc$Guide)

write.table(edit_table_qc, "Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_goodqualitysites_ALL_20240715.txt",sep = ",", quote = F, col.names = T)

sites_to_lookat <- unique(edit_table_qc$Chrom_Pos)
length(unique(sites_to_lookat))

edit_table_sites_to_lookat <- as.data.frame(edit_table_cb_avgter[edit_table_cb_avgter$Chrom_Pos %in% sites_to_lookat,])
edit_table_sites_to_lookat <- edit_table_sites_to_lookat[edit_table_sites_to_lookat$TotalEditedReads >= 10,]
table(edit_table_sites_to_lookat$Guide)

    x <- list(NTC_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "NTC", "Chrom_Pos"],
    AAVS1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "AAVS1", "Chrom_Pos"],
    ADAR1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "ADAR1", "Chrom_Pos"],
    ADARB1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "ADARB1", "Chrom_Pos"],
    ADARB2_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "ADARB2", "Chrom_Pos"])
    p1 <- ggVennDiagram:::ggVennDiagram(x, label_alpha = 0) + ggtitle("Alu sites") +   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

    x <- list(NTC_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "NTC", "Chrom_Pos"],
    AAVS1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "AAVS1", "Chrom_Pos"],
    ADAR1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "ADAR1", "Chrom_Pos"],
    ADARB1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "ADARB1", "Chrom_Pos"],
    ADARB2_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "ADARB2", "Chrom_Pos"])
    p2 <- ggVennDiagram:::ggVennDiagram(x, label_alpha = 0) + ggtitle("Repetitive sites") +   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

    x <- list(NTC_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "NTC", "Chrom_Pos"],
    AAVS1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "AAVS1", "Chrom_Pos"],
    ADAR1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "ADAR1", "Chrom_Pos"],
    ADARB1_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "ADARB1", "Chrom_Pos"],
    ADARB2_sites = edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "ADARB2", "Chrom_Pos"])
    p3 <- ggVennDiagram:::ggVennDiagram(x, label_alpha = 0) + ggtitle("Non-Repetitive sites") +   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

pdf("Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_x10y5_venn_20230928.pdf", width = 18, height = 10)
    p1+p2+p3
dev.off()

sites_all1 <- Reduce(intersect, list(edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "NTC", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "AAVS1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "ADAR1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "ADARB1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "ALU" & edit_table_sites_to_lookat$Guide == "ADARB2", "Chrom_Pos"]))
    length(sites_all1)
    edit_table_sites_to_lookat_noNTCunique1 <- edit_table_sites_to_lookat[edit_table_sites_to_lookat$Chrom_Pos %in% sites_all1,]
  #write.table(sites_all1, "Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_goodqualitysites_ALU_20230928.txt", sep = ",", quote = F, col.names = F)
#sites_all1 <- read.table("Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_goodqualitysites_ALU_20230928.txt", sep = ",", quote = "", header = 1)

    sites_all2 <- Reduce(intersect, list(edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "NTC", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "AAVS1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "ADAR1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "ADARB1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "REP" & edit_table_sites_to_lookat$Guide == "ADARB2", "Chrom_Pos"]))
    length(sites_all2)
    edit_table_sites_to_lookat_noNTCunique2 <- edit_table_sites_to_lookat[edit_table_sites_to_lookat$Chrom_Pos %in% sites_all2,]
  #write.table(sites_all2, "Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_goodqualitysites_REP_20230928.txt", sep = ",", quote = F, col.names = F)
#sites_all2 <- read.table("Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_goodqualitysites_REP_20230928.txt", sep = ",", quote = "", header = 1)

    sites_all3 <- Reduce(intersect, list(edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "NTC", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "AAVS1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "ADAR1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "ADARB1", "Chrom_Pos"],
    edit_table_sites_to_lookat[edit_table_sites_to_lookat$Alu == "NONREP" & edit_table_sites_to_lookat$Guide == "ADARB2", "Chrom_Pos"]))
    length(sites_all3)
    edit_table_sites_to_lookat_noNTCunique3 <- edit_table_sites_to_lookat[edit_table_sites_to_lookat$Chrom_Pos %in% sites_all3,]
  #write.table(sites_all1, "Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_goodqualitysites_NONREP_20230928.txt", sep = ",", quote = F, col.names = F)
#sites_all3 <- read.table("Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_goodqualitysites_NONREP_20230928.txt", sep = ",", quote = "", header = 1)

    edit_table_sites_to_lookat_noNTCunique <- do.call(rbind, list(edit_table_sites_to_lookat_noNTCunique1,edit_table_sites_to_lookat_noNTCunique2,edit_table_sites_to_lookat_noNTCunique3))
write.table(edit_table_sites_to_lookat_noNTCunique, "Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_site_specific_editing_table_20231009.txt", sep = ",", quote = F, col.names = T)
edit_table_sites_to_lookat_noNTCunique <- read.table("Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_site_specific_editing_table_20231009.txt", quote = "", sep = ",", header = T)
## add in effect size

cohensd_sitespecific <- data.frame(matrix(NA, nrow = 3, ncol = 3, dimnames = list(c("ADAR1","ADARB1","ADARB2"), paste0(c("ALU","REP","NONREP"), rep("_effsize", 3)))))

for (i in row.names(cohensd_sitespecific)){
  for (j in colnames(cohensd_sitespecific)){
    region_tmp <- strsplit(j, split = "_")[[1]][1]
  cohensd_sitespecific[i,j] <- effsize:::cohen.d(edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == region_tmp & edit_table_sites_to_lookat_noNTCunique$Guide == "AAVS1", "EditLevel"],
edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == region_tmp & edit_table_sites_to_lookat_noNTCunique$Guide == i, "EditLevel"], pooled = T, na.rm = T)$estimate
  }
}

cohensd_sitespecific
write.table(cohensd_sitespecific, "Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_site_specific_editing_effsize_comptoAAVS1_20231009.txt", sep = ",", quote = F, col.names = F)


    means_alu <- aggregate(ColEditLevel ~  Guide, edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "ALU",], mean)
    means_alu[,2] <- round(means_alu[,2],3)

    means_rep <- aggregate(ColEditLevel ~  Guide, edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "REP",], mean)
    means_rep[,2] <- round(means_rep[,2],3)

    means_nr <- aggregate(ColEditLevel ~  Guide, edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "NONREP",], mean)
    means_nr[,2] <- round(means_nr[,2],3)

pdf("Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_x10y5_scatter_colEditLevel_20230928.pdf", width = 18, height = 10)
   plot_grid(
        ggplot(edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "ALU" ,], aes(x=Guide, y=ColEditLevel, fill=Guide)) + geom_jitter(colour = "grey", shape=16, size = 1,
                                                            position=position_jitter(0.2)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = c("#FF8000","#0000FF","#FF0000","#00C000","#AD07E3")) + coord_cartesian(ylim=c(0, 0.3)) + 
                                                            stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) + 
                                                            geom_text(data = means_alu, aes(label = ColEditLevel, y = ColEditLevel + 0.005)) + ggtitle("Alu Repeats"),
        ggplot(edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "REP" ,], aes(x=Guide, y=ColEditLevel, fill=Guide)) + geom_jitter(colour = "grey", shape=16, size = 1,
                                                            position=position_jitter(0.2)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = c("#FF8000","#0000FF","#FF0000","#00C000","#AD07E3")) + coord_cartesian(ylim=c(0, 0.3)) + 
                                                            stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) + 
                                                            geom_text(data = means_rep, aes(label = ColEditLevel, y = ColEditLevel + 0.005)) + ggtitle("Repetitive Regions"),
        ggplot(edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "NONREP",], aes(x=Guide, y=ColEditLevel, fill=Guide)) + geom_jitter(colour = "grey", shape=16, size = 1,
                                                            position=position_jitter(0.2)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = c("#FF8000","#0000FF","#FF0000","#00C000","#AD07E3")) + coord_cartesian(ylim=c(0, 0.3)) + 
                                                            stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) + 
                                                            geom_text(data = means_nr, aes(label = ColEditLevel, y = ColEditLevel + 0.005)) + ggtitle("Non-repetitive Regions"),
        ncol = 3)
dev.off()

# pdf("Figure3_ADARKO_R1/KO_4ter_supervised_v10_collapseTer_x10y5_scatter_colEditLevel_violin_20230928.pdf", width = 18, height = 10)
#    plot_grid(
#         ggplot(edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "ALU" ,], aes(x=Guide, y=ColEditLevel, fill=Guide)) + geom_jitter(colour = "grey", shape=16, size = 1,
#                                                             position=position_jitter(0.2)) + geom_violin() + theme_classic() + scale_fill_manual(values = c("#FF8000","#0000FF","#FF0000","#00C000","#AD07E3")) + coord_cartesian(ylim=c(0, 0.3)) + 
#                                                             stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) + 
#                                                             geom_text(data = means_alu, aes(label = ColEditLevel, y = ColEditLevel + 0.005)) + ggtitle("Alu Repeats"),
#         ggplot(edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "REP" ,], aes(x=Guide, y=ColEditLevel, fill=Guide)) + geom_jitter(colour = "grey", shape=16, size = 1,
#                                                             position=position_jitter(0.2)) + geom_violin() + theme_classic() + scale_fill_manual(values = c("#FF8000","#0000FF","#FF0000","#00C000","#AD07E3")) + coord_cartesian(ylim=c(0, 0.3)) + 
#                                                             stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) + 
#                                                             geom_text(data = means_rep, aes(label = ColEditLevel, y = ColEditLevel + 0.005)) + ggtitle("Repetitive Regions"),
#         ggplot(edit_table_sites_to_lookat_noNTCunique[edit_table_sites_to_lookat_noNTCunique$Alu == "NONREP",], aes(x=Guide, y=ColEditLevel, fill=Guide)) + geom_jitter(colour = "grey", shape=16, size = 1,
#                                                             position=position_jitter(0.2)) + geom_violin() + theme_classic() + scale_fill_manual(values = c("#FF8000","#0000FF","#FF0000","#00C000","#AD07E3")) + coord_cartesian(ylim=c(0, 0.3)) + 
#                                                             stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) + 
#                                                             geom_text(data = means_nr, aes(label = ColEditLevel, y = ColEditLevel + 0.005)) + ggtitle("Non-repetitive Regions"),
#         ncol = 3)
# dev.off()

    anova_data_editingrate <- edit_table_sites_to_lookat_noNTCunique[which(edit_table_sites_to_lookat_noNTCunique$Alu == "ALU") ,c("ColEditLevel","Guide")]
    # Compute the analysis of variance
    res.aov <- aov(ColEditLevel ~ Guide, data = anova_data_editingrate)
    # Summary of the analysis
    summary(res.aov)
    #library(knitr)
    #kable(summary(res.aov), digits = 3)
    #write.table(summary(res.aov), file = 'Figure3_ADARKO_R1/KO_4ter_supervised_editingrate_aov_summary.txt', quote = F, sep = ',')
    # to know which means are different
    TukeyHSD(res.aov)

    anova_data_editingrate <- edit_table_sites_to_lookat_noNTCunique[which(edit_table_sites_to_lookat_noNTCunique$Alu == "REP") ,c("ColEditLevel","Guide")]
    # Compute the analysis of variance
    res.aov <- aov(ColEditLevel ~ Guide, data = anova_data_editingrate)
    # Summary of the analysis
    summary(res.aov)
    #library(knitr)
    #kable(summary(res.aov), digits = 3)
    #write.table(summary(res.aov), file = 'Figure3_ADARKO_R1/KO_4ter_supervised_editingrate_aov_summary.txt', quote = F, sep = ',')
    # to know which means are different
    TukeyHSD(res.aov)

    anova_data_editingrate <- edit_table_sites_to_lookat_noNTCunique[which(edit_table_sites_to_lookat_noNTCunique$Alu == "NONREP") ,c("ColEditLevel","Guide")]
    # Compute the analysis of variance
    res.aov <- aov(ColEditLevel ~ Guide, data = anova_data_editingrate)
    # Summary of the analysis
    summary(res.aov)
    #library(knitr)
    #kable(summary(res.aov), digits = 3)
    #write.table(summary(res.aov), file = 'Figure3_ADARKO_R1/KO_4ter_supervised_editingrate_aov_summary.txt', quote = F, sep = ',')
    # to know which means are different
    TukeyHSD(res.aov)


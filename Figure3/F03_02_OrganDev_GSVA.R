### Newest GSVA pipeline requires R > 4.3.2 for gsvaParam to work
### installing R=4.3.2 in local home directory and will run everything in this script command line

# /media/Home_Raid1_Voyager/sammi/R/4.2.3/bin/R
# install.packages("remotes")
# library(remotes)
# install_github("rcastelo/GSVA")

# test code below works YAY!!
library(GSVA)
p <- 10000 ## number of genes
n <- 30    ## number of samples
## simulate expression values from a standard Gaussian distribution
X <- matrix(rnorm(p*n), nrow=p,
            dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
## sample gene set sizes
gs <- as.list(sample(10:100, size=100, replace=TRUE))
## sample gene sets
gs <- lapply(gs, function(n, p)
                   paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))
gsvaPar <- gsvaParam(X, gs)
gsva.es <- gsva(gsvaPar, verbose=FALSE)

# install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("GSEABase")
library(org.Hs.eg.db)
goannot <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")
head(goannot)
table(goannot$ONTOLOGY) # BP - Biological Process, Cellular Component, Molecular Function
genesbygo <- split(goannot$ENTREZID, goannot$GO)
length(genesbygo)
head(genesbygo)

library(GSEABase)
library(GSVAdata)

data(c2BroadSets)
class(c2BroadSets)

rpkm <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.txt")
colnames(rpkm)[which(colnames(rpkm) == "Testis.Senior.41")] <- "Testis.senior.41"
#write.table(rpkm, "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.Updated.txt")
#rpkm <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.Updated.txt")
rpkm_entrezids <- mapIds(org.Hs.eg.db, keys = row.names(rpkm), keytype="ENSEMBL", column = "ENTREZID")
names(rpkm_entrezids) <- NULL


log_rpkm <- as.matrix(log2(rpkm+1), nrow = dim(rpkm)[1], ncol = dim(rpkm)[2])
row.names(log_rpkm) <- rpkm_entrezids

Parlogrpkm <- gsvaParam(log_rpkm,
                         c2BroadSets, minSize=5, maxSize=500)
esrnaseq <- gsva(Parlogrpkm)

# colnames(esrnaseq)[which(colnames(esrnaseq) == "Testis.Senior.41")] <- "Testis.senior.41"
write.table(esrnaseq, "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_c2BS_20240227.txt", sep = ",", row.names = T, col.names = T, quote = F)
esrnaseq <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_c2BS_20240227.txt", sep = ",", row.names = 1, header = T, quote = "")


C2gs <- getGmt("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/GeneSets_msigdb/c2.all.v2023.2.Hs.entrez.gmt",
                 collectionType=BroadCollection(category="c2"),
                 geneIdType=EntrezIdentifier())

Parlogrpkm <- gsvaParam(log_rpkm,
                         C2gs, minSize=5, maxSize=Inf)
esrnaseq <- gsva(Parlogrpkm)

colnames(esrnaseq)[which(colnames(esrnaseq) == "Testis.Senior.41")] <- "Testis.senior.41"
write.table(esrnaseq, "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_C2gs_maxInf_20240229.txt", sep = ",", row.names = T, col.names = T, quote = F)
esrnaseq <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_C2gs_maxInf_20240229.txt", sep = ",", row.names = 1, header = T, quote = "")


OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1.csv', sep = ',', header = T)
OrganDev_meta <- OrganDev_meta[OrganDev_meta$Used.library. == "Yes",]

head(OrganDev_meta)
OrganDev_meta$Developmental_stage_sum <- NA
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc'), "Developmental_stage_sum"] <- "Early Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc'), "Developmental_stage_sum"] <- "Late Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('newborn','infant','toddler','school','youngTeenager','teenager'), "Developmental_stage_sum"] <- "Newborn - Teenager"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'), "Developmental_stage_sum"] <- "Adult - Senior"

OrganDev_meta$Stage <- NA
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Early Gestation","Late Gestation"), "Stage"] <- "Prenatal"
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Newborn - Teenager","Adult - Senior"), "Stage"] <- "Postnatal"

OrganDev_meta$Developmental.stage <- factor(OrganDev_meta$Developmental.stage, levels = c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc',
'11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc',
'newborn','infant','toddler','school','youngTeenager','teenager',
'oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'
))

OrganDev_meta <- OrganDev_meta[order(OrganDev_meta$Organ, OrganDev_meta$DaysPostConception),]
colnames(rpkm)
dim(rpkm)

dim(OrganDev_meta)

OrganDev_meta$sample_rpkm <- NA
 

colnames(rpkm)[grep("Brain",colnames(rpkm))]
OrganDev_meta[OrganDev_meta$Organ == "Forebrain/Cerebrum","Developmental.stage"]
OrganDev_meta[OrganDev_meta$Organ == "Forebrain/Cerebrum","sample_rpkm"] <- paste0("Brain.",OrganDev_meta[OrganDev_meta$Organ == "Forebrain/Cerebrum","Developmental.stage"],".",sapply(strsplit(colnames(rpkm)[grep("Brain",colnames(rpkm))], split = "[.]"), "[[",3))
table(OrganDev_meta[OrganDev_meta$Organ == "Forebrain/Cerebrum","Used.library."])

table(sapply(strsplit(colnames(rpkm)[grep("Cerebellum",colnames(rpkm))], split = "[.]"),"[[",2))
table(OrganDev_meta[OrganDev_meta$Organ == "Hindbrain/Cerebellum","Developmental.stage"])
# 1 newborn sample missing
OrganDev_meta[OrganDev_meta$Organ == "Hindbrain/Cerebellum" & OrganDev_meta$Developmental.stage == "newborn", ]
OrganDev_meta[OrganDev_meta$Organ == "Hindbrain/Cerebellum","sample_rpkm"] <- paste0("Cerebellum.",OrganDev_meta[OrganDev_meta$Organ == "Hindbrain/Cerebellum","Developmental.stage"],".",sapply(strsplit(colnames(rpkm)[grep("Cerebellum",colnames(rpkm))], split = "[.]"), "[[",3))

OrganDev_meta[OrganDev_meta$Organ == "Heart","sample_rpkm"] <- paste0("Heart.",OrganDev_meta[OrganDev_meta$Organ == "Heart","Developmental.stage"],".",sapply(strsplit(colnames(rpkm)[grep("Heart",colnames(rpkm))], split = "[.]"), "[[",3))
OrganDev_meta[OrganDev_meta$Organ == "Liver","sample_rpkm"] <- paste0("Liver.",OrganDev_meta[OrganDev_meta$Organ == "Liver","Developmental.stage"],".",sapply(strsplit(colnames(rpkm)[grep("Liver",colnames(rpkm))], split = "[.]"), "[[",3))
OrganDev_meta[OrganDev_meta$Organ == "Kidney","sample_rpkm"] <- paste0("Kidney.",OrganDev_meta[OrganDev_meta$Organ == "Kidney","Developmental.stage"],".",sapply(strsplit(colnames(rpkm)[grep("Kidney",colnames(rpkm))], split = "[.]"), "[[",3))
OrganDev_meta[OrganDev_meta$Organ == "Testis","sample_rpkm"] <- paste0("Testis.",OrganDev_meta[OrganDev_meta$Organ == "Testis","Developmental.stage"],".",sapply(strsplit(colnames(rpkm)[grep("Testis",colnames(rpkm))], split = "[.]"), "[[",3))

OrganDev_meta$sample_rpkm

setdiff(OrganDev_meta$sample_rpkm, colnames(rpkm))
setdiff(colnames(rpkm), OrganDev_meta$sample_rpkm)

# drop ovary from all tables
OrganDev_meta <- OrganDev_meta[OrganDev_meta$Organ != "Ovary",]
rpkm <- rpkm[, grep("Ovary", colnames(rpkm), invert = T)]


aei_table_df <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_OrganDev_v2/output_index_20230607/OrganDev_v1_outputindex/EditingIndex.csv', header = T, sep =',')
row.names(aei_table_df) <- aei_table_df$Sample
aei_table_df$Library <- sapply(strsplit(row.names(aei_table_df), split = "sTS"), "[[",1)
aei_table_df <- aei_table_df[aei_table_df$Library %in% OrganDev_meta$library.ID, ]
dim(aei_table_df)

OrganDev_meta <- OrganDev_meta[OrganDev_meta$library.ID %in% aei_table_df$Library,]
aei_table_df <- aei_table_df[match(OrganDev_meta$library.ID, aei_table_df$Library),]

OrganDev_meta$AEI <- aei_table_df$A2GEditingIndex
write.table(OrganDev_meta, 'Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1_withAEI.csv', sep = ',', row.names = T, col.names = T, quote = F)
OrganDev_meta <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1_withAEI.csv', sep = ',', row.names = 1, header = T, quote = "")

# esrnaseq <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_c2BS_20240227.txt", sep = ",", row.names = 1, header = T, quote = "")
# C2
esrnaseq <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_C2gs_maxInf_20240229.txt", sep = ",", row.names = 1, header = T, quote = "")

esrnaseq <- esrnaseq[, grep("Ovary", colnames(esrnaseq), invert = T)]

esrnaseq <- esrnaseq[, colnames(esrnaseq) %in% OrganDev_meta$sample_rpkm]
esrnaseq <- esrnaseq[, match(OrganDev_meta$sample_rpkm, colnames(esrnaseq))]

# generate a plot of esrnaseq score vs. AEI

pathway <- c("REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE")
pathway <- c("REACTOME_INNATE_IMMUNE_SYSTEM")

pathway_all <- row.names(esrnaseq)[grep("IMMUNE",row.names(esrnaseq))]
pathway_all <- row.names(esrnaseq)[grep("INTERFERON",row.names(esrnaseq))]

#pathway_all <- c("REACTOME_VIRAL_INFECTION_PATHWAYS","REACTOME_INFECTIOUS_DISEASE","REACTOME_DNA_REPLICATION","REACTOME_SARS_COV_2_INFECTION","REACTOME_SYNTHESIS_OF_DNA")
pathway_all <- c("REACTOME_VIRAL_INFECTION_PATHWAYS","REACTOME_DNA_REPLICATION")
esrnaseq[pathway_all,]

write.table(esrnaseq[pathway_all,], file = "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_C2gs_maxInf_20240229_viralinfecanddnareplic_gs_20240326.txt", quote = F, col.names = T, row.names = T, sep = ",")

lat.order.levels <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc','11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc','newborn','infant','toddler','school','youngTeenager','teenager','oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')
DevelopmentStage.data <- data.frame(
  DevelopmentStage = lat.order.levels,
  DevelopmentStage_num = 1:length(lat.order.levels)
) # metadata obtained

library(ggplot2)
for (pathway in pathway_all){

p <- list()
for (i in names(table(sapply(strsplit(colnames(esrnaseq), split = "[.]"), "[[", 1)))){
    esrnaseq_sub <- esrnaseq[,sapply(strsplit(colnames(esrnaseq), split = "[.]"), "[[", 1) == i]
    organmeta_sub <- OrganDev_meta[OrganDev_meta$sample_rpkm %in% colnames(esrnaseq_sub),]

data_to_plot <- data.frame(gs = t(esrnaseq_sub[pathway,]),
aei = organmeta_sub$AEI,
DevelopmentStage = organmeta_sub$Developmental.stage)
colnames(data_to_plot)[1] <- "gs"

cor_results_tmp <- cor.test(data_to_plot$gs, data_to_plot$aei, method=("pearson"))
annotations <- data.frame(
        xpos = c(1),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_results_tmp$estimate, digits = 3), "\n ", "p = ", format(cor_results_tmp$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust 

data_to_plot <- dplyr::left_join(data_to_plot, DevelopmentStage.data, by= "DevelopmentStage")
data_to_plot <- data_to_plot[,which(colnames(data_to_plot) != "DevelopmentStage")]

min_gs <- floor(min(data_to_plot$gs))
max_gs <- ceiling(max(data_to_plot$gs))

min_aei <- floor(min(data_to_plot$aei))
max_aei <- ceiling(max(data_to_plot$aei))

coeff <- (max_aei - min_aei) / ( max_gs - min_gs) ; coeff


data_to_plot <- reshape:::melt(data_to_plot, id=c("DevelopmentStage_num"))
data_to_plot <- as.data.frame(data_to_plot)
print(data_to_plot)
tail(data_to_plot)

data_to_plot$variable <- factor(data_to_plot$variable, levels = c("gs","aei"))

data_to_plot[data_to_plot$variable == "aei","value"] <- data_to_plot[data_to_plot$variable == "aei","value"] / coeff + min_gs
tail(data_to_plot)


p[[i]] <- ggplot(data_to_plot, aes(x = DevelopmentStage_num, y = value, color = variable)) + geom_point(size=5) + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + scale_y_continuous(
    # Features of the first axis
    name = "Gene Score",
    limits = c(min_gs, max_gs),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~((.-min_gs)*coeff), name="AEI")
  ) + coord_cartesian(xlim = c(1,25)) + 
  geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") + #+ ylim(min_exp, max_exp)
  theme_classic() + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), inherit.aes = FALSE) + ggtitle(paste0("Gene Score and AEI for ", i)) +
  # original font size: 18, larger font size: 36 
theme(axis.text=element_text(size=36,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=36,face="bold"), axis.title.y.right = element_text(size=36,face="bold"),
legend.title = element_text(size=36), #change legend title font size
        legend.text = element_text(size=16) , legend.position = "none") 
}

pdf(paste0("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/plots/C2gs_20240229/",pathway,"_byorgan_largerfonts.pdf"), width = 17, height = 12)
print(p)
dev.off()
}


data_to_plot <- data.frame(gs = t(esrnaseq[pathway,]),
aei = OrganDev_meta$AEI,
DevelopmentStage = OrganDev_meta$Developmental.stage)
colnames(data_to_plot)[1] <- "gs"
data_to_plot <- dplyr::left_join(data_to_plot, DevelopmentStage.data, by= "DevelopmentStage")
data_to_plot <- data_to_plot[,which(colnames(data_to_plot) != "DevelopmentStage")]


data_to_plot <- reshape:::melt(data_to_plot, id=c("DevelopmentStage_num"))
data_to_plot <- as.data.frame(data_to_plot)

library(ggplot2)
pdf(paste0("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/plots/",pathway,".pdf"), width = 10, height = 10)
ggplot(data_to_plot, aes(x = DevelopmentStage_num, y = value, color = variable)) + geom_point(size=5) + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + 
  theme_classic() + ggtitle(paste0("Gene Score and AEI")) + coord_cartesian(xlim = c(1,25)) + 
theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)
dev.off()






AEI <- as.matrix(OrganDev_meta$AEI, dimnames = list(colnames(esrnaseq), "AEI"))
row.names(AEI) <- colnames(esrnaseq)
colnames(AEI) <- "AEI"
cor_results <- cor(t(esrnaseq), AEI)
cor_results <- as.matrix(cor_results)
cor_results <- cor_results[order(cor_results, decreasing = T), ]
write.table(cor_results, "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_C2gs_maxInf_corresults_20240304.txt", sep = ",", row.names = T, col.names = F, quote = F)
cor_results <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/logrpkm_C2gs_maxInf_corresults_20240304.txt", sep = ",", header = F, quote = "")

pathway_all <- c("REACTOME_VIRAL_INFECTION_PATHWAYS","REACTOME_DNA_REPLICATION")

cor_results[cor_results$V1 %in% pathway_all, ]

#### a version of the plot to fit on the same figure


pathway_all <- c("REACTOME_VIRAL_INFECTION_PATHWAYS","REACTOME_DNA_REPLICATION")
esrnaseq[pathway_all,]

library(ggplot2)
for (pathway in pathway_all){

p <- list()
for (i in names(table(sapply(strsplit(colnames(esrnaseq), split = "[.]"), "[[", 1)))){
    esrnaseq_sub <- esrnaseq[,sapply(strsplit(colnames(esrnaseq), split = "[.]"), "[[", 1) == i]
    organmeta_sub <- OrganDev_meta[OrganDev_meta$sample_rpkm %in% colnames(esrnaseq_sub),]

data_to_plot <- data.frame(gs = t(esrnaseq_sub[pathway,]),
aei = organmeta_sub$AEI,
DevelopmentStage = organmeta_sub$Developmental.stage)
colnames(data_to_plot)[1] <- "gs"
data_to_plot <- dplyr::left_join(data_to_plot, DevelopmentStage.data, by= "DevelopmentStage")
data_to_plot <- data_to_plot[,which(colnames(data_to_plot) != "DevelopmentStage")]

min_gs <- floor(min(data_to_plot$gs))
max_gs <- ceiling(max(data_to_plot$gs))

min_aei <- floor(min(data_to_plot$aei))
max_aei <- ceiling(max(data_to_plot$aei))

coeff <- (max_aei - min_aei) / ( max_gs - min_gs) ; coeff


data_to_plot <- reshape:::melt(data_to_plot, id=c("DevelopmentStage_num"))
data_to_plot <- as.data.frame(data_to_plot)
tail(data_to_plot)

data_to_plot$variable <- factor(data_to_plot$variable, levels = c("gs","aei"))

data_to_plot[data_to_plot$variable == "aei","value"] <- data_to_plot[data_to_plot$variable == "aei","value"] / coeff + min_gs
tail(data_to_plot)


p[[i]] <- ggplot(data_to_plot, aes(x = DevelopmentStage_num, y = value, color = variable)) + geom_point(size=5) + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + scale_y_continuous(
    # Features of the first axis
    name = "Gene Score",
    limits = c(min_gs, max_gs),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~((.-min_gs)*coeff), name="AEI")
  ) + 
  theme_classic() + #ggtitle(paste0("Gene Score and AEI for ", i)) + 
  coord_cartesian(xlim = c(1,25)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=36), legend.position = "none") +
  # original font size: 18, larger font size: 36
#theme(axis.text=element_text(size=36,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=36,face="bold"), axis.title.y.right = element_text(size=36,face="bold"),
#legend.title = element_text(size=36), #change legend title font size
#        legend.text = element_text(size=16) , legend.position = "none") + 
        geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)
}

pdf(paste0("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/plots/C2gs_20240229/",pathway,"_byorgan_largerfonts_fixaxis.pdf"), width = 17, height = 7.5)
print(p)
dev.off()
}




#### a version of the plot to have both pathways


pathway_all <- c("REACTOME_VIRAL_INFECTION_PATHWAYS","REACTOME_DNA_REPLICATION")
esrnaseq[pathway_all,]

lat.order.levels <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc','11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc','newborn','infant','toddler','school','youngTeenager','teenager','oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')
DevelopmentStage.data <- data.frame(
  DevelopmentStage = lat.order.levels,
  DevelopmentStage_num = 1:length(lat.order.levels)
) # metadata obtained

library(ggplot2)

p <- list()
for (i in sapply(strsplit(colnames(esrnaseq), split = "[.]"), "[[", 1)){
    esrnaseq_sub <- esrnaseq[,sapply(strsplit(colnames(esrnaseq), split = "[.]"), "[[", 1) == i]
    organmeta_sub <- OrganDev_meta[OrganDev_meta$sample_rpkm %in% colnames(esrnaseq_sub),]

data_to_plot <- data.frame(gs = t(esrnaseq_sub[pathway_all,]),
aei = organmeta_sub$AEI,
DevelopmentStage = organmeta_sub$Developmental.stage)
#colnames(data_to_plot)[1] <- "gs"
data_to_plot <- dplyr::left_join(data_to_plot, DevelopmentStage.data, by= "DevelopmentStage")
data_to_plot <- data_to_plot[,which(colnames(data_to_plot) != "DevelopmentStage")]

min_gs <- floor(min(data_to_plot[,paste0("gs.",pathway_all[1])]))
max_gs <- ceiling(max(data_to_plot[,paste0("gs.",pathway_all[1])]))

min_aei <- floor(min(data_to_plot$aei))
max_aei <- ceiling(max(data_to_plot$aei))

coeff <- (max_aei - min_aei) / ( max_gs - min_gs) ; coeff


data_to_plot <- reshape:::melt(data_to_plot, id=c("DevelopmentStage_num"))
data_to_plot <- as.data.frame(data_to_plot)
tail(data_to_plot)

data_to_plot$variable <- factor(data_to_plot$variable, levels = c(paste0("gs.",pathway_all),"aei"))

write.table(data_to_plot, file = paste0("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/plots/C2gs_20240229/pathway_viralanddna_",i,".txt"), quote = F)

data_to_plot[data_to_plot$variable == "aei","value"] <- data_to_plot[data_to_plot$variable == "aei","value"] / coeff + min_gs
tail(data_to_plot)

print (i)
p[[i]] <- ggplot(data_to_plot, aes(x = DevelopmentStage_num, y = value, color = variable)) + geom_point(size=5) + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + scale_y_continuous(
    # Features of the first axis
    name = "Gene Score",
    limits = c(min_gs, max_gs),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~((.-min_gs)*coeff), name="AEI")
  ) + 
  theme_classic() + #ggtitle(paste0("Gene Score and AEI for ", i)) + 
  coord_cartesian(xlim = c(1,25)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=36), legend.position = "none") +
  # original font size: 18, larger font size: 36
#theme(axis.text=element_text(size=36,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=36,face="bold"), axis.title.y.right = element_text(size=36,face="bold"),
#legend.title = element_text(size=36), #change legend title font size
#        legend.text = element_text(size=16) , legend.position = "none") + 
        geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)
}

pdf(paste0("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure_OrganDevelopment_R1/GSVA/plots/C2gs_20240229/pathway_viralanddna_byorgan_largerfonts.pdf"), width = 17, height = 12)
print(p)
dev.off()


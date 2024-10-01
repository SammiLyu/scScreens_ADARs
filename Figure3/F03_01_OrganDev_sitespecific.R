library(ggplot2)
library(limma)
library(edgeR)

input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_OD_v1/"
edit_table <- list()

all_edit_table <- list.files(input_dir)

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    edit_table_tmp$Subject <- sapply(strsplit(i, split = "sTS"), "[[", 1)
    #edit_table_tmp$Organ <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = ".txt"), "[[", 1)]] <- edit_table_tmp
}

length(edit_table)

OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1.csv', sep = ',', header = T)
head(OrganDev_meta)
table(OrganDev_meta$Developmental.stage,OrganDev_meta$Organ)

for (i in names(edit_table)){
    edit_table[[i]][,"Organ"] <- OrganDev_meta[OrganDev_meta$library.ID == sapply(strsplit(i, split = "sTS"), "[[", 1), "Organ"]
    edit_table[[i]][,"Developmental_stage"] <- OrganDev_meta[OrganDev_meta$library.ID == sapply(strsplit(i, split = "sTS"), "[[", 1), "Developmental.stage"]
    edit_table[[i]][,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sapply(strsplit(i, split = "sTS"), "[[", 1), "Sex"]
}

head(edit_table[[1]])

QC_table <- data.frame(matrix(NA, nrow = length(edit_table), ncol = 3, dimnames = list(names(edit_table),c("prefilter","qcfiltered","snpfiltered"))))
QC_table$prefilter <- sapply(edit_table, dim)[1,]

# dbsnp filter
REDIref <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_hg19_mod_nochr.txt", header = T, sep = "\t")
REDIref$Chrom_Pos <- paste0(REDIref$Region,"_",REDIref$Position)
head(REDIref)
dim(REDIref) # [1] 15681871       27
dim(REDIref[REDIref$dbsnp == "-",]) # [1] 13569494       27
SNP_sites <- REDIref[REDIref$dbsnp != "-","Chrom_Pos"]; head(SNP_sites)
ACALAP_sites <- REDIref[grepl("^AC|^AL|^AP", REDIref$Gene.wgEncodeGencodeBasicV34lift37),"Chrom_Pos"]; head(ACALAP_sites)
MITO_sites <- REDIref[REDIref$Region == "M","Chrom_Pos"]; head(MITO_sites)

REDIref_sub <- REDIref[!REDIref$Chrom_Pos %in% c(SNP_sites,ACALAP_sites,MITO_sites),]
site_count_table <- data.frame(table(REDIref_sub$Func.wgEncodeGencodeBasicV34lift37))

site_count_table$Var1 <- factor(site_count_table$Var1, levels = c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic"))

pdf("Figure_OrganDevelopment_R1/REDIref_hg19_mod_nochr_geniclocations.pdf")
ggplot(data=site_count_table, aes(x=Var1, y=Freq, fill=Var1)) + scale_fill_brewer(palette="Set3") + 
   geom_bar(stat="identity") + theme_minimal() + ggtitle("Genic Locations of Sites in Reference") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) + 
     xlab("Region") + ylab("Number of Sites in log scale")
dev.off()

edit_table_cb <- do.call(rbind, edit_table)
edit_table_cb$Developmental_stage <- factor(edit_table_cb$Developmental_stage, levels = c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc','11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc','newborn','infant','toddler','school','youngTeenager','teenager','oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'))
edit_table_cb$Developmental_stage_sum <- NA
edit_table_cb[edit_table_cb$Developmental_stage %in% c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc'), "Developmental_stage_sum"] <- "Early Gestation"
edit_table_cb[edit_table_cb$Developmental_stage %in% c('11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc'), "Developmental_stage_sum"] <- "Late Gestation"
edit_table_cb[edit_table_cb$Developmental_stage %in% c('newborn','infant','toddler','school','youngTeenager','teenager'), "Developmental_stage_sum"] <- "Newborn - Teenager"
edit_table_cb[edit_table_cb$Developmental_stage %in% c('oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'), "Developmental_stage_sum"] <- "Adult - Senior"

write.table(edit_table_cb,"Figure_OrganDevelopment_R1/edit_table_cb_20231130.txt", sep = ",", col.names = T, row.names = T, quote = F)


hist(edit_table_cb$EditLevel)

### keeping the most recent version, same as v5 but run for all organs prenatal vs. postnatal

head(edit_table_cb)

length(unique(edit_table_EE_cb$Chrom_Pos)) # 1317
length(unique(intersect(edit_table_EE_cb$Chrom_Pos, edit_table_cb$Chrom_Pos))) # 451
length(unique(edit_table_MB_cb$Chrom_Pos)) # 36613
length(unique(intersect(edit_table_MB_cb$Chrom_Pos, edit_table_cb$Chrom_Pos))) # 34513

REDIref <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_hg19_mod_nochr.txt", header = T, sep = "\t")
REDIref$Chrom_Pos <- paste0(REDIref$Region,"_",REDIref$Position)
head(REDIref)
dim(REDIref) # [1] 15681871       27
dim(REDIref[REDIref$dbsnp == "-",]) # [1] 13569494       27
SNP_sites <- REDIref[REDIref$dbsnp != "-","Chrom_Pos"]
head(SNP_sites)

# check if some GRIK2 sites are in SNP_sites or not
GRIK2_sites <- c("6_102337689","6_102337702","6_102372572","6_102372589")
intersect(GRIK2_sites, SNP_sites)

####### Forebrain

# try checking site sparsity and then maybe imputation
edit_table_cb_forebrain <- edit_table_cb[edit_table_cb$Organ == "Forebrain/Cerebrum",]
dim(edit_table_cb_forebrain) # 29738662       16
length(unique(edit_table_cb_forebrain$Chrom_Pos)) # 2542727

#edit_table_cb_forebrain <- edit_table_cb_forebrain[edit_table_cb_forebrain$NonRefEdited >= 3,]
dim(edit_table_cb_forebrain) # 570500       16
edit_table_cb_forebrain <- edit_table_cb_forebrain[edit_table_cb_forebrain$TotalEditedReads >= 5,]
dim(edit_table_cb_forebrain) # 29738662       16

edit_table_cb_forebrain <- edit_table_cb_forebrain[edit_table_cb_forebrain$Chrom != "M",]
edit_table_cb_forebrain <- edit_table_cb_forebrain[!edit_table_cb_forebrain$Chrom_Pos %in% SNP_sites,]
edit_table_cb_forebrain <- edit_table_cb_forebrain[!grepl("^AC|^AL|^AP", edit_table_cb_forebrain$Gene),]

dim(edit_table_cb_forebrain) # 21598648       16
length(unique(edit_table_cb_forebrain$Chrom_Pos)) # 1830517

sample_cutoff = 0.6
editlevel_cutoff = 0.05

sample_count <- length(table(edit_table_cb_forebrain$Subject))
sample_count
highquality_sites_coverage <- names(table(edit_table_cb_forebrain$Chrom_Pos)[which(table(edit_table_cb_forebrain$Chrom_Pos) >= sample_cutoff*sample_count)])
length(highquality_sites_coverage) # 179832

highquality_sites_minediting_table <- aggregate(edit_table_cb_forebrain$EditLevel, list(edit_table_cb_forebrain$Chrom_Pos), FUN=mean)
head(highquality_sites_minediting_table)
highquality_sites_minediting <- highquality_sites_minediting_table[highquality_sites_minediting_table$x >= editlevel_cutoff, "Group.1"]
length(highquality_sites_minediting) # 185396

length(intersect(highquality_sites_coverage, highquality_sites_minediting)) # 5433
highquality_sites_forebrain <- intersect(highquality_sites_coverage, highquality_sites_minediting)

edit_table_cb_forebrain_qc <- edit_table_cb_forebrain[edit_table_cb_forebrain$Chrom_Pos %in% highquality_sites_forebrain, ]
dim(edit_table_cb_forebrain_qc)
summary(edit_table_cb_forebrain_qc)
write.table(edit_table_cb_forebrain_qc,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_forebrain_qc_imputedresults_forebrain_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)


OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1.csv', sep = ',', header = T)

head(OrganDev_meta)
OrganDev_meta$Developmental_stage_sum <- NA
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc'), "Developmental_stage_sum"] <- "Early Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc'), "Developmental_stage_sum"] <- "Late Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('newborn','infant','toddler','school','youngTeenager','teenager'), "Developmental_stage_sum"] <- "Newborn - Teenager"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'), "Developmental_stage_sum"] <- "Adult - Senior"

OrganDev_meta$Stage <- NA
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Early Gestation","Late Gestation"), "Stage"] <- "Prenatal"
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Newborn - Teenager","Adult - Senior"), "Stage"] <- "Postnatal"


sample_groups <- OrganDev_meta$library.ID
sites_of_interest <- highquality_sites_forebrain
editExprs <- data.frame(matrix(NA, nrow = length(sites_of_interest), ncol = length(sample_groups), dimnames = list(sites_of_interest, sample_groups)))
    colnames(editExprs) <- sample_groups
    for (k in 1:nrow(edit_table_cb_forebrain_qc)){
        editExprs_sites <- edit_table_cb_forebrain_qc[k, "Chrom_Pos"]
        editExprs_sample <- edit_table_cb_forebrain_qc[k, "Subject"]
        editExprs[editExprs_sites, as.character(editExprs_sample)] <- edit_table_cb_forebrain_qc[k, "EditLevel"]
    }
dim(editExprs)
highquality_samples <- colnames(is.na(editExprs))[colSums(is.na(editExprs)) < 0.2*nrow(editExprs)]


table(OrganDev_meta[OrganDev_meta$library.ID %in% highquality_samples, "Developmental_stage_sum"])

editExprs_hq <- editExprs[colnames(editExprs) %in% highquality_samples]
dim(editExprs_hq)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(editExprs_hq,2,pMiss)
max(apply(editExprs_hq,2,pMiss))
mean(apply(editExprs_hq,2,pMiss))
apply(editExprs_hq,1,pMiss)
max(apply(editExprs_hq,1,pMiss))
mean(apply(editExprs_hq,1,pMiss))

aggr_plot <- aggr(editExprs_hq, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))


# after dropping the low quality sample, clean up the sites one more time
highquality_sites_new <- names(apply(editExprs_hq,1,pMiss))[apply(editExprs_hq,1,pMiss) < (sample_cutoff*100)]
length(highquality_sites_new)
highquality_sites_new <- setdiff(highquality_sites_new, names(which(rowSums(editExprs_hq, na.rm = T) <= editlevel_cutoff)))
length(highquality_sites_new)

editExprs_hq_new <- editExprs_hq[row.names(editExprs_hq) %in% highquality_sites_new,]

apply(editExprs_hq_new,2,pMiss)
max(apply(editExprs_hq_new,2,pMiss))
mean(apply(editExprs_hq_new,2,pMiss))
apply(editExprs_hq_new,1,pMiss)
max(apply(editExprs_hq_new,1,pMiss))
mean(apply(editExprs_hq_new,1,pMiss))

# try the imputation

library(mice)
md.pattern(editExprs_hq_new)

library(VIM)
aggr_plot <- aggr(editExprs_hq_new, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

colnames(editExprs_hq_new) <- paste0("Library",colnames(editExprs_hq_new))
write.table(editExprs_hq_new,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_preimputedresults_forebrain_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_impute <- mice(editExprs_hq_new,m=5,maxit=30,meth="pmm",seed=500)

editExprs_complete <- complete(editExprs_impute,1)
write.table(editExprs_complete,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_forebrain_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_forebrain_60per.txt", sep = ",", header = T, quote = "")

editExprs_complete[GRIK2_sites,]

densityplot(editExprs_impute)



editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 3, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))


design <- model.matrix(~0+Stage+Sex)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_60per_noADAR.txt", sep = "," ,quote = F, row.names = T, col.names = T)

editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_forebrain_60per.txt", sep = ",", header = T, quote = "")

Organ_ADARexp <- read.table("git_repo/RNAEditing/Figure2/ADARexp_OrganDev.txt", header = T)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 3)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Library, split = "sT"), "[[", 1)
Organ_ADARexp$Organ <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 5)
Organ_ADARexp$Sex <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 7)
Organ_ADARexp$log2TPM <- log2(Organ_ADARexp$TPM)
head(Organ_ADARexp)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 5, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex","ADAR1exp","ADAR2exp"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
        editMeta[sample,"ADAR1exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADAR", "log2TPM"]
        editMeta[sample,"ADAR2exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADARB1", "log2TPM"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))
    ADAR1 <- editMeta$ADAR1exp
    ADAR2 <- editMeta$ADAR2exp


design <- model.matrix(~0+Stage+Sex+ADAR1+ADAR2)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = F, row.names = T, col.names = T)


results_table_sig <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_ADAR1ADAR2_20240202.txt", sep = "," ,quote = "", row.names = 1, header = T)
results_table_sig$Chrom_Pos <- paste0("chr", row.names(results_table_sig))
results_table_sig$Chrom_Pos <- stringr:::str_replace_all(results_table_sig$Chrom_Pos, "_", " ")
results_table_sig$Chrom_Pos <- paste0(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 1), " ", as.numeric(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 2))-1, " ", as.numeric(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 2)))

MB_diffedit_forebrain_hg19 <- read.table("Figure_OrganDevelopment_R1/MB_Forebrain_sites_hg19_20231207.txt")
dim(MB_diffedit_forebrain_hg19)
MB_diffedit_forebrain_hg19$Chrom_Pos <- paste0(MB_diffedit_forebrain_hg19$V1, "_", MB_diffedit_forebrain_hg19$V3)
MB_diffedit_forebrain_hg19$Chrom_Pos <- sapply(strsplit(MB_diffedit_forebrain_hg19$Chrom_Pos, split = "chr"), "[[", 2)
MB_diffedit_forebrain_hg19 <- MB_diffedit_forebrain_hg19[!MB_diffedit_forebrain_hg19$Chrom_Pos %in% c(SNP_sites, ACALAP_sites, MITO_sites),]
dim(MB_diffedit_forebrain_hg19)

length(intersect(row.names(results_table_sig), MB_diffedit_forebrain_hg19$Chrom_Pos))

results_table_sig$Chrom_Pos <- paste0("chr", row.names(results_table_sig))
results_table_sig$Chrom_Pos <- stringr:::str_replace_all(results_table_sig$Chrom_Pos, "_", " ")
results_table_sig$Chrom_Pos <- paste0(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 1), " ", as.numeric(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 2))-1, " ", as.numeric(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 2)))

write.table(results_table_sig$Chrom_Pos, "Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_ADAR1ADAR2_siteshg19_20240202.txt", quote = F, row.names = F, col.names = F)
results_table_hg38 <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_ADAR1ADAR2_siteshg38_20240202.txt")
results_table_hg38$Chrom_Pos <- paste0(results_table_hg38$V1, "_", results_table_hg38$V2)

results_table_sig$Chrom_Pos_hg38 <- results_table_hg38$Chrom_Pos
write.table(results_table_sig, "Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_ADAR1ADAR2_20240202.txt", sep = "," ,quote = F, row.names = T, col.names = T)

results_table_sig1 <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_20240130.txt", sep = "," ,quote = "", row.names = 1, header = T)
results_table_sig2 <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_ADAR1ADAR2_20240202.txt", sep = "," ,quote = "", row.names = 1, header = T)
dim(results_table_sig1) ; dim(results_table_sig2)
length(intersect(row.names(results_table_sig1), row.names(results_table_sig2)))

MB_diffedit_forebrain <- read.table("Figure_OrganDevelopment_R1/MB_DiffEditSites_Forebrain_all.csv", skip = 1, sep = ",", header = 1)
dim(MB_diffedit_forebrain)
MB_diffedit_forebrain_hg19 <- read.table("Figure_OrganDevelopment_R1/MB_Forebrain_sites_hg19_20231207.txt")
MB_diffedit_forebrain_hg19$Chrom_Pos <- paste0(MB_diffedit_forebrain_hg19$V1, "_", MB_diffedit_forebrain_hg19$V3)
MB_diffedit_forebrain_hg19$Chrom_Pos <- sapply(strsplit(MB_diffedit_forebrain_hg19$Chrom_Pos, split = "chr"), "[[", 2)
dim(MB_diffedit_forebrain_hg19)

MB_diffedit_forebrain$Chrom_Pos_hg19 <- MB_diffedit_forebrain_hg19$Chrom_Pos
MB_diffedit_forebrain <- MB_diffedit_forebrain[MB_diffedit_forebrain$adj.P.Val.1 < 0.05,]
dim(MB_diffedit_forebrain)

MB_diffedit_forebrain <- MB_diffedit_forebrain[!MB_diffedit_forebrain$Chrom_Pos_hg19 %in% c(SNP_sites, ACALAP_sites, MITO_sites),]
dim(MB_diffedit_forebrain)

length(intersect(MB_diffedit_forebrain$Chrom_Pos_hg19, row.names(results_table_sig2)))


MB_diffedit_forebrain_hg19 <- read.table("Figure_OrganDevelopment_R1/MB_Forebrain_sites_hg19_20231207.txt")
dim(MB_diffedit_forebrain_hg19)
MB_diffedit_forebrain_hg19 <- MB_diffedit_forebrain_hg19[c(1:4786),]
MB_diffedit_forebrain_hg19$Chrom_Pos <- paste0(MB_diffedit_forebrain_hg19$V1, "_", MB_diffedit_forebrain_hg19$V3)
MB_diffedit_forebrain_hg19$Chrom_Pos <- sapply(strsplit(MB_diffedit_forebrain_hg19$Chrom_Pos, split = "chr"), "[[", 2)
MB_diffedit_forebrain_hg19 <- MB_diffedit_forebrain_hg19[!MB_diffedit_forebrain_hg19$Chrom_Pos %in% c(SNP_sites, ACALAP_sites, MITO_sites),]
dim(MB_diffedit_forebrain_hg19)

length(intersect(row.names(results_table_sig), MB_diffedit_forebrain_hg19$Chrom_Pos))

editExprs_complete <- read.table("Figure_OrganDevelopment_R1/hq_sites_sample_imputedresults_forebrain_v3_20240130.txt", sep = ",", header = T, quote = "")

head(results_table_sig)
table(results_table_sig$ALU)
table(results_table_sig$Region)
table(results_table_sig[results_table_sig$logFC > 0, "ALU"])
results_table_sig[results_table_sig$Region %in% c("exonic"),]
row.names(results_table_sig[results_table_sig$logFC > 0, ])

editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0, ]),]
head(OrganDev_meta)
OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"]
OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"]

rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0, ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])
rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0, ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])

results_table_sig[results_table_sig$Gene == "CYFIP2",]

####### Hindbrain

# try checking site sparsity and then maybe imputation
edit_table_cb_hindbrain <- edit_table_cb[edit_table_cb$Organ == "Hindbrain/Cerebellum",]
dim(edit_table_cb_hindbrain) # 41996360       16
length(unique(edit_table_cb_hindbrain$Chrom_Pos)) # 3043941

#edit_table_cb_hindbrain <- edit_table_cb_hindbrain[edit_table_cb_hindbrain$NonRefEdited >= 3,]
dim(edit_table_cb_hindbrain) # 570500       16
edit_table_cb_hindbrain <- edit_table_cb_hindbrain[edit_table_cb_hindbrain$TotalEditedReads >= 5,]
dim(edit_table_cb_hindbrain) # 41996360       16

edit_table_cb_hindbrain <- edit_table_cb_hindbrain[edit_table_cb_hindbrain$Chrom != "M",]
edit_table_cb_hindbrain <- edit_table_cb_hindbrain[!edit_table_cb_hindbrain$Chrom_Pos %in% SNP_sites,]
edit_table_cb_hindbrain <- edit_table_cb_hindbrain[!grepl("^AC|^AL|^AP", edit_table_cb_hindbrain$Gene),]

dim(edit_table_cb_hindbrain) # 30264149       16
length(unique(edit_table_cb_hindbrain$Chrom_Pos)) # 2169809

sample_cutoff = 0.6
editlevel_cutoff = 0.05


sample_count <- length(table(edit_table_cb_hindbrain$Subject))
sample_count # 59
highquality_sites_coverage <- names(table(edit_table_cb_hindbrain$Chrom_Pos)[which(table(edit_table_cb_hindbrain$Chrom_Pos) >= sample_cutoff*sample_count)])
length(highquality_sites_coverage) # 246287

highquality_sites_minediting_table <- aggregate(edit_table_cb_hindbrain$EditLevel, list(edit_table_cb_hindbrain$Chrom_Pos), FUN=mean)
head(highquality_sites_minediting_table)
highquality_sites_minediting <- highquality_sites_minediting_table[highquality_sites_minediting_table$x >= editlevel_cutoff, "Group.1"]
length(highquality_sites_minediting) # 316331

length(intersect(highquality_sites_coverage, highquality_sites_minediting)) # 11881
highquality_sites_hindbrain <- intersect(highquality_sites_coverage, highquality_sites_minediting)

edit_table_cb_hindbrain_qc <- edit_table_cb_hindbrain[edit_table_cb_hindbrain$Chrom_Pos %in% highquality_sites_hindbrain, ]
dim(edit_table_cb_hindbrain_qc) # 609133
summary(edit_table_cb_hindbrain_qc)
write.table(edit_table_cb_hindbrain_qc,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_hindbrain_qc_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)

OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1.csv', sep = ',', header = T)

head(OrganDev_meta)
OrganDev_meta$Developmental_stage_sum <- NA
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc'), "Developmental_stage_sum"] <- "Early Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc'), "Developmental_stage_sum"] <- "Late Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('newborn','infant','toddler','school','youngTeenager','teenager'), "Developmental_stage_sum"] <- "Newborn - Teenager"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'), "Developmental_stage_sum"] <- "Adult - Senior"

OrganDev_meta$Stage <- NA
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Early Gestation","Late Gestation"), "Stage"] <- "Prenatal"
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Newborn - Teenager","Adult - Senior"), "Stage"] <- "Postnatal"


sample_groups <- OrganDev_meta$library.ID
sites_of_interest <- highquality_sites_hindbrain
editExprs <- data.frame(matrix(NA, nrow = length(sites_of_interest), ncol = length(sample_groups), dimnames = list(sites_of_interest, sample_groups)))
    colnames(editExprs) <- sample_groups
    for (k in 1:nrow(edit_table_cb_hindbrain_qc)){
        editExprs_sites <- edit_table_cb_hindbrain_qc[k, "Chrom_Pos"]
        editExprs_sample <- edit_table_cb_hindbrain_qc[k, "Subject"]
        editExprs[editExprs_sites, as.character(editExprs_sample)] <- edit_table_cb_hindbrain_qc[k, "EditLevel"]
    }
dim(editExprs)
highquality_samples <- colnames(is.na(editExprs))[colSums(is.na(editExprs)) < 0.2*nrow(editExprs)]


table(OrganDev_meta[OrganDev_meta$library.ID %in% highquality_samples, "Developmental_stage_sum"])

editExprs_hq <- editExprs[colnames(editExprs) %in% highquality_samples]
dim(editExprs_hq)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(editExprs_hq,2,pMiss)
max(apply(editExprs_hq,2,pMiss)) # 19.5%
mean(apply(editExprs_hq,2,pMiss)) # 7.37%
apply(editExprs_hq,1,pMiss)
max(apply(editExprs_hq,1,pMiss)) # 29.54%
mean(apply(editExprs_hq,1,pMiss)) # 7.37%

aggr_plot <- aggr(editExprs_hq, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))


# after dropping the low quality sample, clean up the sites one more time
highquality_sites_new <- names(apply(editExprs_hq,1,pMiss))[apply(editExprs_hq,1,pMiss) < (sample_cutoff*100)]
length(highquality_sites_new)
highquality_sites_new <- setdiff(highquality_sites_new, names(which(rowSums(editExprs_hq, na.rm = T) <= editlevel_cutoff)))
length(highquality_sites_new)

editExprs_hq_new <- editExprs_hq[row.names(editExprs_hq) %in% highquality_sites_new,]

apply(editExprs_hq_new,2,pMiss)
max(apply(editExprs_hq_new,2,pMiss))
mean(apply(editExprs_hq_new,2,pMiss))
apply(editExprs_hq_new,1,pMiss)
max(apply(editExprs_hq_new,1,pMiss))
mean(apply(editExprs_hq_new,1,pMiss))

# try the imputation

library(mice)
md.pattern(editExprs_hq_new)

library(VIM)
aggr_plot <- aggr(editExprs_hq_new, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

colnames(editExprs_hq_new) <- paste0("Library",colnames(editExprs_hq_new))
write.table(editExprs_hq_new,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_preimputedresults_hindbrain_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_impute <- mice(editExprs_hq_new,m=5,maxit=30,meth="pmm",seed=500)

editExprs_complete <- complete(editExprs_impute,1)
write.table(editExprs_complete,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_hindbrain_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)

editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_hindbrain_60per.txt", sep = ",", header = T, quote = "")


editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 3, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))


design <- model.matrix(~0+Stage+Sex)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_60per_noADAR.txt", sep = "," ,quote = F, row.names = T, col.names = T)


editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_hindbrain_60per.txt", sep = ",", header = T, quote = "")

densityplot(editExprs_impute)

Organ_ADARexp <- read.table("git_repo/RNAEditing/Figure2/ADARexp_OrganDev.txt", header = T)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 3)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Library, split = "sT"), "[[", 1)
Organ_ADARexp$Organ <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 5)
Organ_ADARexp$Sex <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 7)
Organ_ADARexp$log2TPM <- log2(Organ_ADARexp$TPM)
head(Organ_ADARexp)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 5, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex","ADAR1exp","ADAR2exp"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
        editMeta[sample,"ADAR1exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADAR", "log2TPM"]
        editMeta[sample,"ADAR2exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADARB1", "log2TPM"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))
    ADAR1 <- editMeta$ADAR1exp
    ADAR2 <- editMeta$ADAR2exp

design <- model.matrix(~0+Stage+Sex+ADAR1+ADAR2)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = F, row.names = T, col.names = T)


write.table(results_table_sig, "Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_ADAR1ADAR2_20240202.txt", sep = "," ,quote = F, row.names = T, col.names = T)
results_table_sig <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_ADAR1ADAR2_20240202.txt", sep = "," ,quote = "", row.names = 1, header = T)
results_table_sig$Chrom_Pos <- paste0("chr", row.names(results_table_sig))
results_table_sig$Chrom_Pos <- stringr:::str_replace_all(results_table_sig$Chrom_Pos, "_", " ")
results_table_sig$Chrom_Pos <- paste0(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 1), " ", as.numeric(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 2))-1, " ", as.numeric(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 2)))

#MB_diffedit_hindbrain <- read.table("Figure_OrganDevelopment_R1/MB_DiffEditSites_Hindbrain-temporaleffect.txt")
#MB_diffedit_hindbrain$V2 <- as.numeric(sapply(strsplit(MB_diffedit_hindbrain$V1, split = "_"),"[[", 2))-1
#MB_diffedit_hindbrain$V3 <- sapply(strsplit(MB_diffedit_hindbrain$V1, split = "_"),"[[", 2)
#MB_diffedit_hindbrain$V1 <- sapply(strsplit(MB_diffedit_hindbrain$V1, split = "_"),"[[", 1)

#write.table(MB_diffedit_hindbrain, "Figure_OrganDevelopment_R1/MB_hindbrain_sites_hg38_20231207.txt", quote = F, row.names = F, col.names = F)

MB_diffedit_hindbrain <- read.table("Figure_OrganDevelopment_R1/MB_DiffEditSites_Hindbrain_all.csv", skip = 1, sep = ",", header = 1)
dim(MB_diffedit_hindbrain)
MB_diffedit_hindbrain_hg19 <- read.table("Figure_OrganDevelopment_R1/MB_Hindbrain_sites_hg19_20231207.txt")
MB_diffedit_hindbrain_hg19$Chrom_Pos <- paste0(MB_diffedit_hindbrain_hg19$V1, "_", MB_diffedit_hindbrain_hg19$V3)
MB_diffedit_hindbrain_hg19$Chrom_Pos <- sapply(strsplit(MB_diffedit_hindbrain_hg19$Chrom_Pos, split = "chr"), "[[", 2)
dim(MB_diffedit_hindbrain_hg19)

MB_diffedit_hindbrain$Chrom_Pos_hg19 <- MB_diffedit_hindbrain_hg19$Chrom_Pos
MB_diffedit_hindbrain <- MB_diffedit_hindbrain[MB_diffedit_hindbrain$adj.P.Val.1 < 0.05,]
dim(MB_diffedit_hindbrain)

MB_diffedit_hindbrain <- MB_diffedit_hindbrain[!MB_diffedit_hindbrain$Chrom_Pos_hg19 %in% c(SNP_sites, ACALAP_sites, MITO_sites),]
dim(MB_diffedit_hindbrain)

length(intersect(MB_diffedit_hindbrain$Chrom_Pos_hg19, row.names(results_table_sig2)))



MB_diffedit_hindbrain_hg19 <- read.table("Figure_OrganDevelopment_R1/MB_Hindbrain_sites_hg19_20231207.txt")
dim(MB_diffedit_hindbrain_hg19)
MB_diffedit_hindbrain_hg19$Chrom_Pos <- paste0(MB_diffedit_hindbrain_hg19$V1, "_", MB_diffedit_hindbrain_hg19$V3)
MB_diffedit_hindbrain_hg19$Chrom_Pos <- sapply(strsplit(MB_diffedit_hindbrain_hg19$Chrom_Pos, split = "chr"), "[[", 2)
MB_diffedit_hindbrain_hg19 <- MB_diffedit_hindbrain_hg19[!MB_diffedit_hindbrain_hg19$Chrom_Pos %in% c(SNP_sites, ACALAP_sites, MITO_sites),]
dim(MB_diffedit_hindbrain_hg19)

length(intersect(row.names(results_table_sig), MB_diffedit_hindbrain_hg19$Chrom_Pos))

results_table_sig$Chrom_Pos <- paste0("chr", row.names(results_table_sig))
results_table_sig$Chrom_Pos <- stringr:::str_replace_all(results_table_sig$Chrom_Pos, "_", " ")
results_table_sig$Chrom_Pos <- paste0(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 1), " ", as.numeric(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 2))-1, " ", as.numeric(sapply(strsplit(results_table_sig$Chrom_Pos, split = "[ ]"), "[[", 2)))

write.table(results_table_sig$Chrom_Pos, "Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_ADAR1ADAR2_siteshg19_20240202.txt", quote = F, row.names = F, col.names = F)
results_table_hg38 <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_ADAR1ADAR2_siteshg38_20240202.txt")
results_table_hg38$Chrom_Pos <- paste0(results_table_hg38$V1, "_", results_table_hg38$V2)

results_table_sig$Chrom_Pos_hg38 <- results_table_hg38$Chrom_Pos
write.table(results_table_sig, "Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_ADAR1ADAR2_20240202.txt", sep = "," ,quote = F, row.names = T, col.names = T)


results_table_sig <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240131.txt", sep = "," ,quote = "", row.names = 1, header = T)


results_table_sig1 <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240131.txt", sep = "," ,quote = "", row.names = 1, header = T)
results_table_sig2 <- read.table("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_ADAR1ADAR2_20240202.txt", sep = "," ,quote = "", row.names = 1, header = T)
dim(results_table_sig1) ; dim(results_table_sig2)
length(intersect(row.names(results_table_sig1), row.names(results_table_sig2)))

MB_diffedit_hindbrain_hg19 <- read.table("Figure_OrganDevelopment_R1/MB_Hindbrain_sites_hg19_20231207.txt")
dim(MB_diffedit_hindbrain_hg19)
MB_diffedit_hindbrain_hg19 <- MB_diffedit_hindbrain_hg19[c(1:6662),] # take significant ones
MB_diffedit_hindbrain_hg19$Chrom_Pos <- paste0(MB_diffedit_hindbrain_hg19$V1, "_", MB_diffedit_hindbrain_hg19$V3)
MB_diffedit_hindbrain_hg19$Chrom_Pos <- sapply(strsplit(MB_diffedit_hindbrain_hg19$Chrom_Pos, split = "chr"), "[[", 2)
MB_diffedit_hindbrain_hg19 <- MB_diffedit_hindbrain_hg19[!MB_diffedit_hindbrain_hg19$Chrom_Pos %in% c(SNP_sites, ACALAP_sites, MITO_sites),]
dim(MB_diffedit_hindbrain_hg19)

length(intersect(row.names(results_table_sig), MB_diffedit_hindbrain_hg19$Chrom_Pos))

library(ggplot2)
# all postnatal biased cites
df2 <- data.frame(group=rep(c("Prenatal", "Postnatal"), each=3),
                alu=rep(names(table(results_table_sig[row.names(results_table_sig[results_table_sig$logFC > 0, ]),"ALU"])),2),
                meaneditlevel=c(mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "ALU", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "NONREP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "REP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "ALU", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "NONREP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "REP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]))))
df2$group <- factor(df2$group, level = c("Prenatal", "Postnatal"))
pdf("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240130_postnatalbiasedSites_barplot.pdf")
ggplot(data=df2, aes(x=alu, y=meaneditlevel, fill=group)) +
geom_bar(stat="identity", position=position_dodge())
dev.off()
ggplot(data=df2, aes(x=alu, y=meaneditlevel, fill=group)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=meaneditlevel), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

length(results_table_sig[row.names(results_table_sig[results_table_sig$logFC > 0, ]),"ALU"])
df <- data.frame(group=rep(c("Prenatal", "Postnatal"), each=7803),
                alu=rep(results_table_sig[row.names(results_table_sig[results_table_sig$logFC > 0, ]),"ALU"],2),
                meaneditlevel=c(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "ALU", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "NONREP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "REP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "ALU", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "NONREP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$ALU == "REP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])))
df$group <- factor(df$group, level = c("Prenatal", "Postnatal"))
pdf("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240130_postnatalbiasedSites_boxplot.pdf")
ggplot(df, aes(x=alu, y=meaneditlevel, fill=group)) + 
  geom_boxplot(position=position_dodge(1))
dev.off()

df <- data.frame(group=rep(c("Prenatal", "Postnatal"), each=11),
                alu=rep(names(table(results_table_sig[row.names(results_table_sig[results_table_sig$logFC > 0, ]),"Region"])),2),
                meaneditlevel=c(mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "downstream", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "exonic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "intergenic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "intronic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "ncRNA_exonic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "ncRNA_intronic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "ncRNA_splicing", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "splicing", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "upstream", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "UTR3", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "UTR5", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),

                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "downstream", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "exonic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "intergenic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "intronic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "ncRNA_exonic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "ncRNA_intronic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "ncRNA_splicing", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "splicing", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "upstream", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "UTR3", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC > 0 & results_table_sig$Region == "UTR5", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]))))
df$group <- factor(df$group, level = c("Prenatal", "Postnatal"))
pdf("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240130_postnatalbiasedSites_region_barplot.pdf")
ggplot(df, aes(x=alu, y=meaneditlevel, fill=group)) + 
  geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# all prenatal biased cites
df2 <- data.frame(group=rep(c("Prenatal", "Postnatal"), each=2),
                alu=rep(names(table(results_table_sig[row.names(results_table_sig[results_table_sig$logFC < 0, ]),"ALU"])),2),
                meaneditlevel=c(mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "ALU", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "NONREP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                #mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "REP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "ALU", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "NONREP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]))))
                #mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "REP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]))))
df2$group <- factor(df2$group, level = c("Prenatal", "Postnatal"))
pdf("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240130_prenatalbiasedSites_barplot.pdf")
ggplot(data=df2, aes(x=alu, y=meaneditlevel, fill=group)) +
geom_bar(stat="identity", position=position_dodge())
dev.off()
ggplot(data=df2, aes(x=alu, y=meaneditlevel, fill=group)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=meaneditlevel), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

df <- data.frame(group=rep(c("Prenatal", "Postnatal"), each=119),
                alu=rep(results_table_sig[row.names(results_table_sig[results_table_sig$logFC < 0, ]),"ALU"],2),
                meaneditlevel=c(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "ALU", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "NONREP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                #rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "REP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "ALU", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "NONREP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])))
                #rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$ALU == "REP", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]))
df$group <- factor(df$group, level = c("Prenatal", "Postnatal"))
pdf("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240130_prenatalbiasedSites_boxplot.pdf")
ggplot(df, aes(x=alu, y=meaneditlevel, fill=group)) + 
  geom_boxplot(position=position_dodge(1))
dev.off()

editExprs_complete <- read.table("Figure_OrganDevelopment_R1/hq_sites_sample_imputedresults_hindbrain_v3_20240131.txt", sep = ",", header = T, quote = "")

df <- data.frame(group=rep(c("Prenatal", "Postnatal"), each=11),
                alu=rep(names(table(results_table_sig[row.names(results_table_sig[results_table_sig$logFC > 0, ]),"Region"])),2),
                meaneditlevel=c(mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "downstream", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "exonic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "intergenic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "intronic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "ncRNA_exonic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "ncRNA_intronic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "ncRNA_splicing", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "splicing", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "upstream", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "UTR3", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "UTR5", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])])),

                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "downstream", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "exonic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "intergenic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "intronic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "ncRNA_exonic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "ncRNA_intronic", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "ncRNA_splicing", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "splicing", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "upstream", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "UTR3", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])),
                mean(rowMeans(editExprs_complete[row.names(results_table_sig[results_table_sig$logFC < 0 & results_table_sig$Region == "UTR5", ]), colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]))))
df$group <- factor(df$group, level = c("Prenatal", "Postnatal"))
pdf("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240130_prenatalbiasedSites_region_barplot.pdf")
ggplot(df, aes(x=alu, y=meaneditlevel, fill=group)) + 
  geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




# check for GRIK2 site editing levels
df <- data.frame(group=rep(c("Prenatal", "Postnatal"), each=4),
                site=rep(GRIK2_sites,2),
                meaneditlevel=c(rowMeans(editExprs_complete[GRIK2_sites[1], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIK2_sites[2], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIK2_sites[3], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIK2_sites[4], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIK2_sites[1], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIK2_sites[2], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIK2_sites[3], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIK2_sites[4], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])))
df$group <- factor(df$group, level = c("Prenatal", "Postnatal"))
pdf("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240130_GRIK2sites_barplot.pdf")
ggplot(df, aes(x=site, y=meaneditlevel, fill=group)) + 
    geom_bar(stat="identity", position=position_dodge())
dev.off()

# check for GRIA2 site editing levels
GRIA2_sites <- row.names(results_table_sig[results_table_sig$Gene == "GRIA2", ])
df <- data.frame(group=rep(c("Prenatal", "Postnatal"), each=4),
                site=rep(GRIA2_sites,2),
                meaneditlevel=c(rowMeans(editExprs_complete[GRIA2_sites[1], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIA2_sites[2], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIA2_sites[3], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIA2_sites[4], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Prenatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIA2_sites[1], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIA2_sites[2], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIA2_sites[3], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])]),
                rowMeans(editExprs_complete[GRIA2_sites[4], colnames(editExprs_complete) %in% paste0("Library",OrganDev_meta[OrganDev_meta$Stage == "Postnatal", "library.ID"])])))
df$group <- factor(df$group, level = c("Prenatal", "Postnatal"))
pdf("Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_20240130_GRIA2sites_barplot.pdf")
ggplot(df, aes(x=site, y=meaneditlevel, fill=group)) + 
    geom_bar(stat="identity", position=position_dodge())
dev.off()


####### Heart

# try checking site sparsity and then maybe imputation
edit_table_cb_heart <- edit_table_cb[edit_table_cb$Organ == "Heart",]
dim(edit_table_cb_heart) # 21717163       16
length(unique(edit_table_cb_heart$Chrom_Pos)) # 2136633

#edit_table_cb_heart <- edit_table_cb_heart[edit_table_cb_heart$NonRefEdited >= 3,]
#dim(edit_table_cb_heart) # 570500       16
edit_table_cb_heart <- edit_table_cb_heart[edit_table_cb_heart$TotalEditedReads >= 5,]
dim(edit_table_cb_heart) # 21717163       16

edit_table_cb_heart <- edit_table_cb_heart[edit_table_cb_heart$Chrom != "M",]
edit_table_cb_heart <- edit_table_cb_heart[!edit_table_cb_heart$Chrom_Pos %in% SNP_sites,]
edit_table_cb_heart <- edit_table_cb_heart[!grepl("^AC|^AL|^AP", edit_table_cb_heart$Gene),]

dim(edit_table_cb_heart) # 15940752       16
length(unique(edit_table_cb_heart$Chrom_Pos)) # 1564447

sample_cutoff = 0.6
editlevel_cutoff = 0.05



sample_count <- length(table(edit_table_cb_heart$Subject))
sample_count # 50
highquality_sites_coverage <- names(table(edit_table_cb_heart$Chrom_Pos)[which(table(edit_table_cb_heart$Chrom_Pos) >= (sample_cutoff*sample_count))])
length(highquality_sites_coverage) # 129851

highquality_sites_minediting_table <- aggregate(edit_table_cb_heart$EditLevel, list(edit_table_cb_heart$Chrom_Pos), FUN=mean)
head(highquality_sites_minediting_table)
highquality_sites_minediting <- highquality_sites_minediting_table[highquality_sites_minediting_table$x >= editlevel_cutoff, "Group.1"]
length(highquality_sites_minediting) # 101572

length(intersect(highquality_sites_coverage, highquality_sites_minediting)) # 2055
highquality_sites_heart <- intersect(highquality_sites_coverage, highquality_sites_minediting)

edit_table_cb_heart_qc <- edit_table_cb_heart[edit_table_cb_heart$Chrom_Pos %in% highquality_sites_heart, ]
dim(edit_table_cb_heart_qc) # 85436
summary(edit_table_cb_heart_qc)
write.table(edit_table_cb_heart_qc,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_heart_qc_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)

OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1.csv', sep = ',', header = T)

head(OrganDev_meta)
OrganDev_meta$Developmental_stage_sum <- NA
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc'), "Developmental_stage_sum"] <- "Early Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc'), "Developmental_stage_sum"] <- "Late Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('newborn','infant','toddler','school','youngTeenager','teenager'), "Developmental_stage_sum"] <- "Newborn - Teenager"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'), "Developmental_stage_sum"] <- "Adult - Senior"

OrganDev_meta$Stage <- NA
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Early Gestation","Late Gestation"), "Stage"] <- "Prenatal"
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Newborn - Teenager","Adult - Senior"), "Stage"] <- "Postnatal"


sample_groups <- OrganDev_meta$library.ID
sites_of_interest <- highquality_sites_heart
editExprs <- data.frame(matrix(NA, nrow = length(sites_of_interest), ncol = length(sample_groups), dimnames = list(sites_of_interest, sample_groups)))
    colnames(editExprs) <- sample_groups
    for (k in 1:nrow(edit_table_cb_heart_qc)){
        editExprs_sites <- edit_table_cb_heart_qc[k, "Chrom_Pos"]
        editExprs_sample <- edit_table_cb_heart_qc[k, "Subject"]
        editExprs[editExprs_sites, as.character(editExprs_sample)] <- edit_table_cb_heart_qc[k, "EditLevel"]
    }
dim(editExprs)
highquality_samples <- colnames(is.na(editExprs))[colSums(is.na(editExprs)) < 0.2*nrow(editExprs)]


table(OrganDev_meta[OrganDev_meta$library.ID %in% highquality_samples, "Developmental_stage_sum"])

editExprs_hq <- editExprs[colnames(editExprs) %in% highquality_samples]
dim(editExprs_hq)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(editExprs_hq,2,pMiss)
max(apply(editExprs_hq,2,pMiss)) # 12.8%
mean(apply(editExprs_hq,2,pMiss)) # 3.61%
apply(editExprs_hq,1,pMiss)
max(apply(editExprs_hq,1,pMiss)) # 34.29%
mean(apply(editExprs_hq,1,pMiss)) # 3.61%

aggr_plot <- aggr(editExprs_hq, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))


# after dropping the low quality sample, clean up the sites one more time
highquality_sites_new <- names(apply(editExprs_hq,1,pMiss))[apply(editExprs_hq,1,pMiss) < (sample_cutoff*100)]
length(highquality_sites_new)
highquality_sites_new <- setdiff(highquality_sites_new, names(which(rowSums(editExprs_hq, na.rm = T) <= editlevel_cutoff)))
length(highquality_sites_new)

editExprs_hq_new <- editExprs_hq[row.names(editExprs_hq) %in% highquality_sites_new,]

apply(editExprs_hq_new,2,pMiss)
max(apply(editExprs_hq_new,2,pMiss))
mean(apply(editExprs_hq_new,2,pMiss))
apply(editExprs_hq_new,1,pMiss)
max(apply(editExprs_hq_new,1,pMiss))
mean(apply(editExprs_hq_new,1,pMiss))

# try the imputation

library(mice)
md.pattern(editExprs_hq_new)

library(VIM)
aggr_plot <- aggr(editExprs_hq_new, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

colnames(editExprs_hq_new) <- paste0("Library",colnames(editExprs_hq_new))
write.table(editExprs_hq_new,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_preimputedresults_heart_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_impute <- mice(editExprs_hq_new,m=5,maxit=30,meth="pmm",seed=500)

editExprs_complete <- complete(editExprs_impute,1)
write.table(editExprs_complete,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_heart_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_heart_60per.txt", sep = ",", header = T, quote = "")

densityplot(editExprs_impute)

Organ_ADARexp <- read.table("git_repo/RNAEditing/Figure2/ADARexp_OrganDev.txt", header = T)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 3)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Library, split = "sT"), "[[", 1)
Organ_ADARexp$Organ <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 5)
Organ_ADARexp$Sex <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 7)
Organ_ADARexp$log2TPM <- log2(Organ_ADARexp$TPM)
head(Organ_ADARexp)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 3, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Developmental_stage_sum"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Developmental_stage_sum == sample_devstage, "Sex"]
    }

table(editMeta$Stage)
editMeta[editMeta$Stage == "Newborn - Teenager", "Subject"]

editExprs_complete <- editExprs_complete[, !colnames(editExprs_complete) %in% editMeta[editMeta$Stage == "Newborn - Teenager", "Subject"]]
editMeta <- editMeta[!row.names(editMeta) %in% editMeta[editMeta$Stage == "Newborn - Teenager", "Subject"], ]

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    editMeta$Stage <- stringr:::str_replace_all(editMeta$Stage, " ", "")
    Stage <- factor(editMeta$Stage, levels = c("EarlyGestation","LateGestation"))

design <- model.matrix(~0+Stage+Sex)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StageLateGestation - StageEarlyGestation),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_heart_earlygestationVSlategestation_60per_noADAR.txt", sep = "," ,quote = F, row.names = T, col.names = T)

editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_heart_60per.txt", sep = ",", header = T, quote = "")


editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 5, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex","ADAR1exp","ADAR2exp"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Developmental_stage_sum"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Developmental_stage_sum == sample_devstage, "Sex"]
        editMeta[sample,"ADAR1exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADAR", "log2TPM"]
        editMeta[sample,"ADAR2exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADARB1", "log2TPM"]
    }

editExprs_complete <- editExprs_complete[, !colnames(editExprs_complete) %in% editMeta[editMeta$Stage == "Newborn - Teenager", "Subject"]]
editMeta <- editMeta[!row.names(editMeta) %in% editMeta[editMeta$Stage == "Newborn - Teenager", "Subject"], ]

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    editMeta$Stage <- stringr:::str_replace_all(editMeta$Stage, " ", "")
    Stage <- factor(editMeta$Stage, levels = c("EarlyGestation","LateGestation"))
    ADAR1 <- editMeta$ADAR1exp
    ADAR2 <- editMeta$ADAR2exp

design <- model.matrix(~0+Stage+Sex+ADAR1+ADAR2)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StageLateGestation - StageEarlyGestation),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_heart_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = F, row.names = T, col.names = T)


####### Liver

# try checking site sparsity and then maybe imputation
edit_table_cb_liver <- edit_table_cb[edit_table_cb$Organ == "Liver",]
dim(edit_table_cb_liver) # 14897110       16
length(unique(edit_table_cb_liver$Chrom_Pos)) # 1593780

#edit_table_cb_liver <- edit_table_cb_liver[edit_table_cb_liver$NonRefEdited >= 3,]
#dim(edit_table_cb_liver) # 570500       16
edit_table_cb_liver <- edit_table_cb_liver[edit_table_cb_liver$TotalEditedReads >= 5,]
dim(edit_table_cb_liver) # 14897110       16

edit_table_cb_liver <- edit_table_cb_liver[edit_table_cb_liver$Chrom != "M",]
edit_table_cb_liver <- edit_table_cb_liver[!edit_table_cb_liver$Chrom_Pos %in% SNP_sites,]
edit_table_cb_liver <- edit_table_cb_liver[!grepl("^AC|^AL|^AP", edit_table_cb_liver$Gene),]

dim(edit_table_cb_liver) # 11022945       16
length(unique(edit_table_cb_liver$Chrom_Pos)) # 1171965


sample_cutoff = 0.6
editlevel_cutoff = 0.05

sample_count <- length(table(edit_table_cb_liver$Subject))
sample_count # 50
highquality_sites_coverage <- names(table(edit_table_cb_liver$Chrom_Pos)[which(table(edit_table_cb_liver$Chrom_Pos) > (sample_cutoff*sample_count))])
length(highquality_sites_coverage) # 102479

highquality_sites_minediting_table <- aggregate(edit_table_cb_liver$EditLevel, list(edit_table_cb_liver$Chrom_Pos), FUN=mean)
head(highquality_sites_minediting_table)
highquality_sites_minediting <- highquality_sites_minediting_table[highquality_sites_minediting_table$x >= editlevel_cutoff, "Group.1"]
length(highquality_sites_minediting) # 106370

length(intersect(highquality_sites_coverage, highquality_sites_minediting)) # 2393
highquality_sites_liver <- intersect(highquality_sites_coverage, highquality_sites_minediting)

edit_table_cb_liver_qc <- edit_table_cb_liver[edit_table_cb_liver$Chrom_Pos %in% highquality_sites_liver, ]
dim(edit_table_cb_liver_qc) # 102663
summary(edit_table_cb_liver_qc)
write.table(edit_table_cb_forebrain_qc,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_liver_qc_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)


OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1.csv', sep = ',', header = T)

head(OrganDev_meta)
OrganDev_meta$Developmental_stage_sum <- NA
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc'), "Developmental_stage_sum"] <- "Early Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc'), "Developmental_stage_sum"] <- "Late Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('newborn','infant','toddler','school','youngTeenager','teenager'), "Developmental_stage_sum"] <- "Newborn - Teenager"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'), "Developmental_stage_sum"] <- "Adult - Senior"

OrganDev_meta$Stage <- NA
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Early Gestation","Late Gestation"), "Stage"] <- "Prenatal"
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Newborn - Teenager","Adult - Senior"), "Stage"] <- "Postnatal"


sample_groups <- OrganDev_meta$library.ID
sites_of_interest <- highquality_sites_liver
editExprs <- data.frame(matrix(NA, nrow = length(sites_of_interest), ncol = length(sample_groups), dimnames = list(sites_of_interest, sample_groups)))
    colnames(editExprs) <- sample_groups
    for (k in 1:nrow(edit_table_cb_liver_qc)){
        editExprs_sites <- edit_table_cb_liver_qc[k, "Chrom_Pos"]
        editExprs_sample <- edit_table_cb_liver_qc[k, "Subject"]
        editExprs[editExprs_sites, as.character(editExprs_sample)] <- edit_table_cb_liver_qc[k, "EditLevel"]
    }
dim(editExprs)
highquality_samples <- colnames(is.na(editExprs))[colSums(is.na(editExprs)) < 0.2*nrow(editExprs)]


table(OrganDev_meta[OrganDev_meta$library.ID %in% highquality_samples, "Developmental_stage_sum"])

editExprs_hq <- editExprs[colnames(editExprs) %in% highquality_samples]
dim(editExprs_hq)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(editExprs_hq,2,pMiss)
max(apply(editExprs_hq,2,pMiss)) # 17.5%
mean(apply(editExprs_hq,2,pMiss)) # 8.5%
apply(editExprs_hq,1,pMiss)
max(apply(editExprs_hq,1,pMiss)) # 28.21%
mean(apply(editExprs_hq,1,pMiss)) # 7.37%

aggr_plot <- aggr(editExprs_hq, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))


# after dropping the low quality sample, clean up the sites one more time
highquality_sites_new <- names(apply(editExprs_hq,1,pMiss))[apply(editExprs_hq,1,pMiss) < (sample_cutoff*100)]
length(highquality_sites_new)
highquality_sites_new <- setdiff(highquality_sites_new, names(which(rowSums(editExprs_hq, na.rm = T) <= editlevel_cutoff)))
length(highquality_sites_new)

editExprs_hq_new <- editExprs_hq[row.names(editExprs_hq) %in% highquality_sites_new,]

apply(editExprs_hq_new,2,pMiss)
max(apply(editExprs_hq_new,2,pMiss))
mean(apply(editExprs_hq_new,2,pMiss))
apply(editExprs_hq_new,1,pMiss)
max(apply(editExprs_hq_new,1,pMiss))
mean(apply(editExprs_hq_new,1,pMiss))

# try the imputation

library(mice)
md.pattern(editExprs_hq_new)

library(VIM)
aggr_plot <- aggr(editExprs_hq_new, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

colnames(editExprs_hq_new) <- paste0("Library",colnames(editExprs_hq_new))
write.table(editExprs_hq_new,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_preimputedresults_liver_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_impute <- mice(editExprs_hq_new,m=5,maxit=30,meth="pmm",seed=500)

editExprs_complete <- complete(editExprs_impute,1)
write.table(editExprs_complete,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_liver_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)

editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_liver_60per.txt", sep = ",", header = T, quote = "")

densityplot(editExprs_impute)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 3, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))

design <- model.matrix(~0+Stage+Sex)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_liver_prenatalVSpostnatal_60per_noADAR.txt", sep = "," ,quote = F, row.names = T, col.names = T)


editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_liver_60per.txt", sep = ",", header = T, quote = "")

Organ_ADARexp <- read.table("git_repo/RNAEditing/Figure2/ADARexp_OrganDev.txt", header = T)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 3)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Library, split = "sT"), "[[", 1)
Organ_ADARexp$Organ <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 5)
Organ_ADARexp$Sex <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 7)
Organ_ADARexp$log2TPM <- log2(Organ_ADARexp$TPM)
head(Organ_ADARexp)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 5, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex","ADAR1exp","ADAR2exp"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
        editMeta[sample,"ADAR1exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADAR", "log2TPM"]
        editMeta[sample,"ADAR2exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADARB1", "log2TPM"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))
    ADAR1 <- editMeta$ADAR1exp
    ADAR2 <- editMeta$ADAR2exp


design <- model.matrix(~0+Stage+Sex+ADAR1+ADAR2)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_liver_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = F, row.names = T, col.names = T)


editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 3, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Developmental_stage_sum"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Developmental_stage_sum == sample_devstage, "Sex"]
    }

editMeta[editMeta$Stage == "Newborn - Teenager", "Subject"]

editExprs_complete <- editExprs_complete[, !colnames(editExprs_complete) %in% editMeta[editMeta$Stage %in% c("Newborn - Teenager", "Adult - Senior"), "Subject"]]
editMeta <- editMeta[!row.names(editMeta) %in% editMeta[editMeta$Stage %in% c("Newborn - Teenager", "Adult - Senior"), "Subject"], ]

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    editMeta$Stage <- stringr:::str_replace_all(editMeta$Stage, " ", "")
    Stage <- factor(editMeta$Stage, levels = c("EarlyGestation","LateGestation"))

design <- model.matrix(~0+Stage+Sex)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StageLateGestation - StageEarlyGestation),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig) # -> empty table
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
#write.table(results_table_sig, "Figure_OrganDevelopment_R1/SigDiffEditSites_scriptv5_heart_earlygestationVSlategestation_20240110.txt", sep = "," ,quote = F, row.names = T, col.names = T)


####### Testis

# try checking site sparsity and then maybe imputation
edit_table_cb_testis <- edit_table_cb[edit_table_cb$Organ == "Testis",]
dim(edit_table_cb_testis) # 19548287       16
length(unique(edit_table_cb_testis$Chrom_Pos)) # 2279041

#edit_table_cb_testis <- edit_table_cb_testis[edit_table_cb_testis$NonRefEdited >= 3,]
#dim(edit_table_cb_testis) # 570500       16
edit_table_cb_testis <- edit_table_cb_testis[edit_table_cb_testis$TotalEditedReads >= 5,]
dim(edit_table_cb_testis) # 19548287       16

edit_table_cb_testis <- edit_table_cb_testis[edit_table_cb_testis$Chrom != "M",]
edit_table_cb_testis <- edit_table_cb_testis[!edit_table_cb_testis$Chrom_Pos %in% SNP_sites,]
edit_table_cb_testis <- edit_table_cb_testis[!grepl("^AC|^AL|^AP", edit_table_cb_testis$Gene),]

dim(edit_table_cb_testis) # 14143080       16
length(unique(edit_table_cb_testis$Chrom_Pos)) # 1630575

sample_cutoff = 0.6
editlevel_cutoff = 0.05


sample_count <- length(table(edit_table_cb_testis$Subject))
sample_count # 41
highquality_sites_coverage <- names(table(edit_table_cb_testis$Chrom_Pos)[which(table(edit_table_cb_testis$Chrom_Pos) > (sample_cutoff*sample_count))])
length(highquality_sites_coverage) # 152077

highquality_sites_minediting_table <- aggregate(edit_table_cb_testis$EditLevel, list(edit_table_cb_testis$Chrom_Pos), FUN=mean)
head(highquality_sites_minediting_table)
highquality_sites_minediting <- highquality_sites_minediting_table[highquality_sites_minediting_table$x >= editlevel_cutoff, "Group.1"]
length(highquality_sites_minediting) # 155885

length(intersect(highquality_sites_coverage, highquality_sites_minediting)) # 4486
highquality_sites_testis <- intersect(highquality_sites_coverage, highquality_sites_minediting)

edit_table_cb_testis_qc <- edit_table_cb_testis[edit_table_cb_testis$Chrom_Pos %in% highquality_sites_testis, ]
dim(edit_table_cb_testis_qc) # 159901
summary(edit_table_cb_testis_qc)
write.table(edit_table_cb_testis_qc,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_testis_qc_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)


OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1.csv', sep = ',', header = T)

head(OrganDev_meta)
OrganDev_meta$Developmental_stage_sum <- NA
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc'), "Developmental_stage_sum"] <- "Early Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc'), "Developmental_stage_sum"] <- "Late Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('newborn','infant','toddler','school','youngTeenager','teenager'), "Developmental_stage_sum"] <- "Newborn - Teenager"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'), "Developmental_stage_sum"] <- "Adult - Senior"

OrganDev_meta$Stage <- NA
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Early Gestation","Late Gestation"), "Stage"] <- "Prenatal"
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Newborn - Teenager","Adult - Senior"), "Stage"] <- "Postnatal"


sample_groups <- OrganDev_meta$library.ID
sites_of_interest <- highquality_sites_testis
editExprs <- data.frame(matrix(NA, nrow = length(sites_of_interest), ncol = length(sample_groups), dimnames = list(sites_of_interest, sample_groups)))
    colnames(editExprs) <- sample_groups
    for (k in 1:nrow(edit_table_cb_testis_qc)){
        editExprs_sites <- edit_table_cb_testis_qc[k, "Chrom_Pos"]
        editExprs_sample <- edit_table_cb_testis_qc[k, "Subject"]
        editExprs[editExprs_sites, as.character(editExprs_sample)] <- edit_table_cb_testis_qc[k, "EditLevel"]
    }
dim(editExprs)
highquality_samples <- colnames(is.na(editExprs))[colSums(is.na(editExprs)) < 0.2*nrow(editExprs)]


table(OrganDev_meta[OrganDev_meta$library.ID %in% highquality_samples, "Developmental_stage_sum"])

editExprs_hq <- editExprs[colnames(editExprs) %in% highquality_samples]
dim(editExprs_hq)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(editExprs_hq,2,pMiss)
max(apply(editExprs_hq,2,pMiss)) # 17.8%
mean(apply(editExprs_hq,2,pMiss)) # 7.9%
apply(editExprs_hq,1,pMiss)
max(apply(editExprs_hq,1,pMiss)) # 34.4%
mean(apply(editExprs_hq,1,pMiss)) # 7.9%

aggr_plot <- aggr(editExprs_hq, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))


# after dropping the low quality sample, clean up the sites one more time
highquality_sites_new <- names(apply(editExprs_hq,1,pMiss))[apply(editExprs_hq,1,pMiss) < (sample_cutoff*100)]
length(highquality_sites_new)
highquality_sites_new <- setdiff(highquality_sites_new, names(which(rowSums(editExprs_hq, na.rm = T) <= editlevel_cutoff)))
length(highquality_sites_new)

editExprs_hq_new <- editExprs_hq[row.names(editExprs_hq) %in% highquality_sites_new,]

apply(editExprs_hq_new,2,pMiss)
max(apply(editExprs_hq_new,2,pMiss))
mean(apply(editExprs_hq_new,2,pMiss))
apply(editExprs_hq_new,1,pMiss)
max(apply(editExprs_hq_new,1,pMiss))
mean(apply(editExprs_hq_new,1,pMiss))

# try the imputation

library(mice)
md.pattern(editExprs_hq_new)

library(VIM)
aggr_plot <- aggr(editExprs_hq_new, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

colnames(editExprs_hq_new) <- paste0("Library",colnames(editExprs_hq_new))
write.table(editExprs_hq_new,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_preimputedresults_testis_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_impute <- mice(editExprs_hq_new,m=5,maxit=30,meth="pmm",seed=500)

editExprs_complete <- complete(editExprs_impute,1)
write.table(editExprs_complete,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_testis_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_testis_60per.txt", sep = ",", header = T, quote = "")

densityplot(editExprs_impute)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 3, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))

design <- model.matrix(~0+Stage)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_testis_prenatalVSpostnatal_60per_noADAR.txt", sep = "," ,quote = F, row.names = T, col.names = T)

editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_testis_60per.txt", sep = ",", header = T, quote = "")

Organ_ADARexp <- read.table("git_repo/RNAEditing/Figure2/ADARexp_OrganDev.txt", header = T)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 3)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Library, split = "sT"), "[[", 1)
Organ_ADARexp$Organ <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 5)
Organ_ADARexp$Sex <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 7)
Organ_ADARexp$log2TPM <- log2(Organ_ADARexp$TPM)
head(Organ_ADARexp)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 5, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex","ADAR1exp","ADAR2exp"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
        editMeta[sample,"ADAR1exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADAR", "log2TPM"]
        editMeta[sample,"ADAR2exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADARB1", "log2TPM"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))
    ADAR1 <- editMeta$ADAR1exp
    ADAR2 <- editMeta$ADAR2exp


design <- model.matrix(~0+Stage+ADAR1+ADAR2)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_testis_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = F, row.names = T, col.names = T)


####### Kidney

# try checking site sparsity and then maybe imputation
edit_table_cb_kidney <- edit_table_cb[edit_table_cb$Organ == "Kidney",]
dim(edit_table_cb_kidney) # 22591482       16
length(unique(edit_table_cb_kidney$Chrom_Pos)) # 2140833

#edit_table_cb_kidney <- edit_table_cb_kidney[edit_table_cb_kidney$NonRefEdited >= 3,]
#dim(edit_table_cb_kidney) # 570500       16
edit_table_cb_kidney <- edit_table_cb_kidney[edit_table_cb_kidney$TotalEditedReads >= 5,]
dim(edit_table_cb_kidney) # 22591482       16

edit_table_cb_kidney <- edit_table_cb_kidney[edit_table_cb_kidney$Chrom != "M",]
edit_table_cb_kidney <- edit_table_cb_kidney[!edit_table_cb_kidney$Chrom_Pos %in% SNP_sites,]
edit_table_cb_kidney <- edit_table_cb_kidney[!grepl("^AC|^AL|^AP", edit_table_cb_kidney$Gene),]

dim(edit_table_cb_kidney) # 16477234       16
length(unique(edit_table_cb_kidney$Chrom_Pos)) # 1554842

sample_cutoff = 0.6
editlevel_cutoff = 0.05


sample_count <- length(table(edit_table_cb_kidney$Subject))
sample_count # 40
highquality_sites_coverage <- names(table(edit_table_cb_kidney$Chrom_Pos)[which(table(edit_table_cb_kidney$Chrom_Pos) > (sample_cutoff*sample_count))])
length(highquality_sites_coverage) # 208870

highquality_sites_minediting_table <- aggregate(edit_table_cb_kidney$EditLevel, list(edit_table_cb_kidney$Chrom_Pos), FUN=mean)
head(highquality_sites_minediting_table)
highquality_sites_minediting <- highquality_sites_minediting_table[highquality_sites_minediting_table$x >= editlevel_cutoff, "Group.1"]
length(highquality_sites_minediting) # 139813

length(intersect(highquality_sites_coverage, highquality_sites_minediting)) # 6850
highquality_sites_kidney <- intersect(highquality_sites_coverage, highquality_sites_minediting)

edit_table_cb_kidney_qc <- edit_table_cb_kidney[edit_table_cb_kidney$Chrom_Pos %in% highquality_sites_kidney, ]
dim(edit_table_cb_kidney_qc) # 235702
summary(edit_table_cb_kidney_qc)
write.table(edit_table_cb_kidney_qc,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_kidney_qc_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)


OrganDev_meta <- read.table('Figure_OrganDevelopment_R1/OrganDevelopment_Human_SuppTable_updatedv1.csv', sep = ',', header = T)

head(OrganDev_meta)
OrganDev_meta$Developmental_stage_sum <- NA
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc'), "Developmental_stage_sum"] <- "Early Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc'), "Developmental_stage_sum"] <- "Late Gestation"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('newborn','infant','toddler','school','youngTeenager','teenager'), "Developmental_stage_sum"] <- "Newborn - Teenager"
OrganDev_meta[OrganDev_meta$Developmental.stage %in% c('oldTeenager','youngAdult','youngMidAge','olderMidAge','senior'), "Developmental_stage_sum"] <- "Adult - Senior"

OrganDev_meta$Stage <- NA
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Early Gestation","Late Gestation"), "Stage"] <- "Prenatal"
OrganDev_meta[OrganDev_meta$Developmental_stage_sum %in% c("Newborn - Teenager","Adult - Senior"), "Stage"] <- "Postnatal"


sample_groups <- OrganDev_meta$library.ID
sites_of_interest <- highquality_sites_kidney
editExprs <- data.frame(matrix(NA, nrow = length(sites_of_interest), ncol = length(sample_groups), dimnames = list(sites_of_interest, sample_groups)))
    colnames(editExprs) <- sample_groups
    for (k in 1:nrow(edit_table_cb_kidney_qc)){
        editExprs_sites <- edit_table_cb_kidney_qc[k, "Chrom_Pos"]
        editExprs_sample <- edit_table_cb_kidney_qc[k, "Subject"]
        editExprs[editExprs_sites, as.character(editExprs_sample)] <- edit_table_cb_kidney_qc[k, "EditLevel"]
    }
dim(editExprs)
highquality_samples <- colnames(is.na(editExprs))[colSums(is.na(editExprs)) < 0.2*nrow(editExprs)]


table(OrganDev_meta[OrganDev_meta$library.ID %in% highquality_samples, "Developmental_stage_sum"])
# Early Gestation     Late Gestation Newborn - Teenager 
#                12                 13                  3 


editExprs_hq <- editExprs[colnames(editExprs) %in% highquality_samples]
dim(editExprs_hq)

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(editExprs_hq,2,pMiss)
max(apply(editExprs_hq,2,pMiss)) # 17.1%
mean(apply(editExprs_hq,2,pMiss)) # 4.93%
apply(editExprs_hq,1,pMiss)
max(apply(editExprs_hq,1,pMiss)) # 32.1%
mean(apply(editExprs_hq,1,pMiss)) # 4.93%

aggr_plot <- aggr(editExprs_hq, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))


# after dropping the low quality sample, clean up the sites one more time
highquality_sites_new <- names(apply(editExprs_hq,1,pMiss))[apply(editExprs_hq,1,pMiss) < (sample_cutoff*100)]
length(highquality_sites_new)
highquality_sites_new <- setdiff(highquality_sites_new, names(which(rowSums(editExprs_hq, na.rm = T) <= editlevel_cutoff)))
length(highquality_sites_new)

editExprs_hq_new <- editExprs_hq[row.names(editExprs_hq) %in% highquality_sites_new,]

apply(editExprs_hq_new,2,pMiss)
max(apply(editExprs_hq_new,2,pMiss))
mean(apply(editExprs_hq_new,2,pMiss))
apply(editExprs_hq_new,1,pMiss)
max(apply(editExprs_hq_new,1,pMiss))
mean(apply(editExprs_hq_new,1,pMiss))

# try the imputation

library(mice)
md.pattern(editExprs_hq_new)

library(VIM)
aggr_plot <- aggr(editExprs_hq_new, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

colnames(editExprs_hq_new) <- paste0("Library",colnames(editExprs_hq_new))
write.table(editExprs_hq_new,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_preimputedresults_kidney_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_impute <- mice(editExprs_hq_new,m=5,maxit=30,meth="pmm",seed=500)

editExprs_complete <- complete(editExprs_impute,1)
write.table(editExprs_complete,"Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_kidney_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)
editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_kidney_60per.txt", sep = ",", header = T, quote = "")

densityplot(editExprs_impute)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 3, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))

design <- model.matrix(~0+Stage+Sex)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_kidney_prenatalVSpostnatal_60per_noADAR.txt", sep = "," ,quote = F, row.names = T, col.names = T)


editExprs_complete <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_kidney_60per.txt", sep = ",", header = T, quote = "")

Organ_ADARexp <- read.table("git_repo/RNAEditing/Figure2/ADARexp_OrganDev.txt", header = T)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 3)
Organ_ADARexp$Library <- sapply(strsplit(Organ_ADARexp$Library, split = "sT"), "[[", 1)
Organ_ADARexp$Organ <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 5)
Organ_ADARexp$Sex <- sapply(strsplit(Organ_ADARexp$Sample, split = "_"), "[[", 7)
Organ_ADARexp$log2TPM <- log2(Organ_ADARexp$TPM)
head(Organ_ADARexp)

editMeta <- data.frame(matrix(NA, nrow = ncol(editExprs_complete), ncol = 5, dimnames = list(colnames(editExprs_complete), c("Subject","Stage","Sex","ADAR1exp","ADAR2exp"))))
row.names(editMeta) <- colnames(editExprs_complete)
    editMeta$Subject <- sapply(strsplit(colnames(editExprs_complete), split = "_"), "[[", 1)
    for (sample in row.names(editMeta)){
        sample_subject <- strsplit(sample, "Library")[[1]][2]
        sample_devstage <-  OrganDev_meta[OrganDev_meta$library.ID == sample_subject, "Stage"]
        editMeta[sample,"Stage"] <- sample_devstage
        editMeta[sample,"Sex"] <- OrganDev_meta[OrganDev_meta$library.ID == sample_subject & OrganDev_meta$Stage == sample_devstage, "Sex"]
        editMeta[sample,"ADAR1exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADAR", "log2TPM"]
        editMeta[sample,"ADAR2exp"] <- Organ_ADARexp[Organ_ADARexp$Library == sample_subject & Organ_ADARexp$Gene == "ADARB1", "log2TPM"]
    }

    Subject = as.factor(editMeta$Subject)
    Sex = as.factor(editMeta$Sex)
    Stage <- factor(editMeta$Stage, levels = c("Prenatal","Postnatal"))
    ADAR1 <- editMeta$ADAR1exp
    ADAR2 <- editMeta$ADAR2exp


design <- model.matrix(~0+Stage+Sex+ADAR1+ADAR2)

fit <- lmFit(editExprs_complete,design)
    cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fitDupCor <- eBayes(fit2)

    results_table <- topTable(fitDupCor, coef="DevEffect", n=nrow(editExprs_complete))
    dim(results_table)
    results_table <- results_table[!is.na(results_table$adj.P.Val), ]
    dim(results_table)

results_table_sig <- results_table[results_table$adj.P.Val < 0.05,]
dim(results_table_sig)
for (n in row.names(results_table_sig)){
        results_table_sig[n,"ALU"] <- REDIref[REDIref$Chrom_Pos == n, "type"]
        results_table_sig[n,"Gene"] <- REDIref[REDIref$Chrom_Pos == n, "Gene.wgEncodeGencodeBasicV34lift37"]
        results_table_sig[n,"Repeat"] <- REDIref[REDIref$Chrom_Pos == n, "repeat."]
        results_table_sig[n,"Region"] <- REDIref[REDIref$Chrom_Pos == n, "Func.wgEncodeGencodeBasicV34lift37"]
    }
#results_table_sig <- results_table_sig[which(sapply(strsplit(row.names(results_table_sig), split = "_"), "[[", 1) != "M"),]
write.table(results_table_sig, "Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_kidney_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = F, row.names = T, col.names = T)

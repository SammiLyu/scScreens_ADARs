# plot for each organ (*6) (heart??)
# plot for forebrian & hindbrain combined (done)
# plot for all organs (done)


results_table_sig <- do.call(rbind, list(read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_testis_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_kidney_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T)))
head(results_table_sig)

results_table_sig <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_testis_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T)


results_table_sig <- list(read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_liver_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_kidney_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_testis_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T)
)
names(results_table_sig) <- c("Forebrain","Hindbrain","Liver","Kidney","Testis")

for (i in names(results_table_sig)){
    results_table_sig[[i]]$Chrom_Pos <- row.names(results_table_sig[[i]])
    results_table_sig[[i]]$Organ <- i
}

#results_table_sig <- do.call(rbind, results_table_sig)

# head(results_table_sig[[1]])
# table(results_table_sig[[5]]$Region)
# sapply(sapply(results_table_sig, "[[", "Region"), table)
# unique(unlist(sapply(sapply(sapply(results_table_sig, "[[", "Region"), table), names)))

 site_count_table <- data.frame(region = rep(unique(unlist(sapply(sapply(sapply(results_table_sig, "[[", "Region"), table), names))), length(names(results_table_sig))),
             site_counts = NA, organ = rep(names(results_table_sig), each = length(unique(unlist(sapply(sapply(sapply(results_table_sig, "[[", "Region"), table), names))))))

 for (i in 1:nrow(site_count_table)){
     site_count_table[i, "site_counts"] <- sapply(sapply(results_table_sig, "[[", "Region"), table)[[site_count_table[i, "organ"]]][site_count_table[i, "region"]]
 }
 site_count_table
 site_count_table$region <- factor(site_count_table$region, levels = c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic"))

 pdf("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_forebrainhindbrainlivertestiskidney_prenatalVSpostnatal_60per_ADAR1ADAR2_sitesbyorgan.pdf")
 ggplot(data=site_count_table, aes(x=organ, y=site_counts, fill=region)) + scale_fill_brewer(palette="Set3", drop = F) + 
   geom_bar(stat="identity") + theme_minimal() + ggtitle("Genic Locations of Sites Captured in Each Organ") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) + 
     xlab("Organ") + ylab("Number of Sites")
 dev.off()

# pdf("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_forebrainhindbrainlivertestiskidney_prenatalVSpostnatal_60per_ADAR1ADAR2_sitesbyregion.pdf")
# ggplot(data=site_count_table, aes(x=region, y=site_counts, fill=organ)) + scale_fill_brewer(palette="Set3") + 
#   geom_bar(stat="identity") + theme_minimal() + ggtitle("Genic Locations of Sites Captured in Each Organ") +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
#     axis.title=element_text(size=14,face="bold")) + 
#     xlab("Region") + ylab("Number of Sites")
# dev.off()


#results_table_sig <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_heart_earlygestationVSlategestation_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T)

results_table_sig$logFC <- results_table_sig$logFC * 100

rpkm <- read.table("Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.txt")
head(rpkm)
colnames(rpkm)

# devtools::install_version("dbplyr", version = "2.3.4")
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- data.frame(genes = row.names(rpkm))
genes <- row.names(rpkm)
listAttributes(mart)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
df_m <- merge(df,G_list,by.x="genes",by.y="ensembl_gene_id")

results_table_sig$Gene2 <- NA
for (i in 1:nrow(results_table_sig)){
    if (length(grep(";",results_table_sig[i,"Gene"])) == 1)
    results_table_sig[i, "Gene2"] <- strsplit(results_table_sig[i, "Gene"], split = ";")[[1]][2]
    results_table_sig[i, "Gene"] <- strsplit(results_table_sig[i, "Gene"], split = ";")[[1]][1]
}

# damn this needs to be changed
rpkm_sub <- rpkm[,grep("Testis", colnames(rpkm))]
rpkm_sub <- rpkm
prenatal_group <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc', '11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc')
postnatal_group <- c('newborn','infant','toddler','school','youngTeenager','teenager', 'oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')

# some ensembl id's map to multiple genes
# for now just take average of all expression mapped to the gens with the same name
# underlying assumption is that genes located in different loci may have the same expression -> very likely not true, but hard to distinguish with current data
for (i in 1:nrow(results_table_sig)){
    if (!is.na(results_table_sig[i,"Gene"]) ){
    Gene_sym <- df_m[df_m$hgnc_symbol == results_table_sig[i,"Gene"],"genes"]
    results_table_sig[i, "GeneExprPrenatal"] <- mean(rowMeans(rpkm_sub[Gene_sym, sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% prenatal_group]))
    results_table_sig[i, "GeneExprPostnatal"] <- mean(rowMeans(rpkm_sub[Gene_sym, sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% postnatal_group]))
    results_table_sig[i, "GeneExprlfc"] <- log2(results_table_sig[i, "GeneExprPostnatal"] / results_table_sig[i, "GeneExprPrenatal"])
    }

    if (!is.na(results_table_sig[i,"Gene2"])){
    Gene2_sym <- df_m[df_m$hgnc_symbol == results_table_sig[i,"Gene2"],"genes"]
    results_table_sig[i, "Gene2ExprPrenatal"] <- mean(rowMeans(rpkm_sub[Gene2_sym, sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% prenatal_group]))
    results_table_sig[i, "Gene2ExprPostnatal"] <- mean(rowMeans(rpkm_sub[Gene2_sym, sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% postnatal_group]))
    results_table_sig[i, "Gene2Exprlfc"] <- log2(results_table_sig[i, "Gene2ExprPostnatal"] / results_table_sig[i, "Gene2ExprPrenatal"])
    }
}
head(results_table_sig,20)

#top10 <- row.names(results_table_sig[sort(results_table_sig$logFC, decreasing = T, index.return = T)$ix,])[c(1:10)]
#bot10 <- row.names(results_table_sig[sort(results_table_sig$logFC, decreasing = F, index.return = T)$ix,])[c(1:10)]

library(ggplot2)
data_to_plot <- data.frame(site = c(results_table_sig$Chrom_Pos, results_table_sig$Chrom_Pos),
    exprlfc = c(results_table_sig$GeneExprlfc, results_table_sig$Gene2Exprlfc),
                            editlfc = c(results_table_sig$logFC, results_table_sig$logFC),
                            gene=c(results_table_sig$Gene, results_table_sig$Gene2),
                            region=c(results_table_sig$Region, results_table_sig$Region))
data_to_plot <- data_to_plot[!is.na(data_to_plot$exprlfc),]

write.table(data_to_plot, file = "Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_forebrainhindbrainliverkidneytestis_prenatalVSpostnatal_60per_ADAR1ADAR2_summarytablebeforeplotting.txt", sep = ",", row.names = F, quote = F)

library(dplyr)
data_to_plot <- data_to_plot %>%
  group_by(gene) %>%
  summarise(meanexprlfc = mean(exprlfc), meaneditlfc = mean(editlfc), region = region, n = n()) %>%
  distinct()




#mycolors <- c("postnatally biased"="blue", "NO"="grey", "prenatally biased"="red")
#data_to_plot$rel <- "NO"
#data_to_plot$rel[data_to_plot$site %in% top10] <- "postnatally biased"
#data_to_plot$rel[data_to_plot$site %in% bot10] <- "prenatally biased"

data_to_plot_noninf <- data_to_plot[is.finite(rowSums(data_to_plot[c(2,3)])),]

cor_results <- cor.test(data_to_plot_noninf$meanexprlfc, data_to_plot_noninf$meaneditlfc, method = "pearson")

p1 <- ggplot(data_to_plot, aes(x=meanexprlfc, y=meaneditlfc)) + geom_point() + 
#scale_color_manual(values=mycolors) +
      theme_minimal() + xlim(-5,5) + ylim(-20,60) + geom_smooth(method = 'lm', fullrange = T) +
      geom_hline(yintercept = 0) + geom_vline (xintercept = 0)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_results$estimate, digits = 3), "\n ", "p = ", format(cor_results$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust
p1 <- p1 + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
ggtitle("Correlating Editing Rate Change to Expression Change Per Gene") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("Postnatal/Prenatal log2FC") + ylab("Average Delta Editing Rate/Gene (%)")
p1

cor_results_table <- data.frame(region =  c("All Regions Combined", names(table(data_to_plot$region))),
                                coeff = NA,
                                lower_conf = NA,
                                upper_conf = NA,
                                p_value = NA)
for (i in cor_results_table$region){
    # filter table first to remove regions with fewer than 4 genes
    if (i != "All Regions Combined" & table(data_to_plot$region)[i] < 4){
        cor_results_table <- cor_results_table[-which(cor_results_table$region == i),]
    }
}

for (i in cor_results_table$region){
    if (i == "All Regions Combined"){
        data_to_plot_sub <- data_to_plot[is.finite(rowSums(data_to_plot[c(2,3)])),]
    }
    else {
        data_to_plot_sub <- data_to_plot[data_to_plot$region == i,]
        data_to_plot_sub <- data_to_plot_sub[is.finite(rowSums(data_to_plot_sub[c(2,3)])),]
    }
    cor_results <- cor.test(data_to_plot_sub$meanexprlfc, data_to_plot_sub$meaneditlfc, method = "pearson")
    cor_results_table[cor_results_table$region == i, "coeff"] <- cor_results$estimate
    cor_results_table[cor_results_table$region == i, "lower_conf"] = cor_results$conf.int[1]
    cor_results_table[cor_results_table$region == i, "upper_conf"] = cor_results$conf.int[2]
    cor_results_table[cor_results_table$region == i, "p_value"] = cor_results$p.value
}

cor_results_table <- cor_results_table[order(cor_results_table$coeff, decreasing = F),]
cor_results_table$region <- factor(cor_results_table$region, levels = cor_results_table$region)

write.table(cor_results_table, file = "Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_forebrainhindbrainliverkidneytestis_prenatalVSpostnatal_60per_ADAR1ADAR2_figs8style.txt", sep = ",", row.names = F, quote = F)

p2 <- ggplot(cor_results_table) +
    geom_bar(aes(x=region, y=coeff, fill = region), stat="identity", alpha=0.7) +
    geom_errorbar( aes(x=region, ymin=lower_conf, ymax=upper_conf), width=0.2, colour="black", alpha=1, size=0.5) +
    theme_minimal() + ggtitle("Pearson Correlation Coefficient for Average Delta Editing vs. Gene Exp") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("Region") + ylab("Pearson Correlation Coefficient")
library(cowplot)

pdf("Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_forebrainhindbrainliverkidneytestis_prenatalVSpostnatal_60per_ADAR1ADAR2_figs8style.pdf", width = 15, height = 10)
plot_grid(p1,p2)
dev.off()


# drop regions with fewer than 20 genes
table(data_to_plot$region)

for (i in cor_results_table$region){
    # filter table first to remove regions with fewer than 20 genes
    if (i != "All Regions Combined" & table(data_to_plot$region)[i] < 20){
        cor_results_table <- cor_results_table[-which(cor_results_table$region == i),]
    }
}
write.table(cor_results_table, file = "Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_forebrainhindbrainliverkidneytestis_prenatalVSpostnatal_60per_ADAR1ADAR2_figs8style_min20genes.txt", sep = ",", row.names = F, quote = F)

p2 <- ggplot(cor_results_table) +
    geom_bar(aes(x=region, y=coeff, fill = region), stat="identity", alpha=0.7) +
    geom_errorbar( aes(x=region, ymin=lower_conf, ymax=upper_conf), width=0.2, colour="black", alpha=1, size=0.5) +
    theme_minimal() + ggtitle("Pearson Correlation Coefficient for Average Delta Editing vs. Gene Exp") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("Region") + ylab("Pearson Correlation Coefficient")
library(cowplot)

pdf("Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_forebrainhindbrainliverkidneytestis_prenatalVSpostnatal_60per_ADAR1ADAR2_figs8style_min20genes.pdf", width = 15, height = 10)
plot_grid(p1,p2)
dev.off()

for (i in cor_results_table$region){
    # filter table first to remove regions with fewer than 30 genes
    if (i != "All Regions Combined" & table(data_to_plot$region)[i] < 30){
        cor_results_table <- cor_results_table[-which(cor_results_table$region == i),]
    }
}
write.table(cor_results_table, file = "Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_forebrainhindbrainliverkidneytestis_prenatalVSpostnatal_60per_ADAR1ADAR2_figs8style_min30genes.txt", sep = ",", row.names = F, quote = F)

p2 <- ggplot(cor_results_table) +
    geom_bar(aes(x=region, y=coeff, fill = region), stat="identity", alpha=0.7) +
    geom_errorbar( aes(x=region, ymin=lower_conf, ymax=upper_conf), width=0.2, colour="black", alpha=1, size=0.5) +
    theme_minimal() + ggtitle("Pearson Correlation Coefficient for Average Delta Editing vs. Gene Exp") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("Region") + ylab("Pearson Correlation Coefficient")
library(cowplot)

pdf("Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_forebrainhindbrainliverkidneytestis_prenatalVSpostnatal_60per_ADAR1ADAR2_figs8style_min30genes.pdf", width = 15, height = 10)
plot_grid(p1,p2)
dev.off()


# for liver, additional scatter plot for sites located in intronic regions only
results_table_sig <- read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_liver_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T)

results_table_sig$logFC <- results_table_sig$logFC * 100

rpkm <- read.table("Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.txt")
head(rpkm)
colnames(rpkm)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- data.frame(genes = row.names(rpkm))
genes <- row.names(rpkm)
listAttributes(mart)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
df_m <- merge(df,G_list,by.x="genes",by.y="ensembl_gene_id")

results_table_sig <- results_table_sig[results_table_sig$Region == "intronic",]
results_table_sig$Chrom_Pos <- row.names(results_table_sig)

results_table_sig$Gene2 <- NA
for (i in 1:nrow(results_table_sig)){
    if (length(grep(";",results_table_sig[i,"Gene"])) == 1)
    results_table_sig[i, "Gene2"] <- strsplit(results_table_sig[i, "Gene"], split = ";")[[1]][2]
    results_table_sig[i, "Gene"] <- strsplit(results_table_sig[i, "Gene"], split = ";")[[1]][1]
}

rpkm_sub <- rpkm[,grep("Liver", colnames(rpkm))]
rpkm_sub <- rpkm
prenatal_group <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc', '11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc')
postnatal_group <- c('newborn','infant','toddler','school','youngTeenager','teenager', 'oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')

# some ensembl id's map to multiple genes
# for now just take average of all expression mapped to the gens with the same name
# underlying assumption is that genes located in different loci may have the same expression -> very likely not true, but hard to distinguish with current data
for (i in 1:nrow(results_table_sig)){
    if (!is.na(results_table_sig[i,"Gene"]) ){
    Gene_sym <- df_m[df_m$hgnc_symbol == results_table_sig[i,"Gene"],"genes"]
    results_table_sig[i, "GeneExprPrenatal"] <- mean(rowMeans(rpkm_sub[Gene_sym, sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% prenatal_group]))
    results_table_sig[i, "GeneExprPostnatal"] <- mean(rowMeans(rpkm_sub[Gene_sym, sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% postnatal_group]))
    results_table_sig[i, "GeneExprlfc"] <- log2(results_table_sig[i, "GeneExprPostnatal"] / results_table_sig[i, "GeneExprPrenatal"])
    }

    if (!is.na(results_table_sig[i,"Gene2"])){
    Gene2_sym <- df_m[df_m$hgnc_symbol == results_table_sig[i,"Gene2"],"genes"]
    results_table_sig[i, "Gene2ExprPrenatal"] <- mean(rowMeans(rpkm_sub[Gene2_sym, sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% prenatal_group]))
    results_table_sig[i, "Gene2ExprPostnatal"] <- mean(rowMeans(rpkm_sub[Gene2_sym, sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% postnatal_group]))
    results_table_sig[i, "Gene2Exprlfc"] <- log2(results_table_sig[i, "Gene2ExprPostnatal"] / results_table_sig[i, "Gene2ExprPrenatal"])
    }
}
head(results_table_sig,20)

library(ggplot2)
data_to_plot <- data.frame(site = c(results_table_sig$Chrom_Pos, results_table_sig$Chrom_Pos),
    exprlfc = c(results_table_sig$GeneExprlfc, results_table_sig$Gene2Exprlfc),
                            editlfc = c(results_table_sig$logFC, results_table_sig$logFC),
                            gene=c(results_table_sig$Gene, results_table_sig$Gene2),
                            region=c(results_table_sig$Region, results_table_sig$Region))
data_to_plot <- data_to_plot[!is.na(data_to_plot$exprlfc),]

write.table(data_to_plot, file = "Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_kidney_prenatalVSpostnatal_60per_ADAR1ADAR2_summarytablebeforeplotting_intronicregiononly.txt", sep = ",", row.names = F, quote = F)

library(dplyr)
data_to_plot <- data_to_plot %>%
  group_by(gene) %>%
  summarise(meanexprlfc = mean(exprlfc), meaneditlfc = mean(editlfc), region = region, n = n()) %>%
  distinct()

data_to_plot_noninf <- data_to_plot[is.finite(rowSums(data_to_plot[c(2,3)])),]

cor_results <- cor.test(data_to_plot_noninf$meanexprlfc, data_to_plot_noninf$meaneditlfc, method = "pearson")

p1 <- ggplot(data_to_plot, aes(x=meanexprlfc, y=meaneditlfc)) + geom_point() + 
#scale_color_manual(values=mycolors) +
      theme_minimal() + xlim(-5,5) + ylim(-20,60) + geom_smooth(method = 'lm', fullrange = T) +
      geom_hline(yintercept = 0) + geom_vline (xintercept = 0)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_results$estimate, digits = 3), "\n ", "p = ", format(cor_results$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust
p1 <- p1 + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
ggtitle("Correlating Editing Rate Change to Expression Change Per Gene") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("Postnatal/Prenatal log2FC") + ylab("Average Delta Editing Rate/Gene (%)")
p1

pdf("Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_liver_prenatalVSpostnatal_60per_ADAR1ADAR2_figs8style_intronicregionly.pdf", width = 7.5, height = 10)
p1
dev.off()






transparent_theme <- theme(
 axis.title.x = element_blank(),
 axis.title.y = element_blank(),
 axis.text.x = element_blank(), 
 axis.text.y = element_blank(),
 axis.ticks = element_blank(),
 panel.grid = element_blank(),
 axis.line = element_blank(),
 panel.background = element_rect(fill = "transparent",colour = NA),
 plot.background = element_rect(fill = "transparent",colour = NA))

p2 <- ggplot(data_to_plot, aes(x = exprlfc)) + geom_boxplot() + transparent_theme + xlim(-7,7)
p3 <- ggplot(data_to_plot, aes(x = editlfc)) + geom_boxplot() + coord_flip() + transparent_theme + xlim(-0.2,0.6)

pdf("Figure_OrganDevelopment_R1/Corr_Plots/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_70per_ADAR1ADAR2.pdf", width = 10, height = 10)
cowplot:::plot_grid(p3,p1,NA,p2, nrow= 2, ncol = 2, rel_widths = c(1,15), rel_heights = c(20,1))
dev.off()

setdiff(results_table_sig$Gene, df_m$hgnc_symbol)
setdiff(results_table_sig$Gene2, df_m$hgnc_symbol)

setdiff(df_m$hgnc_symbol, results_table_sig$Gene)




## plot editing levels and expression levels across entire development period for a specific gene (consensus genes?)
## version updated by Feb28 / Rerun with other genes Mar26

edit_level_imputed <- list(read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_forebrain_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_hindbrain_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_liver_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_kidney_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_testis_60per.txt", sep = "," ,quote = "", row.names = 1, header = T)
)
names(edit_level_imputed) <- c("Forebrain","Hindbrain","Liver","Kidney","Testis")

REDIref <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_hg19_mod_nochr.txt", header = T, sep = "\t")
REDIref$Chrom_Pos <- paste0(REDIref$Region,"_",REDIref$Position)


rpkm <- read.table("Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.txt")
colnames(rpkm)[which(colnames(rpkm) == "Testis.Senior.41")] <- "Testis.senior.41"

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- data.frame(genes = row.names(rpkm))
genes <- row.names(rpkm)
#listAttributes(mart)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
df_m <- merge(df,G_list,by.x="genes",by.y="ensembl_gene_id")

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

lat.order.levels <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc','11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc','newborn','infant','toddler','school','youngTeenager','teenager','oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')
DevelopmentStage.data <- data.frame(
  DevelopmentStage = lat.order.levels,
  DevelopmentStage_num = 1:length(lat.order.levels)
) # metadata obtained

# genes checked: 
# VPS41, NUP43, ZDHHC20, PRKCSH, PSMD12, MAVS, PPIA, SNRPD3, EIF2AK2, PSMB2
gene_to_check <- c("VPS41", "NUP43", "ZDHHC20", "PRKCSH", "PSMD12", "MAVS", "PPIA", "SNRPD3", "EIF2AK2", "PSMB2")
gene_to_check <- c("EIF2AK2", "MAVS", "GINS1", "H2AZ2")


for (gene in gene_to_check){

#gene <- "VPS41"
#edit_level_imputed <- list(read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_forebrain_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
#read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_hindbrain_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
#read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_liver_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
#read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_kidney_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
#read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_testis_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T)
#)
#names(edit_level_imputed) <- c("Forebrain","Hindbrain","Liver","Kidney","Testis")

#edit_level_imputed_sub <- list()
#for (i in names(edit_level_imputed)){
#    edit_level_imputed_sub[[i]] <- edit_level_imputed[[i]][edit_level_imputed[[i]]$Gene == gene,]
#}

sites_of_interest <- REDIref[grep(gene, REDIref$Gene.wgEncodeGencodeBasicV34lift37),"Chrom_Pos"]

edit_level_imputed_sub <- list()
for (i in names(edit_level_imputed)){
    edit_level_imputed_sub[[i]] <- edit_level_imputed[[i]][row.names(edit_level_imputed[[i]]) %in% sites_of_interest,]
}

edit_level_imputed_sub_colMean <- list()
for (i in names(edit_level_imputed_sub)){
    edit_level_imputed_sub_colMean[[i]] <- colMeans(edit_level_imputed_sub[[i]]) # mean edit level of gene obtained
}


Gene_sym <- df_m[df_m$hgnc_symbol == gene,"genes"]

rpkm_gene <- rpkm[Gene_sym,]
row.names(rpkm_gene) <- "GeneExpr"
rpkm_gene["DevelopmentStage", ] <- sapply(strsplit(colnames(rpkm_gene), split = "[.]"), "[[", 2) # expression level of gene obtained


library(ggplot2)
library(cowplot)

p <- list()
for (i in names(edit_level_imputed)){
    edit_to_plot <- data.frame(edit_level_imputed_sub_colMean[[i]])
    colnames(edit_to_plot) <- "EditLevel"
    for (j in row.names(edit_to_plot)){
        edit_to_plot[j, "DevelopmentStage"] <- OrganDev_meta[paste0("Library",OrganDev_meta$library.ID) == j, "Developmental.stage"]
    }
    if (i == "Forebrain"){
        expr_to_plot <- rpkm_gene[,grep("Brain", colnames(rpkm))]
    }
    else if (i == "Hindbrain"){
        expr_to_plot <- rpkm_gene[,grep("Cerebellum", colnames(rpkm))]
    }
    else {
        expr_to_plot <- rpkm_gene[,grep(i, colnames(rpkm))]
    }
    edit_to_plot <- dplyr::left_join(edit_to_plot, DevelopmentStage.data, by= "DevelopmentStage")
    expr_to_plot <- dplyr::left_join(data.frame(t(expr_to_plot)), DevelopmentStage.data, by= "DevelopmentStage")
    expr_to_plot$GeneExpr <- log2(as.numeric(expr_to_plot$GeneExpr))

    p[[paste0(i, "_Edit")]] <- ggplot(edit_to_plot, aes(x = DevelopmentStage_num, y = EditLevel)) + geom_point(size=5, color = "darkblue") + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + 
  theme_classic() + ggtitle(paste0("Mean Editing Level in ",i)) + coord_cartesian(xlim = c(1,25)) + 
theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)

    p[[paste0(i, "_Expr")]] <- ggplot(expr_to_plot, aes(x = DevelopmentStage_num, y = GeneExpr)) + geom_point(size=5, color = "darkblue") + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + 
  theme_classic() + ggtitle(paste0("Gene Expression Level in ",i)) + coord_cartesian(xlim = c(1,25)) + 
theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)

edit_expr_to_plot <- dplyr::inner_join(edit_to_plot, expr_to_plot, by= "DevelopmentStage")

    cor_results <- cor.test(edit_expr_to_plot$GeneExpr, edit_expr_to_plot$EditLevel, method = "pearson")
    annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_results$estimate, digits = 3), "\n ", "p = ", format(cor_results$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust 

    p[[paste0(i, "_Corr")]] <- ggplot(edit_expr_to_plot, aes(x=GeneExpr, y=EditLevel)) + geom_point(size = 5, color = "darkblue") + 
      theme_classic() + geom_smooth(method = 'lm', fullrange = T, se = FALSE) + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
ggtitle("Correlating Editing Rate Change to Expression Change") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 18),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("log2RPKM") + ylab("Average Editing Rate")
}

pdf(paste0("Figure_OrganDevelopment_R1/GeneSpecificPlots/",gene,"_editingplusgeneexpr_20240326.pdf"), width = 30, height = 40)
print(plot_grid(plotlist = p, nrow = 5, ncol = 3))
dev.off()

}




## plot editing levels and expression levels across entire development period for a specific gene (consensus genes?)

# genes checked: 
# PSMD12 - "17_65334753","17_65334856","17_65334867"
# EIF2AK2 - "2_37328075","37328082"
edit_level_imputed <- list(read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_forebrain_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_hindbrain_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_liver_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_kidney_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_testis_60per.txt", sep = "," ,quote = "", row.names = 1, header = T)
)
names(edit_level_imputed) <- c("Forebrain","Hindbrain","Liver","Kidney","Testis")

edit_level_imputed_sub <- list()
for (i in names(edit_level_imputed)){
    edit_level_imputed_sub[[i]] <- edit_level_imputed[[i]][row.names(edit_level_imputed[[i]]) %in% c("2_37328075","37328082"),]
}

edit_level_imputed_sub_colMean <- list()
for (i in names(edit_level_imputed_sub)){
    edit_level_imputed_sub_colMean[[i]] <- colMeans(edit_level_imputed_sub[[i]]) # mean edit level of gene obtained
}


rpkm <- read.table("Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.txt")

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- data.frame(genes = row.names(rpkm))
genes <- row.names(rpkm)
listAttributes(mart)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
df_m <- merge(df,G_list,by.x="genes",by.y="ensembl_gene_id")

Gene_sym <- df_m[df_m$hgnc_symbol == "EIF2AK2","genes"]

rpkm_gene <- rpkm[Gene_sym,]
row.names(rpkm_gene) <- "GeneExpr"
rpkm_gene["DevelopmentStage", ] <- sapply(strsplit(colnames(rpkm_gene), split = "[.]"), "[[", 2) # expression level of gene obtained


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

lat.order.levels <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc','11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc','newborn','infant','toddler','school','youngTeenager','teenager','oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')
DevelopmentStage.data <- data.frame(
  DevelopmentStage = lat.order.levels,
  DevelopmentStage_num = 1:length(lat.order.levels)
) # metadata obtained


p <- list()
for (i in names(edit_level_imputed)){
    edit_to_plot <- data.frame(edit_level_imputed_sub_colMean[[i]])
    colnames(edit_to_plot) <- "EditLevel"
    for (j in row.names(edit_to_plot)){
        edit_to_plot[j, "DevelopmentStage"] <- OrganDev_meta[paste0("Library",OrganDev_meta$library.ID) == j, "Developmental.stage"]
    }
    if (i == "Forebrain"){
        expr_to_plot <- rpkm_gene[,grep("Brain", colnames(rpkm))]
    } else if (i == "Hindbrain"){
        expr_to_plot <- rpkm_gene[,grep("Cerebellum", colnames(rpkm))]
    } else {
        expr_to_plot <- rpkm_gene[,grep(i, colnames(rpkm))]
    }
    edit_to_plot <- dplyr::left_join(edit_to_plot, DevelopmentStage.data, by= "DevelopmentStage")
    expr_to_plot <- dplyr::left_join(data.frame(t(expr_to_plot)), DevelopmentStage.data, by= "DevelopmentStage")
    expr_to_plot$GeneExpr <- log2(as.numeric(expr_to_plot$GeneExpr))

    p[[paste0(i, "_Edit")]] <- ggplot(edit_to_plot, aes(x = DevelopmentStage_num, y = EditLevel)) + geom_point(size=5, color = "darkblue") + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + 
  theme_classic() + ggtitle(paste0("Mean Editing Level in ",i)) + coord_cartesian(xlim = c(1,25)) + 
theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)

    p[[paste0(i, "_Expr")]] <- ggplot(expr_to_plot, aes(x = DevelopmentStage_num, y = GeneExpr)) + geom_point(size=5, color = "darkblue") + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + 
  theme_classic() + ggtitle(paste0("Gene Expression Level in ",i)) + coord_cartesian(xlim = c(1,25)) + 
theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)

edit_expr_to_plot <- dplyr::inner_join(edit_to_plot, expr_to_plot, by= "DevelopmentStage")

    cor_results <- cor.test(edit_expr_to_plot$GeneExpr, edit_expr_to_plot$EditLevel, method = "pearson")
    annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_results$estimate, digits = 3), "\n ", "p = ", format(cor_results$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust 

    p[[paste0(i, "_Corr")]] <- ggplot(edit_expr_to_plot, aes(x=GeneExpr, y=EditLevel)) + geom_point(size = 5, color = "darkblue") + 
      theme_classic() + geom_smooth(method = 'lm', fullrange = T, se = FALSE) + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
ggtitle("Correlating Editing Rate Change to Expression Change") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 18),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("log2RPKM") + ylab("Average Editing Rate")
}

pdf("Figure_OrganDevelopment_R1/GeneSpecificPlots/EIF2AK2_editingplusgeneexpr.pdf", width = 30, height = 40)
plot_grid(plotlist = p, nrow = 5, ncol = 3)
dev.off()









    metadata_df_sub <- reshape2:::melt(results_table_sig_sub[results_table_sig_sub$Organ == i,c("ADAR1exp","ADAR2exp","AEI_A2GEditing","DevelopmentStage_num")], id.vars="DevelopmentStage_num")
    metadata_df_sub[metadata_df_sub$variable == "AEI_A2GEditing","value"] <- metadata_df_sub[metadata_df_sub$variable == "AEI_A2GEditing","value"] * coeff + min_exp

ggplot(metadata_df_sub, aes(x = DevelopmentStage_num, y = value)) + geom_point(size=7.5, aes(color = variable, shape = variable)) + scale_color_manual(values = c('#619CFF','#00BA38','#F8766D')) + scale_shape_manual(values = c(21,21,16)) +
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



### March 21 request - plotting logDiffEdit vs. logFC for convergence sites
results_table_sig <- list(read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_forebrain_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_hindbrain_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_liver_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_kidney_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/SigDiffEditSites_scriptv5_testis_prenatalVSpostnatal_60per_ADAR1ADAR2.txt", sep = "," ,quote = "", row.names = 1, header = T)
)
names(results_table_sig) <- c("Forebrain","Hindbrain","Liver","Kidney","Testis")

for (i in names(results_table_sig)){
    results_table_sig[[i]]$Chrom_Pos <- row.names(results_table_sig[[i]])
    results_table_sig[[i]]$Organ <- i
}

convergence_sites <- read.table("Figure_OrganDevelopment_R1/Corr_Plots/convergence_sites.csv", skip = 1, sep = ",", quote = "")
colnames(convergence_sites) <- c("Site","Gene","Region")
head(convergence_sites)

for (i in names(results_table_sig)){
    results_table_sig[[i]] <- results_table_sig[[i]][results_table_sig[[i]]$Chrom_Pos %in% convergence_sites$Site,]
}

sapply(results_table_sig, dim)

results_table_toplot <- do.call(rbind, results_table_sig)
results_table_toplot$logFC <- results_table_toplot$logFC * 100

rpkm <- read.table("Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.txt")
head(rpkm)
colnames(rpkm)

# devtools::install_version("dbplyr", version = "2.3.4")
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- data.frame(genes = row.names(rpkm))
genes <- row.names(rpkm)
listAttributes(mart)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
df_m <- merge(df,G_list,by.x="genes",by.y="ensembl_gene_id")

results_table_toplot$Gene2 <- NA
for (i in 1:nrow(results_table_toplot)){
    if (length(grep(";",results_table_toplot[i,"Gene"])) == 1)
    results_table_toplot[i, "Gene2"] <- strsplit(results_table_toplot[i, "Gene"], split = ";")[[1]][2]
    results_table_toplot[i, "Gene"] <- strsplit(results_table_toplot[i, "Gene"], split = ";")[[1]][1]
}

rpkm_sub <- rpkm
prenatal_group <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc', '11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc')
postnatal_group <- c('newborn','infant','toddler','school','youngTeenager','teenager', 'oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')

for (i in 1:nrow(results_table_toplot)){
    Organ_tmp <- results_table_toplot[i, "Organ"]
    if (Organ_tmp == "Forebrain"){Organ_tmp <- "Brain"} else if (Organ_tmp == "Hindbrain"){Organ_tmp <- "Cerebellum"}
    if (!is.na(results_table_toplot[i,"Gene"]) ){
    Gene_sym <- df_m[df_m$hgnc_symbol == results_table_toplot[i,"Gene"],"genes"]
    results_table_toplot[i, "GeneExprPrenatal"] <- mean(rowMeans(rpkm_sub[Gene_sym, (sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% prenatal_group) & (sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",1) == Organ_tmp)]))
    results_table_toplot[i, "GeneExprPostnatal"] <- mean(rowMeans(rpkm_sub[Gene_sym, (sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% postnatal_group) & (sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",1) == Organ_tmp)]))
    results_table_toplot[i, "GeneExprlfc"] <- log2(results_table_toplot[i, "GeneExprPostnatal"] / results_table_toplot[i, "GeneExprPrenatal"])
    }

    if (!is.na(results_table_toplot[i,"Gene2"])){
    Gene2_sym <- df_m[df_m$hgnc_symbol == results_table_toplot[i,"Gene2"],"genes"]
    results_table_toplot[i, "Gene2ExprPrenatal"] <- mean(rowMeans(rpkm_sub[Gene2_sym, (sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% prenatal_group) & (sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",1) == Organ_tmp)]))
    results_table_toplot[i, "Gene2ExprPostnatal"] <- mean(rowMeans(rpkm_sub[Gene2_sym, (sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",2) %in% postnatal_group) & (sapply(strsplit(colnames(rpkm_sub), split = "[.]"),"[[",1) == Organ_tmp)]))
    results_table_toplot[i, "Gene2Exprlfc"] <- log2(results_table_toplot[i, "Gene2ExprPostnatal"] / results_table_toplot[i, "Gene2ExprPrenatal"])
    }
}
head(results_table_toplot,20)
results_table_toplot[results_table_toplot$Chrom_Pos == "19_11561185",]

results_table_toplot$Organ_Chrom_Pos <- paste0(results_table_toplot$Organ, "_", results_table_toplot$Chrom_Pos)
data_to_plot <- data.frame(site = c(results_table_toplot$Organ_Chrom_Pos, results_table_toplot$Organ_Chrom_Pos),
    exprlfc = c(results_table_toplot$GeneExprlfc, results_table_toplot$Gene2Exprlfc),
                            editlfc = c(results_table_toplot$logFC, results_table_toplot$logFC),
                            gene=c(results_table_toplot$Gene, results_table_toplot$Gene2),
                            region=c(results_table_toplot$Region, results_table_toplot$Region),
                            organ=c(results_table_toplot$Organ, results_table_toplot$Organ))
data_to_plot <- data_to_plot[!is.na(data_to_plot$exprlfc),]
data_to_plot[data_to_plot$gene == "MCM4",]

min(data_to_plot$exprlfc) ; max(data_to_plot$exprlfc)
min(data_to_plot$editlfc) ; max(data_to_plot$editlfc)

cor_results <- cor.test(data_to_plot$exprlfc, data_to_plot$editlfc, method = "pearson")

library(ggplot2)
p1 <- ggplot(data_to_plot, aes(x=exprlfc, y=editlfc)) + geom_point() + 
#scale_color_manual(values=mycolors) +
      theme_minimal() + xlim(-5,5) + ylim(-20,80) + geom_smooth(method = 'lm', fullrange = T) +
      geom_hline(yintercept = 0) + geom_vline (xintercept = 0)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_results$estimate, digits = 3), "\n ", "p = ", format(cor_results$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust
p1 <- p1 + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
ggtitle("Correlating Editing Rate Change to Expression Change Per Site Per Organ for Convergent Genes") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("Postnatal/Prenatal log2FC") + ylab("Average Delta Editing Rate/Gene (%)")
p1

pdf("Figure_OrganDevelopment_R1/Corr_Plots/CorrelatingEditingRatetoExpression_PerSitePerOrgan_ConvergentGenes.pdf", width = 10, height = 8)
p1
dev.off()

library(dplyr)
data_to_plot$gene_organ <- paste0(data_to_plot$gene, "_", data_to_plot$organ)
data_to_plot <- data_to_plot %>%
  group_by(gene_organ) %>%
  summarise(meanexprlfc = mean(exprlfc), meaneditlfc = mean(editlfc), region = region, n = n()) %>%
  distinct()

min(data_to_plot$meanexprlfc) ; max(data_to_plot$meanexprlfc)
min(data_to_plot$meaneditlfc) ; max(data_to_plot$meaneditlfc)

cor_results <- cor.test(data_to_plot$meanexprlfc, data_to_plot$meaneditlfc, method = "pearson")

p1 <- ggplot(data_to_plot, aes(x=meanexprlfc, y=meaneditlfc)) + geom_point() + 
#scale_color_manual(values=mycolors) +
      theme_minimal() + xlim(-5,5) + ylim(-20,80) + geom_smooth(method = 'lm', fullrange = T) +
      geom_hline(yintercept = 0) + geom_vline (xintercept = 0)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_results$estimate, digits = 3), "\n ", "p = ", format(cor_results$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust
p1 <- p1 + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
ggtitle("Correlating Editing Rate Change to Expression Change Per Gene Per Organ for Convergent Genes") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("Postnatal/Prenatal log2FC") + ylab("Average Delta Editing Rate/Gene (%)")
p1

pdf("Figure_OrganDevelopment_R1/Corr_Plots/CorrelatingEditingRatetoExpression_PerGenePerOrgan_ConvergentGenes.pdf", width = 10, height = 8)
p1
dev.off()


### April 4 request - plotting logDiffEdit vs. logFC for intronic sites in selected genes for Liver only

edit_level_imputed <- list(read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_forebrain_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_hindbrain_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_liver_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_kidney_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_imputedresults_testis_60per.txt", sep = "," ,quote = "", row.names = 1, header = T)
)
names(edit_level_imputed) <- c("Forebrain","Hindbrain","Liver","Kidney","Testis")

REDIref <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_hg19_mod_nochr.txt", header = T, sep = "\t")
REDIref$Chrom_Pos <- paste0(REDIref$Region,"_",REDIref$Position)


rpkm <- read.table("Figure_OrganDevelopment_R1/OrganDev.Human.RPKM.txt")
colnames(rpkm)[which(colnames(rpkm) == "Testis.Senior.41")] <- "Testis.senior.41"

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- data.frame(genes = row.names(rpkm))
genes <- row.names(rpkm)
#listAttributes(mart)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
df_m <- merge(df,G_list,by.x="genes",by.y="ensembl_gene_id")

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

lat.order.levels <- c('4wpc','5wpc','6wpc','7wpc','8wpc','9wpc','10wpc','11wpc','12wpc','13wpc','16wpc','18wpc','19wpc','20wpc','newborn','infant','toddler','school','youngTeenager','teenager','oldTeenager','youngAdult','youngMidAge','olderMidAge','senior')
DevelopmentStage.data <- data.frame(
  DevelopmentStage = lat.order.levels,
  DevelopmentStage_num = 1:length(lat.order.levels)
) # metadata obtained

gene_to_check <- c("SREBF1","HMGCS1","FDPS","FADS2","EBP")
ensg_to_check <- c("ENSG00000072310","ENSG00000112972","ENSG00000160752","ENSG00000134824","ENSG00000147155")
names(gene_to_check) <- ensg_to_check

REDIref_intronic <- REDIref[REDIref$Func.wgEncodeGencodeBasicV34lift37 == "intronic",]

for (gene in gene_to_check){

#gene <- "VPS41"
#edit_level_imputed <- list(read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_forebrain_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
#read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_hindbrain_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
#read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_liver_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
#read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_kidney_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T),
#read.table("Figure_OrganDevelopment_R1/20240214_allsigtables/hq_sites_sample_edit_table_cb_testis_qc_60per.txt", sep = "," ,quote = "", row.names = 1, header = T)
#)
#names(edit_level_imputed) <- c("Forebrain","Hindbrain","Liver","Kidney","Testis")

#edit_level_imputed_sub <- list()
#for (i in names(edit_level_imputed)){
#    edit_level_imputed_sub[[i]] <- edit_level_imputed[[i]][edit_level_imputed[[i]]$Gene == gene,]
#}

sites_of_interest <- REDIref_intronic[grep(gene, REDIref_intronic$Gene.wgEncodeGencodeBasicV34lift37),"Chrom_Pos"]

edit_level_imputed_sub <- list()
for (i in names(edit_level_imputed)){
    edit_level_imputed_sub[[i]] <- edit_level_imputed[[i]][row.names(edit_level_imputed[[i]]) %in% sites_of_interest,]
}

edit_level_imputed_sub_colMean <- list()
for (i in names(edit_level_imputed_sub)){
    edit_level_imputed_sub_colMean[[i]] <- colMeans(edit_level_imputed_sub[[i]]) # mean edit level of gene obtained
}


#Gene_sym <- df_m[df_m$hgnc_symbol == gene,"genes"]
Gene_sym <- names(which(gene_to_check == gene ))

rpkm_gene <- rpkm[Gene_sym,]
row.names(rpkm_gene) <- "GeneExpr"
rpkm_gene["DevelopmentStage", ] <- sapply(strsplit(colnames(rpkm_gene), split = "[.]"), "[[", 2) # expression level of gene obtained


library(ggplot2)
library(cowplot)

p <- list()
for (i in names(edit_level_imputed)){
    if (i == "Liver") {
    edit_to_plot <- data.frame(edit_level_imputed_sub_colMean[[i]])
    colnames(edit_to_plot) <- "EditLevel"
    for (j in row.names(edit_to_plot)){
        edit_to_plot[j, "DevelopmentStage"] <- OrganDev_meta[paste0("Library",OrganDev_meta$library.ID) == j, "Developmental.stage"]
    }
    if (i == "Forebrain"){
        expr_to_plot <- rpkm_gene[,grep("Brain", colnames(rpkm))]
    } else if (i == "Hindbrain"){
        expr_to_plot <- rpkm_gene[,grep("Cerebellum", colnames(rpkm))]
    } else {
        expr_to_plot <- rpkm_gene[,grep(i, colnames(rpkm))]
    }
    edit_to_plot <- dplyr::left_join(edit_to_plot, DevelopmentStage.data, by= "DevelopmentStage")
    expr_to_plot <- dplyr::left_join(data.frame(t(expr_to_plot)), DevelopmentStage.data, by= "DevelopmentStage")
    expr_to_plot$GeneExpr <- log2(as.numeric(expr_to_plot$GeneExpr))

    p[[paste0(i, "_Edit")]] <- ggplot(edit_to_plot, aes(x = DevelopmentStage_num, y = EditLevel)) + geom_point(size=5, color = "darkblue") + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + 
  theme_classic() + ggtitle(paste0("Mean Editing Level in ",i)) + coord_cartesian(xlim = c(1,25)) + 
theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)

    p[[paste0(i, "_Expr")]] <- ggplot(expr_to_plot, aes(x = DevelopmentStage_num, y = GeneExpr)) + geom_point(size=5, color = "darkblue") + 
    geom_smooth(linewidth=1.5, se = F) + 
  scale_x_continuous(
    breaks=1:length(lat.order.levels),
    labels=lat.order.levels
  ) + 
  theme_classic() + ggtitle(paste0("Gene Expression Level in ",i)) + coord_cartesian(xlim = c(1,25)) + 
theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16)) + geom_vline(xintercept = DevelopmentStage.data[DevelopmentStage.data$DevelopmentStage %in% c("11wpc","newborn","oldTeenager"),"DevelopmentStage_num"], linetype="dotdash") #+ ylim(min_exp, max_exp)

edit_expr_to_plot <- dplyr::inner_join(edit_to_plot, expr_to_plot, by= "DevelopmentStage")

    
    cor_results <- cor.test(edit_expr_to_plot$GeneExpr, edit_expr_to_plot$EditLevel, method = "pearson")
    annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_results$estimate, digits = 3), "\n ", "p = ", format(cor_results$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust 

    p[[paste0(i, "_Corr")]] <- ggplot(edit_expr_to_plot, aes(x=GeneExpr, y=EditLevel)) + geom_point(size = 5, color = "darkblue") + 
      theme_classic() + geom_smooth(method = 'lm', fullrange = T, se = FALSE) + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
ggtitle("Correlating Editing Rate Change to Expression Change") +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 20), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 18),
    axis.title=element_text(size=14,face="bold")) + 
    xlab("log2RPKM") + ylab("Average Editing Rate")
}
}
pdf(paste0("Figure_OrganDevelopment_R1/GeneSpecificPlots/",gene,"_editingplusgeneexpr_Liver_intronic_20240404.pdf"), width = 30, height = 8)
print(plot_grid(plotlist = p, nrow = 1, ncol = 3))
dev.off()

}

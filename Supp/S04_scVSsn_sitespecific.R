# read in sc and sn site specific editing
# carry out a bunch of qc
# benchmark # of sites by repeat and genic location for each cell type

# v2 adds in some normalization: norm by cell count, norm by reads, norm by cell count & reads

# for i in /media/NAS1/Sammi_NAS1/RNA_editing_Archive/dedupped_bam/dedupped_bam_H1_20230102/*.bam; do echo $i >> /media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_H1_20230102_reads.txt; samtools view -c $i >> /media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_H1_20230102_reads.txt ; done
# for i in /media/NAS1/Sammi_NAS1/RNA_editing_Archive/dedupped_bam/dedupped_bam_TS_20230125/*.bam; do echo $i >> /media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_TS_20230125_reads.txt; samtools view -c $i >> /media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_TS_20230125_reads.txt ; done


## reference table
# dbsnp filter
REDIref <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_GRCh38.txt", header = T, sep = "\t")
REDIref$Chrom_Pos <- paste0(REDIref$Region,"_",REDIref$Position)
head(REDIref)
dim(REDIref) # [1] 15681871       27
dim(REDIref[REDIref$dbsnp == "-",]) # [1] 13569494       27
SNP_sites <- REDIref[REDIref$dbsnp != "-","Chrom_Pos"]; head(SNP_sites)
ACALAP_sites <- REDIref[grepl("^AC|^AL|^AP", REDIref$Gene.wgEncodeGencodeBasicV34),"Chrom_Pos"]; head(ACALAP_sites); length(ACALAP_sites)
#length(REDIref[grep("AC[[:digit:]]|AL[[:digit:]]|AP[[:digit:]]", REDIref$Gene.wgEncodeGencodeBasicV34),"Chrom_Pos"])
MITO_sites <- REDIref[REDIref$Region == "GRCh38_chrM","Chrom_Pos"]; head(MITO_sites)


## read in sc site specific editing
input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_H1_v2/"
edit_table <- list()

all_edit_table <- list.files(input_dir)

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$CellType <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
}

length(edit_table) # 59
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sc <- edit_table
#edit_table_sc <- do.call(rbind, edit_table_sc)
#dim(edit_table_sc)
#write.table(edit_table_sc,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sc_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sc_noACALAP <- edit_table_sc[!grepl("^AC|^AL|^AP", edit_table_sc$Gene),]
#write.table(edit_table_sc_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sc_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

length(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))

# just drop Hematopoietic & Kidney-Prog, only keep samples with at least 10 cells in all teratoma
region_order <-  c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic")
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
                    "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle")
#setdiff(celltype_order,unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))

sc_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/WT_4ter_annotv1_cellcounttable.txt",sep = ",", header = 1, quote = "")
sc_read_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_H1_20230102_reads.txt",sep = ",", header = F, quote = "")
sc_read_count <- as.data.frame(matrix(sc_read_count$V1, ncol = 2, byrow = TRUE))
colnames(sc_read_count) <- c("Sample","Reads")
sc_read_count$Sample <- sapply(strsplit(sc_read_count$Sample, split = "/"), "[[", 8)
sc_read_count$Sample <- sapply(strsplit(sc_read_count$Sample, split = "_dedup"), "[[", 1)

site_count_table <- data.frame(region = rep(region_order, length(celltype_order)),
             mean_site_counts = NA, sd_site_counts = NA, n = NA, mean_site_counts_normToCellCount = NA, sd_site_counts_normToCellCount = NA, 
             mean_site_counts_normToReadCount = NA, sd_site_counts_normToReadCount = NA, mean_site_counts_normToCellCountReadCount = NA, sd_site_counts_normToCellCountReadCount = NA,
             celltype = rep(celltype_order, each = length(region_order)))

 for (j in 1:nrow(site_count_table)){
    edit_table_list_tmp <- edit_table[sapply(strsplit(names(edit_table), split = "_"), "[[", 1) == site_count_table[j,"celltype"]]
    sapply(edit_table_list_tmp,dim)
    site_count_tmp <- data.frame(region = rep(region_order, 4),
             site_counts = NA, cell_counts = NA, read_counts = NA,
             sample = rep(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 2)), each = length(region_order)))
    for (k in 1:nrow(site_count_tmp)){
        site_count_tmp[k,"cell_counts"] <- sc_cell_count[site_count_table[j,"celltype"], site_count_tmp[k,"sample"]]
        site_count_tmp[k,"read_counts"] <- as.numeric(sc_read_count[sc_read_count$Sample == paste0(site_count_table[j,"celltype"], "_", site_count_tmp[k,"sample"]), "Reads"])
        if (site_count_tmp[k,"region"] %in% names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)[sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)), split = "_"), "[[",2) == site_count_tmp[k,"sample"]][[1]])){
        site_count_tmp[k,"site_counts"] <- lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)[sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)), split = "_"), "[[",2) == site_count_tmp[k,"sample"]][[1]][[site_count_tmp[k,"region"]]]
        }
    }
    
     site_count_table[j, "mean_site_counts"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"], na.rm = T)
     site_count_table[j, "n"] <- length(which(!is.na(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"])))
     site_count_table[j, "sd_site_counts"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"], na.rm = F)
     site_count_table[j, "sd_site_counts_normToCellCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToReadCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCountReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"]), na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCountReadCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"]), na.rm = T) / sqrt(site_count_table[j, "n"])
 }
site_count_table

site_count_table_sc <- site_count_table
site_count_table_sc$protocol <- "sc"


#  for (j in 1:nrow(site_count_table)){
#     edit_table_list_tmp <- edit_table[sapply(strsplit(names(edit_table), split = "_"), "[[", 1) == site_count_table[j,"celltype"]]
#     site_count_tmp <- sapply(sapply(edit_table_list_tmp, "[[", "Region"),table)
#     site_count_table[j, "mean_site_counts"] <-  mean(site_count_tmp[site_count_table[j,"region"],], na.rm = T)
#     site_count_table[j, "n"] <- length(which(!is.na(site_count_tmp[site_count_table[j,"region"],])))
#     site_count_table[j, "sd_site_counts"] <- sd(site_count_tmp[site_count_table[j,"region"],], na.rm = T) / sqrt(site_count_table[j, "n"])
 
#  }

# sample_cutoff = 0.6
# editlevel_cutoff = 0.05

# sample_count <- length(table(edit_table_sc_noACALAP$Teratoma))
# sample_count
# highquality_sites_coverage <- names(table(edit_table_sc_noACALAP$Chrom_Pos)[which(table(edit_table_sc_noACALAP$Chrom_Pos) >= sample_cutoff*sample_count)])
# length(highquality_sites_coverage) # 179832

# highquality_sites_minediting_table <- aggregate(edit_table_sc_noACALAP$EditLevel, list(edit_table_sc_noACALAP$Chrom_Pos), FUN=mean)
# head(highquality_sites_minediting_table)
# highquality_sites_minediting <- highquality_sites_minediting_table[highquality_sites_minediting_table$x >= editlevel_cutoff, "Group.1"]
# length(highquality_sites_minediting) # 185396

# length(intersect(highquality_sites_coverage, highquality_sites_minediting)) # 5433
# highquality_sites_sc <- intersect(highquality_sites_coverage, highquality_sites_minediting)

# edit_table_sc_noACALAP_qc <- edit_table_sc_noACALAP[edit_table_sc_noACALAP$Chrom_Pos %in% highquality_sites_sc, ]
# dim(edit_table_sc_noACALAP_qc)
# summary(edit_table_sc_noACALAP_qc)
# write.table(edit_table_sc_noACALAP_qc,"Figure1_WildtypeH1_R1/20240307_scVSsn/hq_sites_sample_edit_table_sc_minEdit3_totEdit5_noM_noSNP_noACALAP_qc_60per.txt", sep = ",", col.names = T, row.names = T, quote = F)



## read in sn site specific editing, keep wk10 only

input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_TS_v1/"
edit_table <- list()

all_edit_table <- list.files(input_dir)
all_edit_table <- all_edit_table[sapply(strsplit(all_edit_table, split = "_"),"[[",2) == "wk10"]

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    if (dim(edit_table_tmp)[1] == 0){
        print("empty file ignored")
    }
    else {
    edit_table_tmp$Age <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 3)
    edit_table_tmp$Section <- sapply(strsplit(i, split = "_"), "[[", 4)
    edit_table_tmp$CellType <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
    }
}

length(edit_table) # 83
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sn <- edit_table
#edit_table_sn <- do.call(rbind, edit_table_sn)
#dim(edit_table_sn)
#write.table(edit_table_sn,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sn_noACALAP <- edit_table_sn[!grepl("^AC|^AL|^AP", edit_table_sn$Gene),]
#write.table(edit_table_sn_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

# drop tables that are empty after filtering
edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

length(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))

# dropping cell types including "Pancreatic-Prog","Definitive-Endoderm","Hematopoietic","Kidney-Prog"
region_order <-  c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic")
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
                    "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle")
setdiff(celltype_order,unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))
setdiff(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)), celltype_order)

# right now considering each teratoma and each experiment individually, this assumption may not hold true

sn_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_celltype_dist_col_20230130.txt",sep = ",", header = 1, quote = "")
sn_cell_count <- sn_cell_count[,sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",1) == "wk10"]
sn_read_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_TS_20230125_reads.txt",sep = ",", header = F, quote = "")
sn_read_count <- as.data.frame(matrix(sn_read_count$V1, ncol = 2, byrow = TRUE))
colnames(sn_read_count) <- c("Sample","Reads")
sn_read_count$Sample <- sapply(strsplit(sn_read_count$Sample, split = "/"), "[[", 8)
sn_read_count$Sample <- sapply(strsplit(sn_read_count$Sample, split = "_dedup"), "[[", 1)
sn_read_count <- sn_read_count[sapply(strsplit(sn_read_count$Sample, split = "_"),"[[",2) == "wk10", ]

site_count_table <- data.frame(region = rep(region_order, length(celltype_order)),
             mean_site_counts = NA, sd_site_counts = NA, n = NA, mean_site_counts_normToCellCount = NA, sd_site_counts_normToCellCount = NA, 
             mean_site_counts_normToReadCount = NA, sd_site_counts_normToReadCount = NA, mean_site_counts_normToCellCountReadCount = NA, sd_site_counts_normToCellCountReadCount = NA,
             celltype = rep(celltype_order, each = length(region_order)))

 for (j in 1:nrow(site_count_table)){
    edit_table_list_tmp <- edit_table[sapply(strsplit(names(edit_table), split = "_"), "[[", 1) == site_count_table[j,"celltype"]]
    sapply(edit_table_list_tmp,dim)
    site_count_tmp <- data.frame(region = rep(region_order, length(unique(paste0(sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 3), "_", sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 4))))),
             site_counts = NA, cell_counts = NA, read_counts = NA,
             sample = rep(unique(paste0(sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 3), "_", sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 4))), each = length(region_order)))
    for (k in 1:nrow(site_count_tmp)){
        sample_tmp <- paste0(sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)), split = "_"), "[[",3), "_", sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)), split = "_"), "[[",4))
        site_count_tmp[k,"cell_counts"] <- sn_cell_count[site_count_table[j,"celltype"], paste0(sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",2), "_", sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",3)) == site_count_tmp[k,"sample"]]
        site_count_tmp[k,"read_counts"] <- as.numeric(sn_read_count[sn_read_count$Sample == paste0(site_count_table[j,"celltype"], "_wk10_", site_count_tmp[k,"sample"]), "Reads"])
        if (site_count_tmp[k,"region"] %in% names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)[sample_tmp == site_count_tmp[k,"sample"]][[1]])){
        site_count_tmp[k,"site_counts"] <- lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)[sample_tmp == site_count_tmp[k,"sample"]][[1]][[site_count_tmp[k,"region"]]]
        }
    }
    
     site_count_table[j, "mean_site_counts"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"], na.rm = T)
     site_count_table[j, "n"] <- length(which(!is.na(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"])))
     site_count_table[j, "sd_site_counts"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])
 
     site_count_table[j, "mean_site_counts_normToReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToReadCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCountReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"]), na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCountReadCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"]), na.rm = T) / sqrt(site_count_table[j, "n"])
 }
site_count_table

site_count_table_sn <- site_count_table

site_count_table_sn$protocol <- "sn"

site_count_table_sc <- site_count_table_sc[site_count_table_sc$celltype %in% celltype_order, ]
site_count_table_cb <- rbind(site_count_table_sc, site_count_table_sn)

site_count_table_cb$region <- factor(site_count_table_cb$region, levels = c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic"))
site_count_table_cb$celltype <- factor(site_count_table_cb$celltype, levels = celltype_order)

library(ggplot2)
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_geniclocation_stackedbar.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts, fill=region)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Genic Locations of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Stage") + ylab("Number of RNA Editing Sites / Sample")
dev.off()


pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_geniclocation_stackedbar_normbyCellCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToCellCount, fill=region)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Genic Locations of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Stage") + ylab("Number of RNA Editing Sites / Cell")
dev.off()


pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_geniclocation_stackedbar_normbyReadCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToReadCount, fill=region)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Genic Locations of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Stage") + ylab("Number of RNA Editing Sites / Read")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_geniclocation_stackedbar_normbyCellCountsReadCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToCellCountReadCount, fill=region)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Genic Locations of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Stage") + ylab("Number of RNA Editing Sites per Cell per Read")
dev.off()


### drop multiome libraries

input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_TS_v1/"
edit_table <- list()

all_edit_table <- list.files(input_dir)
all_edit_table <- all_edit_table[sapply(strsplit(all_edit_table, split = "_"),"[[",2) == "wk10"]
all_edit_table <- all_edit_table[sapply(strsplit(all_edit_table, split = "_"),"[[",4) %in% c(1,2)]

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    if (dim(edit_table_tmp)[1] == 0){
        print("empty file ignored")
    }
    else {
    edit_table_tmp$Age <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 3)
    edit_table_tmp$Section <- sapply(strsplit(i, split = "_"), "[[", 4)
    edit_table_tmp$CellType <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
    }
}

length(edit_table) # 83
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sn <- edit_table
#edit_table_sn <- do.call(rbind, edit_table_sn)
#dim(edit_table_sn)
#write.table(edit_table_sn,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sn_noACALAP <- edit_table_sn[!grepl("^AC|^AL|^AP", edit_table_sn$Gene),]
#write.table(edit_table_sn_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

# drop tables that are empty after filtering
edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

length(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))

# dropping cell types including "Pancreatic-Prog","Definitive-Endoderm","Hematopoietic","Kidney-Prog"
region_order <-  c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic")
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
                    "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle")
setdiff(celltype_order,unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))
setdiff(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)), celltype_order)

# right now considering each teratoma and each experiment individually, this assumption may not hold true

sn_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_celltype_dist_col_20230130.txt",sep = ",", header = 1, quote = "")
sn_cell_count <- sn_cell_count[,sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",1) == "wk10"]
sn_cell_count <- sn_cell_count[,sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",3) %in% c(1,2)]
sn_read_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_TS_20230125_reads.txt",sep = ",", header = F, quote = "")
sn_read_count <- as.data.frame(matrix(sn_read_count$V1, ncol = 2, byrow = TRUE))
colnames(sn_read_count) <- c("Sample","Reads")
sn_read_count$Sample <- sapply(strsplit(sn_read_count$Sample, split = "/"), "[[", 8)
sn_read_count$Sample <- sapply(strsplit(sn_read_count$Sample, split = "_dedup"), "[[", 1)
sn_read_count <- sn_read_count[sapply(strsplit(sn_read_count$Sample, split = "_"),"[[",2) == "wk10", ]
sn_read_count <- sn_read_count[sapply(strsplit(sn_read_count$Sample, split = "_"),"[[",4) %in% c(1,2), ]


site_count_table <- data.frame(region = rep(region_order, length(celltype_order)),
             mean_site_counts = NA, sd_site_counts = NA, n = NA, mean_site_counts_normToCellCount = NA, sd_site_counts_normToCellCount = NA, 
             mean_site_counts_normToReadCount = NA, sd_site_counts_normToReadCount = NA, mean_site_counts_normToCellCountReadCount = NA, sd_site_counts_normToCellCountReadCount = NA,
             celltype = rep(celltype_order, each = length(region_order)))

 for (j in 1:nrow(site_count_table)){
    edit_table_list_tmp <- edit_table[sapply(strsplit(names(edit_table), split = "_"), "[[", 1) == site_count_table[j,"celltype"]]
    sapply(edit_table_list_tmp,dim)
    site_count_tmp <- data.frame(region = rep(region_order, length(unique(paste0(sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 3), "_", sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 4))))),
             site_counts = NA, cell_counts = NA, read_counts = NA,
             sample = rep(unique(paste0(sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 3), "_", sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 4))), each = length(region_order)))
    for (k in 1:nrow(site_count_tmp)){
        sample_tmp <- paste0(sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)), split = "_"), "[[",3), "_", sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)), split = "_"), "[[",4))
        site_count_tmp[k,"cell_counts"] <- sn_cell_count[site_count_table[j,"celltype"], paste0(sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",2), "_", sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",3)) == site_count_tmp[k,"sample"]]
        site_count_tmp[k,"read_counts"] <- as.numeric(sn_read_count[sn_read_count$Sample == paste0(site_count_table[j,"celltype"], "_wk10_", site_count_tmp[k,"sample"]), "Reads"])
        if (site_count_tmp[k,"region"] %in% names(lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)[sample_tmp == site_count_tmp[k,"sample"]][[1]])){
        site_count_tmp[k,"site_counts"] <- lapply(sapply(edit_table_list_tmp, "[[", "Region"),table)[sample_tmp == site_count_tmp[k,"sample"]][[1]][[site_count_tmp[k,"region"]]]
        }
    }
    
     site_count_table[j, "mean_site_counts"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"], na.rm = T)
     site_count_table[j, "n"] <- length(which(!is.na(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"])))
     site_count_table[j, "sd_site_counts"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])
 
     site_count_table[j, "mean_site_counts_normToReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToReadCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCountReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"]), na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCountReadCount"] <- sd(site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$region == site_count_table[j,"region"]),"read_counts"]), na.rm = T) / sqrt(site_count_table[j, "n"])
 }
site_count_table

site_count_table_sn <- site_count_table

site_count_table_sn$protocol <- "sn"

site_count_table_sc <- site_count_table_sc[site_count_table_sc$celltype %in% celltype_order, ]
site_count_table_cb <- rbind(site_count_table_sc, site_count_table_sn)

site_count_table_cb$region <- factor(site_count_table_cb$region, levels = c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic"))
site_count_table_cb$celltype <- factor(site_count_table_cb$celltype, levels = celltype_order)

library(ggplot2)
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_geniclocation_stackedbar_dropmultiome.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts, fill=region)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Genic Locations of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Stage") + ylab("Number of RNA Editing Sites / Sample")
dev.off()


pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_geniclocation_stackedbar_dropmultiome_normbyCellCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToCellCount, fill=region)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Genic Locations of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Stage") + ylab("Number of RNA Editing Sites / Cell")
dev.off()


pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_geniclocation_stackedbar_dropmultiome_normbyReadCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToReadCount, fill=region)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Genic Locations of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Stage") + ylab("Number of RNA Editing Sites / Read")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_geniclocation_stackedbar_dropmultiome_normbyCellCountsReadCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToCellCountReadCount, fill=region)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Genic Locations of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Stage") + ylab("Number of RNA Editing Sites per Cell per Read")
dev.off()








## read in sc site specific editing
input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_H1_v2/"
edit_table <- list()

all_edit_table <- list.files(input_dir)

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$CellType <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
}

length(edit_table) # 59
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sc <- edit_table
#edit_table_sc <- do.call(rbind, edit_table_sc)
#dim(edit_table_sc)
#write.table(edit_table_sc,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sc_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sc_noACALAP <- edit_table_sc[!grepl("^AC|^AL|^AP", edit_table_sc$Gene),]
#write.table(edit_table_sc_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sc_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

length(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))

# just drop Hematopoietic & Kidney-Prog, only keep samples with at least 10 cells in all teratoma
alu_order <-  c("ALU","REP","NONREP")
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
                    "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle")
#setdiff(celltype_order,unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))

sc_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/WT_4ter_annotv1_cellcounttable.txt",sep = ",", header = 1, quote = "")
sc_read_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_H1_20230102_reads.txt",sep = ",", header = F, quote = "")
sc_read_count <- as.data.frame(matrix(sc_read_count$V1, ncol = 2, byrow = TRUE))
colnames(sc_read_count) <- c("Sample","Reads")
sc_read_count$Sample <- sapply(strsplit(sc_read_count$Sample, split = "/"), "[[", 8)
sc_read_count$Sample <- sapply(strsplit(sc_read_count$Sample, split = "_dedup"), "[[", 1)

site_count_table <- data.frame(alu = rep(alu_order, length(celltype_order)),
             mean_site_counts = NA, sd_site_counts = NA, n = NA, mean_site_counts_normToCellCount = NA, sd_site_counts_normToCellCount = NA, 
             mean_site_counts_normToReadCount = NA, sd_site_counts_normToReadCount = NA,  mean_site_counts_normToCellCountReadCount = NA, sd_site_counts_normToCellCountReadCount = NA,
             celltype = rep(celltype_order, each = length(alu_order)))

 for (j in 1:nrow(site_count_table)){
    edit_table_list_tmp <- edit_table[sapply(strsplit(names(edit_table), split = "_"), "[[", 1) == site_count_table[j,"celltype"]]
    sapply(edit_table_list_tmp,dim)
    site_count_tmp <- data.frame(alu = rep(alu_order, 4),
             site_counts = NA, cell_counts = NA, read_counts = NA,
             sample = rep(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 2)), each = length(alu_order)))
    for (k in 1:nrow(site_count_tmp)){
        site_count_tmp[k,"cell_counts"] <- sc_cell_count[site_count_table[j,"celltype"], site_count_tmp[k,"sample"]]
        site_count_tmp[k,"read_counts"] <- as.numeric(sc_read_count[sc_read_count$Sample == paste0(site_count_table[j,"celltype"], "_", site_count_tmp[k,"sample"]), "Reads"])
        if (site_count_tmp[k,"alu"] %in% names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)[sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)), split = "_"), "[[",2) == site_count_tmp[k,"sample"]][[1]])){
        site_count_tmp[k,"site_counts"] <- lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)[sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)), split = "_"), "[[",2) == site_count_tmp[k,"sample"]][[1]][[site_count_tmp[k,"alu"]]]
        }
    }
    
     site_count_table[j, "mean_site_counts"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"], na.rm = T)
     site_count_table[j, "n"] <- length(which(!is.na(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"])))
     site_count_table[j, "sd_site_counts"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])
 
     site_count_table[j, "mean_site_counts_normToReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToReadCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCountReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"]), na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCountReadCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"]), na.rm = T) / sqrt(site_count_table[j, "n"])

 }
site_count_table

site_count_table_sc <- site_count_table
site_count_table_sc$protocol <- "sc"


## read in sn site specific editing, keep wk10 only

input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_TS_v1/"
edit_table <- list()

all_edit_table <- list.files(input_dir)
all_edit_table <- all_edit_table[sapply(strsplit(all_edit_table, split = "_"),"[[",2) == "wk10"]

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    if (dim(edit_table_tmp)[1] == 0){
        print("empty file ignored")
    }
    else {
    edit_table_tmp$Age <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 3)
    edit_table_tmp$Section <- sapply(strsplit(i, split = "_"), "[[", 4)
    edit_table_tmp$CellType <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
    }
}

length(edit_table) # 83
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sn <- edit_table
#edit_table_sn <- do.call(rbind, edit_table_sn)
#dim(edit_table_sn)
#write.table(edit_table_sn,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sn_noACALAP <- edit_table_sn[!grepl("^AC|^AL|^AP", edit_table_sn$Gene),]
#write.table(edit_table_sn_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

# drop tables that are empty after filtering
edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

length(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))

# dropping cell types including "Pancreatic-Prog","Definitive-Endoderm","Hematopoietic","Kidney-Prog"
alu_order <-  c("ALU","REP","NONREP")
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
                    "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle")
setdiff(celltype_order,unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))
setdiff(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)), celltype_order)

# right now considering each teratoma and each experiment individually, this assumption may not hold true

sn_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_celltype_dist_col_20230130.txt",sep = ",", header = 1, quote = "")
sn_cell_count <- sn_cell_count[,sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",1) == "wk10"]
sn_read_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_TS_20230125_reads.txt",sep = ",", header = F, quote = "")
sn_read_count <- as.data.frame(matrix(sn_read_count$V1, ncol = 2, byrow = TRUE))
colnames(sn_read_count) <- c("Sample","Reads")
sn_read_count$Sample <- sapply(strsplit(sn_read_count$Sample, split = "/"), "[[", 8)
sn_read_count$Sample <- sapply(strsplit(sn_read_count$Sample, split = "_dedup"), "[[", 1)
sn_read_count <- sn_read_count[sapply(strsplit(sn_read_count$Sample, split = "_"),"[[",2) == "wk10", ]

site_count_table <- data.frame(alu = rep(alu_order, length(celltype_order)),
             mean_site_counts = NA, sd_site_counts = NA, n = NA, mean_site_counts_normToCellCount = NA, sd_site_counts_normToCellCount = NA,
             mean_site_counts_normToReadCount = NA, sd_site_counts_normToReadCount = NA, mean_site_counts_normToCellCountReadCount = NA, sd_site_counts_normToCellCountReadCount = NA,
             celltype = rep(celltype_order, each = length(alu_order)))

 for (j in 1:nrow(site_count_table)){
    edit_table_list_tmp <- edit_table[sapply(strsplit(names(edit_table), split = "_"), "[[", 1) == site_count_table[j,"celltype"]]
    sapply(edit_table_list_tmp,dim)
    site_count_tmp <- data.frame(alu = rep(alu_order, length(unique(paste0(sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 3), "_", sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 4))))),
             site_counts = NA, cell_counts = NA,
             sample = rep(unique(paste0(sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 3), "_", sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 4))), each = length(alu_order)))
    for (k in 1:nrow(site_count_tmp)){
        sample_tmp <- paste0(sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)), split = "_"), "[[",3), "_", sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)), split = "_"), "[[",4))
        site_count_tmp[k,"cell_counts"] <- sn_cell_count[site_count_table[j,"celltype"], paste0(sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",2), "_", sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",3)) == site_count_tmp[k,"sample"]]
        site_count_tmp[k,"read_counts"] <- as.numeric(sn_read_count[sn_read_count$Sample == paste0(site_count_table[j,"celltype"], "_wk10_", site_count_tmp[k,"sample"]), "Reads"])
        if (site_count_tmp[k,"alu"] %in% names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)[sample_tmp == site_count_tmp[k,"sample"]][[1]])){
        site_count_tmp[k,"site_counts"] <- lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)[sample_tmp == site_count_tmp[k,"sample"]][[1]][[site_count_tmp[k,"alu"]]]
        }
    }
    
     site_count_table[j, "mean_site_counts"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"], na.rm = T)
     site_count_table[j, "n"] <- length(which(!is.na(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"])))
     site_count_table[j, "sd_site_counts"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])
 
     site_count_table[j, "mean_site_counts_normToCellCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToReadCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCountReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"]), na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCountReadCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"]), na.rm = T) / sqrt(site_count_table[j, "n"])

 }
site_count_table

site_count_table_sn <- site_count_table

site_count_table_sn$protocol <- "sn"

site_count_table_sc <- site_count_table_sc[site_count_table_sc$celltype %in% celltype_order, ]
site_count_table_cb <- rbind(site_count_table_sc, site_count_table_sn)

site_count_table_cb$alu <- factor(site_count_table_cb$alu, levels = c("ALU","REP","NONREP"))
site_count_table_cb$celltype <- factor(site_count_table_cb$celltype, levels = celltype_order)

library(ggplot2)
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_repeat_stackedbar.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts, fill=alu)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Repeat Categories of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Protocol") + ylab("Number of RNA Editing Sites / Sample")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_repeat_stackedbar_normbyCellCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToCellCount, fill=alu)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Repeat Categories of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Protocol") + ylab("Number of RNA Editing Sites / Cell")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_repeat_stackedbar_normbyReadCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToReadCount, fill=alu)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Repeat Categories of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Protocol") + ylab("Number of RNA Editing Sites / Read")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_repeat_stackedbar_normbyCellCountReadCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToCellCountReadCount, fill=alu)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Repeat Categories of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Protocol") + ylab("Number of RNA Editing Sites per Cell per Read")
dev.off()

## drop multiome libraries

input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_TS_v1/"
edit_table <- list()

all_edit_table <- list.files(input_dir)
all_edit_table <- all_edit_table[sapply(strsplit(all_edit_table, split = "_"),"[[",2) == "wk10"]
all_edit_table <- all_edit_table[sapply(strsplit(all_edit_table, split = "_"),"[[",4) %in% c(1,2)]

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    if (dim(edit_table_tmp)[1] == 0){
        print("empty file ignored")
    }
    else {
    edit_table_tmp$Age <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 3)
    edit_table_tmp$Section <- sapply(strsplit(i, split = "_"), "[[", 4)
    edit_table_tmp$CellType <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
    }
}

length(edit_table) # 83
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sn <- edit_table
#edit_table_sn <- do.call(rbind, edit_table_sn)
#dim(edit_table_sn)
#write.table(edit_table_sn,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sn_noACALAP <- edit_table_sn[!grepl("^AC|^AL|^AP", edit_table_sn$Gene),]
#write.table(edit_table_sn_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

# drop tables that are empty after filtering
edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

length(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))

# dropping cell types including "Pancreatic-Prog","Definitive-Endoderm","Hematopoietic","Kidney-Prog"
alu_order <-  c("ALU","REP","NONREP")
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
                    "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle")
setdiff(celltype_order,unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)))
setdiff(unique(sapply(strsplit(names(edit_table), split = "_"), "[[", 1)), celltype_order)

# right now considering each teratoma and each experiment individually, this assumption may not hold true

sn_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_celltype_dist_col_20230130.txt",sep = ",", header = 1, quote = "")
sn_cell_count <- sn_cell_count[,sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",1) == "wk10"]
sn_cell_count <- sn_cell_count[,sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",3) %in% c(1,2)]
sn_read_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/20240307_scVSsn/dedupped_bam_TS_20230125_reads.txt",sep = ",", header = F, quote = "")
sn_read_count <- as.data.frame(matrix(sn_read_count$V1, ncol = 2, byrow = TRUE))
colnames(sn_read_count) <- c("Sample","Reads")
sn_read_count$Sample <- sapply(strsplit(sn_read_count$Sample, split = "/"), "[[", 8)
sn_read_count$Sample <- sapply(strsplit(sn_read_count$Sample, split = "_dedup"), "[[", 1)
sn_read_count <- sn_read_count[sapply(strsplit(sn_read_count$Sample, split = "_"),"[[",2) == "wk10", ]
sn_read_count <- sn_read_count[sapply(strsplit(sn_read_count$Sample, split = "_"),"[[",4) %in% c(1,2), ]

site_count_table <- data.frame(alu = rep(alu_order, length(celltype_order)),
             mean_site_counts = NA, sd_site_counts = NA, n = NA, mean_site_counts_normToCellCount = NA, sd_site_counts_normToCellCount = NA,
             mean_site_counts_normToReadCount = NA, sd_site_counts_normToReadCount = NA, mean_site_counts_normToCellCountReadCount = NA, sd_site_counts_normToCellCountReadCount = NA,
             celltype = rep(celltype_order, each = length(alu_order)))

 for (j in 1:nrow(site_count_table)){
    edit_table_list_tmp <- edit_table[sapply(strsplit(names(edit_table), split = "_"), "[[", 1) == site_count_table[j,"celltype"]]
    sapply(edit_table_list_tmp,dim)
    site_count_tmp <- data.frame(alu = rep(alu_order, length(unique(paste0(sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 3), "_", sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 4))))),
             site_counts = NA, cell_counts = NA,
             sample = rep(unique(paste0(sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 3), "_", sapply(strsplit(names(edit_table_list_tmp), split = "_"), "[[", 4))), each = length(alu_order)))
    for (k in 1:nrow(site_count_tmp)){
        sample_tmp <- paste0(sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)), split = "_"), "[[",3), "_", sapply(strsplit(names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)), split = "_"), "[[",4))
        site_count_tmp[k,"cell_counts"] <- sn_cell_count[site_count_table[j,"celltype"], paste0(sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",2), "_", sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",3)) == site_count_tmp[k,"sample"]]
        site_count_tmp[k,"read_counts"] <- as.numeric(sn_read_count[sn_read_count$Sample == paste0(site_count_table[j,"celltype"], "_wk10_", site_count_tmp[k,"sample"]), "Reads"])
        if (site_count_tmp[k,"alu"] %in% names(lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)[sample_tmp == site_count_tmp[k,"sample"]][[1]])){
        site_count_tmp[k,"site_counts"] <- lapply(sapply(edit_table_list_tmp, "[[", "Alu"),table)[sample_tmp == site_count_tmp[k,"sample"]][[1]][[site_count_tmp[k,"alu"]]]
        }
    }
    
     site_count_table[j, "mean_site_counts"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"], na.rm = T)
     site_count_table[j, "n"] <- length(which(!is.na(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"])))
     site_count_table[j, "sd_site_counts"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])
 
     site_count_table[j, "mean_site_counts_normToCellCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"], na.rm = T)
     site_count_table[j, "sd_site_counts_normToReadCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"], na.rm = T) / sqrt(site_count_table[j, "n"])

     site_count_table[j, "mean_site_counts_normToCellCountReadCount"] <-  mean(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"]), na.rm = T)
     site_count_table[j, "sd_site_counts_normToCellCountReadCount"] <- sd(site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"site_counts"] / (site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"cell_counts"] * site_count_tmp[which(site_count_tmp$alu == site_count_table[j,"alu"]),"read_counts"]), na.rm = T) / sqrt(site_count_table[j, "n"])

 }
site_count_table

site_count_table_sn <- site_count_table

site_count_table_sn$protocol <- "sn"

site_count_table_sc <- site_count_table_sc[site_count_table_sc$celltype %in% celltype_order, ]
site_count_table_cb <- rbind(site_count_table_sc, site_count_table_sn)

site_count_table_cb$alu <- factor(site_count_table_cb$alu, levels = c("ALU","REP","NONREP"))
site_count_table_cb$celltype <- factor(site_count_table_cb$celltype, levels = celltype_order)

library(ggplot2)
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_repeat_stackedbar_dropmultiome.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts, fill=alu)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Repeat Categories of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Protocol") + ylab("Number of RNA Editing Sites / Sample")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_repeat_stackedbar_dropmultiome_normbyCellCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToCellCount, fill=alu)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Repeat Categories of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Protocol") + ylab("Number of RNA Editing Sites / Cell")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_repeat_stackedbar_dropmultiome_normbyReadCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToReadCount, fill=alu)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Repeat Categories of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Protocol") + ylab("Number of RNA Editing Sites / Read")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecific_repeat_stackedbar_dropmultiome_normbyCellCountReadCount.pdf", width = 25, height = 15)
ggplot(data=site_count_table_cb, aes(x=protocol, y=mean_site_counts_normToCellCountReadCount, fill=alu)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Repeat Categories of Sites Captured") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~celltype,scales = "free_x") +
     xlab("Protocol") + ylab("Number of RNA Editing Sites per Cell per Read")
dev.off()




#### cell count bar plot + percentage distribution plot

celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Pancreatic-Prog","Definitive-Endoderm","Gut-Epi","Hematopoietic","Adipogenic-MSC-Fib",
                     "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle", "Kidney-Prog")

sc_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/WT_4ter_annotv1_cellcounttable.txt",sep = ",", header = 1, quote = "")

sn_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_celltype_dist_col_20230130.txt",sep = ",", header = 1, quote = "")
sn_cell_count <- sn_cell_count[,sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",1) == "wk10"]


data_to_plot_sc <- reshape:::melt(as.matrix(sc_cell_count))
colnames(data_to_plot_sc) <- c("cell_type","sample","cell_count")
data_to_plot_sc$protocol <- "sc"
data_to_plot_sn <- reshape:::melt(as.matrix(sn_cell_count))
colnames(data_to_plot_sn) <- c("cell_type","sample","cell_count")
data_to_plot_sn$protocol <- "sn"

data_to_plot_cb <- rbind(data_to_plot_sc, data_to_plot_sn)
data_to_plot_cb$cell_type <- factor(data_to_plot_cb$cell_type, levels = celltype_order)

cell_count_table <- data.frame(cell_type = rep(names(table(data_to_plot_cb$cell_type)), 2),
             mean_cell_counts = NA, sd_cell_counts = NA, n = NA, protocol = rep(c("sc","sn"), each = length(table(data_to_plot_cb$cell_type))))

for (i in 1:nrow(cell_count_table)){
    cell_count_table[i,"mean_cell_counts"] <- mean(data_to_plot_cb[data_to_plot_cb$cell_type == cell_count_table[i,"cell_type"] & data_to_plot_cb$protocol ==  cell_count_table[i,"protocol"], "cell_count"])
    cell_count_table[i,"n"] <- length(data_to_plot_cb[data_to_plot_cb$cell_type == cell_count_table[i,"cell_type"] & data_to_plot_cb$protocol ==  cell_count_table[i,"protocol"], "cell_count"])
    cell_count_table[i,"sd_cell_counts"] <- sd(data_to_plot_cb[data_to_plot_cb$cell_type == cell_count_table[i,"cell_type"] & data_to_plot_cb$protocol ==  cell_count_table[i,"protocol"], "cell_count"]) / sqrt(cell_count_table[i,"n"])

}
cell_count_table
cell_count_table$cell_type <- factor(cell_count_table$cell_type, levels = celltype_order)


library(ggplot2)
ggplot(data=cell_count_table, aes(x=protocol, y=mean_cell_counts, fill=protocol)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Cell Count from single cell RNA-seq vs. single nucleus RNA-seq") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~cell_type,scales = "free_x") +
     xlab("Protocol") + ylab("Number of Cells per Library")

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_CellCount.pdf", width = 27, height = 15)
ggplot(data=cell_count_table, aes(x=protocol, y=mean_cell_counts, fill=protocol)) + scale_fill_brewer(palette="Set3") + geom_bar(stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=mean_cell_counts-sd_cell_counts, ymax=mean_cell_counts+sd_cell_counts), width=.2, color = 'red',
                 position=position_dodge(.9)) +
    geom_jitter(data=data_to_plot_cb, aes(protocol, cell_count), stat="identity", color = "darkgrey",  position = position_jitterdodge(jitter.height = .1, jitter.width = .2)) + 
     theme_minimal() + ggtitle("Cell Count from single cell RNA-seq vs. single nucleus RNA-seq") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~cell_type,scales = "free_x") + 
     xlab("Protocol") + ylab("Number of Cells per Library")
dev.off()





sc_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure1_WildtypeH1_R1/WT_4ter_annotv1_cellcounttable.txt",sep = ",", header = 1, quote = "")

sn_cell_count <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure2_TimeSeries_R1/TS_all_RNA_celltype_dist_col_20230130.txt",sep = ",", header = 1, quote = "")
sn_cell_count <- sn_cell_count[,sapply(strsplit(colnames(sn_cell_count), split = "_"),"[[",1) == "wk10"]

sc_cell_count <- sweep(sc_cell_count,2,colSums(sc_cell_count),`/`) * 100
sn_cell_count <- sweep(sn_cell_count,2,colSums(sn_cell_count),`/`) * 100


data_to_plot_sc <- reshape:::melt(as.matrix(sc_cell_count))
colnames(data_to_plot_sc) <- c("cell_type","sample","cell_count")
data_to_plot_sc$protocol <- "sc"
data_to_plot_sn <- reshape:::melt(as.matrix(sn_cell_count))
colnames(data_to_plot_sn) <- c("cell_type","sample","cell_count")
data_to_plot_sn$protocol <- "sn"

data_to_plot_cb <- rbind(data_to_plot_sc, data_to_plot_sn)
data_to_plot_cb$cell_type <- factor(data_to_plot_cb$cell_type, levels = celltype_order)

cell_count_table <- data.frame(cell_type = rep(names(table(data_to_plot_cb$cell_type)), 2),
             mean_cell_counts = NA, sd_cell_counts = NA, n = NA, protocol = rep(c("sc","sn"), each = length(table(data_to_plot_cb$cell_type))))

for (i in 1:nrow(cell_count_table)){
    cell_count_table[i,"mean_cell_counts"] <- mean(data_to_plot_cb[data_to_plot_cb$cell_type == cell_count_table[i,"cell_type"] & data_to_plot_cb$protocol ==  cell_count_table[i,"protocol"], "cell_count"])
    cell_count_table[i,"n"] <- length(data_to_plot_cb[data_to_plot_cb$cell_type == cell_count_table[i,"cell_type"] & data_to_plot_cb$protocol ==  cell_count_table[i,"protocol"], "cell_count"])
    cell_count_table[i,"sd_cell_counts"] <- sd(data_to_plot_cb[data_to_plot_cb$cell_type == cell_count_table[i,"cell_type"] & data_to_plot_cb$protocol ==  cell_count_table[i,"protocol"], "cell_count"]) / sqrt(cell_count_table[i,"n"])

}
cell_count_table
cell_count_table$cell_type <- factor(cell_count_table$cell_type, levels = celltype_order)

library(ggplot2)
ggplot(data=cell_count_table, aes(x=protocol, y=mean_cell_counts, fill=protocol)) + scale_fill_brewer(palette="Set3") + 
    geom_bar(stat="identity", width = 0.9) + theme_minimal() + ggtitle("Cell Percentage Distribution from single cell RNA-seq vs. single nucleus RNA-seq") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~cell_type,scales = "free_x") +
     xlab("Protocol") + ylab("Percentage of Cells within Library")

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_PercentageWithinLibrary.pdf", width = 27, height = 15)
ggplot(data=cell_count_table, aes(x=protocol, y=mean_cell_counts, fill=protocol)) + scale_fill_brewer(palette="Set3") + geom_bar(stat="identity", width = 0.9) + 
  geom_errorbar(aes(ymin=mean_cell_counts-sd_cell_counts, ymax=mean_cell_counts+sd_cell_counts), width=.2, color = 'red',
                 position=position_dodge(.9)) +
    geom_jitter(data=data_to_plot_cb, aes(protocol, cell_count), stat="identity", color = "darkgrey",  position = position_jitterdodge(jitter.height = .1, jitter.width = .2)) + 
     theme_minimal() + ggtitle("Cell Percentage Distribution from single cell RNA-seq vs. single nucleus RNA-seq") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), axis.text.y = element_text( vjust = 0.5, hjust=1, size = 12),
     axis.title=element_text(size=14,face="bold")) +
  facet_grid(~cell_type,scales = "free_x") + 
     xlab("Protocol") + ylab("Percentage of Cells within Library")
dev.off()



### -------------------------------------------------------------------
### contrasting RNA editing levels between sc and sn 
### AEI + site specific

# AEI
sc_AEI <- read.table('git_repo/RNAEditing/Figure3/AEI_H1_v1_EditingIndex.csv', header = T, sep =',')
sn_AEI <- read.table('git_repo/RNAEditing/Figure4/AEI_TS_v1_EditingIndex.txt', header = T, sep ='\t')

head(sc_AEI)
sc_AEI$celltype <- sapply(strsplit(sc_AEI$Sample, split = "_"), "[[", 1)
sc_AEI$teratoma <- sapply(strsplit(sc_AEI$Sample, split = "_"), "[[", 2)

head(sn_AEI)
sn_AEI <- sn_AEI[sn_AEI$Age == "wk10",]

sc_vs_sn_AEI <- data.frame(celltype = intersect(names(table(sc_AEI$celltype)), names(table(sn_AEI$CellType))),
sc_AEI_mean = NA, sc_AEI_sd = NA, sn_AEI_mean = NA, sn_AEI_sd = NA)

for (i in 1:nrow(sc_vs_sn_AEI)) {
    sc_vs_sn_AEI[i, "sc_AEI_mean"] <- mean(as.numeric(sc_AEI[sc_AEI$celltype == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"]))
    sc_vs_sn_AEI[i, "sc_AEI_sd"] <- sd(as.numeric(sc_AEI[sc_AEI$celltype == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"])) / sqrt(length(sc_AEI[sc_AEI$celltype == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"]))

    sc_vs_sn_AEI[i, "sn_AEI_mean"] <- mean(as.numeric(sn_AEI[sn_AEI$CellType == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"]))
    sc_vs_sn_AEI[i, "sn_AEI_sd"] <- sd(as.numeric(sn_AEI[sn_AEI$CellType == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"])) / sqrt(length(sn_AEI[sn_AEI$CellType == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"]))

}

cor_test_result <- cor.test(sc_vs_sn_AEI$sc_AEI_mean, sc_vs_sn_AEI$sn_AEI_mean)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_test_result$estimate, digits = 3), "\n ", "p = ", format(cor_test_result$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust

# same color scheme as Fig 5C
cluster_col <- c("#0000FF","#FF0000","#00C000","#F07CAB","#000080",
                 "#AD07E3","#FF8000","#000000","#90BFF9","#C0C0FF",
                 "#D30B94","#A00000","#F2B77C","#00FF00","#94641F")
names(cluster_col) <- sc_vs_sn_AEI$celltype
cluster_col

library(ggplot2)
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_AEI_withgreytrendline.pdf", width = 10, height = 8)
ggplot(sc_vs_sn_AEI, aes(x = sc_AEI_mean, y = sn_AEI_mean, color = celltype)) + geom_point(size=5) + scale_color_manual(values = cluster_col) + 
        geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific AEI (sc vs. sn)")) + coord_cartesian(xlim = c(0.3,0.6), ylim = c(0.5, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean AEI for single cell Samples") + ylab("Mean AEI for single nucleus Samples")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_AEI_withouttrendline.pdf", width = 10, height = 8)
ggplot(sc_vs_sn_AEI, aes(x = sc_AEI_mean, y = sn_AEI_mean, color = celltype)) + geom_point(size=5) + scale_color_manual(values = cluster_col) + 
        #geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific AEI (sc vs. sn)")) + coord_cartesian(xlim = c(0.3,0.6), ylim = c(0.5, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean AEI for single cell Samples") + ylab("Mean AEI for single nucleus Samples")
dev.off()

# AEI dropping muscle and kidney-prog
sc_AEI <- read.table('git_repo/RNAEditing/Figure3/AEI_H1_v1_EditingIndex.csv', header = T, sep =',')
sn_AEI <- read.table('git_repo/RNAEditing/Figure4/AEI_TS_v1_EditingIndex.txt', header = T, sep ='\t')

head(sc_AEI)
sc_AEI$celltype <- sapply(strsplit(sc_AEI$Sample, split = "_"), "[[", 1)
sc_AEI$teratoma <- sapply(strsplit(sc_AEI$Sample, split = "_"), "[[", 2)

head(sn_AEI)
sn_AEI <- sn_AEI[sn_AEI$Age == "wk10",]

sc_vs_sn_AEI <- data.frame(celltype = intersect(names(table(sc_AEI$celltype)), names(table(sn_AEI$CellType))),
sc_AEI_mean = NA, sc_AEI_sd = NA, sn_AEI_mean = NA, sn_AEI_sd = NA)

for (i in 1:nrow(sc_vs_sn_AEI)) {
    sc_vs_sn_AEI[i, "sc_AEI_mean"] <- mean(as.numeric(sc_AEI[sc_AEI$celltype == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"]))
    sc_vs_sn_AEI[i, "sc_AEI_sd"] <- sd(as.numeric(sc_AEI[sc_AEI$celltype == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"])) / sqrt(length(sc_AEI[sc_AEI$celltype == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"]))

    sc_vs_sn_AEI[i, "sn_AEI_mean"] <- mean(as.numeric(sn_AEI[sn_AEI$CellType == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"]))
    sc_vs_sn_AEI[i, "sn_AEI_sd"] <- sd(as.numeric(sn_AEI[sn_AEI$CellType == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"])) / sqrt(length(sn_AEI[sn_AEI$CellType == sc_vs_sn_AEI[i, "celltype"], "A2GEditingIndex"]))

}
sc_vs_sn_AEI <- sc_vs_sn_AEI[!sc_vs_sn_AEI$celltype %in% c("Kidney-Prog", "Muscle"),]

cor_test_result <- cor.test(sc_vs_sn_AEI$sc_AEI_mean, sc_vs_sn_AEI$sn_AEI_mean)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_test_result$estimate, digits = 3), "\n ", "p = ", format(cor_test_result$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust

# same color scheme as Fig 5C
cluster_col <- c("#0000FF","#FF0000","#00C000","#F07CAB","#000080",
                 "#FF8000","#90BFF9","#C0C0FF",
                 "#D30B94","#A00000","#F2B77C","#00FF00","#94641F")
names(cluster_col) <- sc_vs_sn_AEI$celltype
cluster_col

library(ggplot2)
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_AEI_droppingKidneyMuscle_withgreytrendline.pdf", width = 10, height = 8)
ggplot(sc_vs_sn_AEI, aes(x = sc_AEI_mean, y = sn_AEI_mean, color = celltype)) + geom_point(size=5) + scale_color_manual(values = cluster_col) + 
        geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific AEI (sc vs. sn)")) + coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean AEI for single cell Samples") + ylab("Mean AEI for single nucleus Samples")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_AEI_droppingKidneyMuscle_withouttrendline.pdf", width = 10, height = 8)
ggplot(sc_vs_sn_AEI, aes(x = sc_AEI_mean, y = sn_AEI_mean, color = celltype)) + geom_point(size=5) + scale_color_manual(values = cluster_col) + 
        #geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific AEI (sc vs. sn)")) + coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean AEI for single cell Samples") + ylab("Mean AEI for single nucleus Samples")
dev.off()

###### site specific version series 1 by cell type
## reference table
# dbsnp filter
REDIref <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_GRCh38.txt", header = T, sep = "\t")
REDIref$Chrom_Pos <- paste0(REDIref$Region,"_",REDIref$Position)
head(REDIref)
dim(REDIref) # [1] 15681871       27
dim(REDIref[REDIref$dbsnp == "-",]) # [1] 13569494       27
SNP_sites <- REDIref[REDIref$dbsnp != "-","Chrom_Pos"]; head(SNP_sites)
ACALAP_sites <- REDIref[grepl("^AC|^AL|^AP", REDIref$Gene.wgEncodeGencodeBasicV34),"Chrom_Pos"]; head(ACALAP_sites); length(ACALAP_sites)
#length(REDIref[grep("AC[[:digit:]]|AL[[:digit:]]|AP[[:digit:]]", REDIref$Gene.wgEncodeGencodeBasicV34),"Chrom_Pos"])
MITO_sites <- REDIref[REDIref$Region == "GRCh38_chrM","Chrom_Pos"]; head(MITO_sites)


## read in sc site specific editing
input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_H1_v2/"
edit_table <- list()

all_edit_table <- list.files(input_dir)

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$CellType <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
}

length(edit_table) # 59
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sc <- edit_table
#edit_table_sc <- do.call(rbind, edit_table_sc)
#dim(edit_table_sc)
#write.table(edit_table_sc,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sc_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sc_noACALAP <- edit_table_sc[!grepl("^AC|^AL|^AP", edit_table_sc$Gene),]
#write.table(edit_table_sc_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sc_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

edit_table_sc <- edit_table

## read in sn site specific editing, keep wk10 only

input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_TS_v1/"
edit_table <- list()

all_edit_table <- list.files(input_dir)
all_edit_table <- all_edit_table[sapply(strsplit(all_edit_table, split = "_"),"[[",2) == "wk10"]

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    if (dim(edit_table_tmp)[1] == 0){
        print("empty file ignored")
    }
    else {
    edit_table_tmp$Age <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 3)
    edit_table_tmp$Section <- sapply(strsplit(i, split = "_"), "[[", 4)
    edit_table_tmp$CellType <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
    }
}

length(edit_table) # 83
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sn <- edit_table
#edit_table_sn <- do.call(rbind, edit_table_sn)
#dim(edit_table_sn)
#write.table(edit_table_sn,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sn_noACALAP <- edit_table_sn[!grepl("^AC|^AL|^AP", edit_table_sn$Gene),]
#write.table(edit_table_sn_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

# drop tables that are empty after filtering
edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

edit_table_sn <- edit_table

## for each of the cell type, for each of the three repeats/for each of the genic regions, locate the same sites that are showing up in the sc and sn, scatter plot
region_order <-  c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic")
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
                    "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle")
# if x = 10, then have to drop Chondrogenic-MSC-Fib
#celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
#                    "Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle")
# if x = 20, then have to drop MyoFib
#celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Adipogenic-MSC-Fib",
#                    "Cycling-MSC-Fib","MSC-Fib","Muscle","Pericytes","Smooth-Muscle")
alu_order <- c("ALU","REP","NONREP")

# intersect(names(table(sapply(strsplit(names(edit_table_sc), split = "_"), "[[", 1))), names(table(sapply(strsplit(names(edit_table_sn), split = "_"), "[[", 1))))

# run 1 version NOT splitting the sites into repeats or genic locations first
data_to_plot <- list()
for (celltype_tmp in celltype_order) {
    print(celltype_tmp)
    sc_table_tmp <- edit_table_sc[sapply(strsplit(names(edit_table_sc), split = "_"), "[[", 1) == celltype_tmp]
    sn_table_tmp <- edit_table_sn[sapply(strsplit(names(edit_table_sn), split = "_"), "[[", 1) == celltype_tmp]

    sc_table_tmp <- do.call(rbind, sc_table_tmp)
    sn_table_tmp <- do.call(rbind, sn_table_tmp)

    inter_sites_tmp <- intersect(sc_table_tmp$Chrom_Pos, sn_table_tmp$Chrom_Pos)
    print(length(inter_sites_tmp))

    sc_table_tmp <- sc_table_tmp[sc_table_tmp$Chrom_Pos %in% inter_sites_tmp,]
    sn_table_tmp <- sn_table_tmp[sn_table_tmp$Chrom_Pos %in% inter_sites_tmp,]

    library(dplyr)
    data_to_plot_sc <- sc_table_tmp %>%
        group_by(Chrom_Pos) %>%
        summarise(meanEditLevel = mean(EditLevel)) %>%
        distinct()
    data_to_plot_sn <- sn_table_tmp %>%
        group_by(Chrom_Pos) %>%
        summarise(meanEditLevel = mean(EditLevel)) %>%
        distinct()

    data_to_plot_tmp <- data.frame(site = inter_sites_tmp, meanEdit_sc = NA, meanEdit_sn = NA, celltype = NA)
    for (i in 1:nrow(data_to_plot_tmp)){
        data_to_plot_tmp[i, "meanEdit_sc"] <- data_to_plot_sc[data_to_plot_sc$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
        data_to_plot_tmp[i, "meanEdit_sn"] <- data_to_plot_sn[data_to_plot_sn$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
    }
    data_to_plot_tmp$celltype <- celltype_tmp
    data_to_plot[[celltype_tmp]] <- data_to_plot_tmp
}
data_to_plot <- do.call(rbind, data_to_plot)

data_to_plot <- list()
for (celltype_tmp in celltype_order) {
    print(celltype_tmp)
    sc_table_tmp <- edit_table_sc[sapply(strsplit(names(edit_table_sc), split = "_"), "[[", 1) == celltype_tmp]
    sn_table_tmp <- edit_table_sn[sapply(strsplit(names(edit_table_sn), split = "_"), "[[", 1) == celltype_tmp]

    sc_table_tmp <- do.call(rbind, sc_table_tmp)
    sn_table_tmp <- do.call(rbind, sn_table_tmp)

    inter_sites_tmp <- intersect(sc_table_tmp$Chrom_Pos, sn_table_tmp$Chrom_Pos)
    print(length(inter_sites_tmp))

    sc_table_tmp <- sc_table_tmp[sc_table_tmp$Chrom_Pos %in% inter_sites_tmp,]
    sn_table_tmp <- sn_table_tmp[sn_table_tmp$Chrom_Pos %in% inter_sites_tmp,]

    library(dplyr)
    data_to_plot_sc <- sc_table_tmp %>%
        group_by(Chrom_Pos) %>%
        summarise(meanEditLevel = mean(EditLevel), n=n()) %>%
        distinct()
    data_to_plot_sn <- sn_table_tmp %>%
        group_by(Chrom_Pos) %>%
        summarise(meanEditLevel = mean(EditLevel), n=n()) %>%
        distinct()

    data_to_plot_tmp <- data.frame(site = inter_sites_tmp, meanEdit_sc = NA, meanEdit_sn = NA, n_sc = NA, n_sn = NA, celltype = NA)
    for (i in 1:nrow(data_to_plot_tmp)){
        data_to_plot_tmp[i, "meanEdit_sc"] <- data_to_plot_sc[data_to_plot_sc$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
        data_to_plot_tmp[i, "meanEdit_sn"] <- data_to_plot_sn[data_to_plot_sn$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
        data_to_plot_tmp[i, "n_sc"] <- data_to_plot_sc[data_to_plot_sc$Chrom_Pos == data_to_plot_tmp[i, "site"], "n"]
        data_to_plot_tmp[i, "n_sn"] <- data_to_plot_sn[data_to_plot_sn$Chrom_Pos == data_to_plot_tmp[i, "site"], "n"]

    }
    data_to_plot_tmp$celltype <- celltype_tmp
    data_to_plot[[celltype_tmp]] <- data_to_plot_tmp
}
data_to_plot <- do.call(rbind, data_to_plot)
head(data_to_plot)
dim(data_to_plot) # 21313
dim(data_to_plot[data_to_plot$n_sc >= 2 & data_to_plot$n_sn >= 2,])
data_to_plot <- data_to_plot[data_to_plot$n_sc >= 2 & data_to_plot$n_sn >= 2,]
data_to_plot <- data_to_plot[data_to_plot$n_sc >= 3 & data_to_plot$n_sn >= 3,]


#cluster_col <- c("#0000FF","#FF0000","#00C000","#F07CAB",
#                 "#FF8000","#000000","#90BFF9","#C0C0FF",
#                 "#D30B94","#A00000","#F2B77C","#00FF00","#94641F")
# if x = 10, then drop a color
#cluster_col <- c("#0000FF","#00C000","#F07CAB",
#                 "#FF8000","#000000","#90BFF9","#C0C0FF",
#                 "#D30B94","#A00000","#F2B77C","#00FF00","#94641F")
# if x = 20, then drop a color
#cluster_col <- c("#0000FF","#00C000","#F07CAB",
#                 "#FF8000","#000000","#C0C0FF",
#                 "#D30B94","#A00000","#F2B77C","#00FF00","#94641F")
# if at least 50% sample sees the site, then drop colors
cluster_col <- c("#0000FF","#00C000","#F07CAB",
                 "#FF8000","#000000","#C0C0FF",
                 "#D30B94","#A00000","#F2B77C")
names(cluster_col) <- names(table(data_to_plot$celltype))
cluster_col

cor_test_result <- cor.test(data_to_plot$meanEdit_sc, data_to_plot$meanEdit_sn)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_test_result$estimate, digits = 3), "\n ", "p = ", format(cor_test_result$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust

library(ggplot2)
# first plot all cell types on the same plot, different color
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecificediting_withouttrendline_x5y3_atleast3samples.pdf", width = 10, height = 8)
ggplot(data_to_plot, aes(x = meanEdit_sc, y = meanEdit_sn, color = celltype)) + geom_point(size=2) + scale_color_manual(values = cluster_col) + 
        #geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific Site Specific Editing (sc vs. sn)")) + #coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean Editing Level Per Site for single cell Samples") + ylab("Mean Editing Level Per Site for single nucleus Samples")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecificediting_withgreytrendline.pdf", width = 10, height = 8)
ggplot(data_to_plot, aes(x = meanEdit_sc, y = meanEdit_sn, color = celltype)) + geom_point(size=2) + scale_color_manual(values = cluster_col) + 
        geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific Site Specific Editing (sc vs. sn)")) + #coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean Editing Level Per Site for single cell Samples") + ylab("Mean Editing Level Per Site for single nucleus Samples")
dev.off()


# run 1 version NOT splitting the sites into repeats or genic locations first, remove 100% edited sites

data_to_plot <- list()
for (celltype_tmp in celltype_order) {
    print(celltype_tmp)
    sc_table_tmp <- edit_table_sc[sapply(strsplit(names(edit_table_sc), split = "_"), "[[", 1) == celltype_tmp]
    sn_table_tmp <- edit_table_sn[sapply(strsplit(names(edit_table_sn), split = "_"), "[[", 1) == celltype_tmp]

    sc_table_tmp <- do.call(rbind, sc_table_tmp)
    sn_table_tmp <- do.call(rbind, sn_table_tmp)

    sc_table_tmp <- sc_table_tmp[sc_table_tmp$EditLevel != 1,]
    sn_table_tmp <- sn_table_tmp[sn_table_tmp$EditLevel != 1,]

    inter_sites_tmp <- intersect(sc_table_tmp$Chrom_Pos, sn_table_tmp$Chrom_Pos)
    print(length(inter_sites_tmp))

    sc_table_tmp <- sc_table_tmp[sc_table_tmp$Chrom_Pos %in% inter_sites_tmp,]
    sn_table_tmp <- sn_table_tmp[sn_table_tmp$Chrom_Pos %in% inter_sites_tmp,]

    library(dplyr)
    data_to_plot_sc <- sc_table_tmp %>%
        group_by(Chrom_Pos) %>%
        summarise(meanEditLevel = mean(EditLevel)) %>%
        distinct()
    data_to_plot_sn <- sn_table_tmp %>%
        group_by(Chrom_Pos) %>%
        summarise(meanEditLevel = mean(EditLevel)) %>%
        distinct()

    data_to_plot_tmp <- data.frame(site = inter_sites_tmp, meanEdit_sc = NA, meanEdit_sn = NA, celltype = NA)
    for (i in 1:nrow(data_to_plot_tmp)){
        data_to_plot_tmp[i, "meanEdit_sc"] <- data_to_plot_sc[data_to_plot_sc$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
        data_to_plot_tmp[i, "meanEdit_sn"] <- data_to_plot_sn[data_to_plot_sn$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
    }
    data_to_plot_tmp$celltype <- celltype_tmp
    data_to_plot[[celltype_tmp]] <- data_to_plot_tmp
}
data_to_plot <- do.call(rbind, data_to_plot)

cluster_col <- c("#0000FF","#FF0000","#00C000","#F07CAB",
                 "#FF8000","#000000","#90BFF9","#C0C0FF",
                 "#D30B94","#A00000","#F2B77C","#00FF00","#94641F")
names(cluster_col) <- names(table(data_to_plot$celltype))
cluster_col

cor_test_result <- cor.test(data_to_plot$meanEdit_sc, data_to_plot$meanEdit_sn)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_test_result$estimate, digits = 3), "\n ", "p = ", format(cor_test_result$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust

library(ggplot2)
# first plot all cell types on the same plot, different color
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecificediting_no100sites_withouttrendline.pdf", width = 10, height = 8)
ggplot(data_to_plot, aes(x = meanEdit_sc, y = meanEdit_sn, color = celltype)) + geom_point(size=2) + scale_color_manual(values = cluster_col) + 
        #geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific Site Specific Editing (sc vs. sn)")) + #coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean Editing Level Per Site for single cell Samples") + ylab("Mean Editing Level Per Site for single nucleus Samples")
dev.off()

pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecificediting_no100sites_withgreytrendline.pdf", width = 10, height = 8)
ggplot(data_to_plot, aes(x = meanEdit_sc, y = meanEdit_sn, color = celltype)) + geom_point(size=2) + scale_color_manual(values = cluster_col) + 
        geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific Site Specific Editing (sc vs. sn)")) + #coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean Editing Level Per Site for single cell Samples") + ylab("Mean Editing Level Per Site for single nucleus Samples")
dev.off()





# second plot, let's break each cell type down, but not separating out the Alu's and the genic locations
p <- list()
for (celltype_tmp in celltype_order) {
    data_to_plot_tmp <- data_to_plot[data_to_plot$celltype == celltype_tmp,]

    print(dim(data_to_plot_tmp))

    if (dim(data_to_plot_tmp)[1] > 3){
            cor_test_result <- cor.test(data_to_plot_tmp$meanEdit_sc, data_to_plot_tmp$meanEdit_sn)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_test_result$estimate, digits = 3), "\n ", "p = ", format(cor_test_result$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust

    p[[celltype_tmp]] <- ggplot(data_to_plot_tmp, aes(x = meanEdit_sc, y = meanEdit_sn, color = celltype)) + geom_point(size=2) + scale_color_manual(values = cluster_col) + 
        #geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific Site Specific Editing (sc vs. sn) in ", celltype_tmp, " (n = ", dim(data_to_plot_tmp)[1], ")")) + #coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16), legend.position = "none")  + xlab("Mean Editing Level Per Site for single cell Samples") + ylab("Mean Editing Level Per Site for single nucleus Samples")
    }
}
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecificediting_foreachcelltype_withouttrendline.pdf", width = 8, height = 8)
p
dev.off()


# run 2 version splitting the sites into genic locations
p <- list()
for (celltype_tmp in celltype_order) {
    print(celltype_tmp)
    for (region_tmp in region_order){
        sc_table_tmp <- edit_table_sc[sapply(strsplit(names(edit_table_sc), split = "_"), "[[", 1) == celltype_tmp]
        sn_table_tmp <- edit_table_sn[sapply(strsplit(names(edit_table_sn), split = "_"), "[[", 1) == celltype_tmp]

        sc_table_tmp <- do.call(rbind, sc_table_tmp)
        sn_table_tmp <- do.call(rbind, sn_table_tmp)

        sc_table_tmp <- sc_table_tmp[sc_table_tmp$Region == region_tmp, ]
        sn_table_tmp <- sn_table_tmp[sn_table_tmp$Region == region_tmp, ]

        inter_sites_tmp <- intersect(sc_table_tmp$Chrom_Pos, sn_table_tmp$Chrom_Pos)
        print(region_tmp)
        print(length(inter_sites_tmp))

        if (length(inter_sites_tmp) > 3){
            sc_table_tmp <- sc_table_tmp[sc_table_tmp$Chrom_Pos %in% inter_sites_tmp,]
            sn_table_tmp <- sn_table_tmp[sn_table_tmp$Chrom_Pos %in% inter_sites_tmp,]

            library(dplyr)
            data_to_plot_sc <- sc_table_tmp %>%
                group_by(Chrom_Pos) %>%
                summarise(meanEditLevel = mean(EditLevel)) %>%
                distinct()
            data_to_plot_sn <- sn_table_tmp %>%
                group_by(Chrom_Pos) %>%
                summarise(meanEditLevel = mean(EditLevel)) %>%
                distinct()

            data_to_plot_tmp <- data.frame(site = inter_sites_tmp, meanEdit_sc = NA, meanEdit_sn = NA, celltype = NA)
            for (i in 1:nrow(data_to_plot_tmp)){
                data_to_plot_tmp[i, "meanEdit_sc"] <- data_to_plot_sc[data_to_plot_sc$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
                data_to_plot_tmp[i, "meanEdit_sn"] <- data_to_plot_sn[data_to_plot_sn$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
            }
            data_to_plot_tmp$celltype <- celltype_tmp

            cor_test_result <- cor.test(data_to_plot_tmp$meanEdit_sc, data_to_plot_tmp$meanEdit_sn)
            annotations <- data.frame(
                    xpos = c(-Inf),
                    ypos =  c(Inf),
                    annotateText = c(paste0(" R = ", format(cor_test_result$estimate, digits = 3), "\n ", "p = ", format(cor_test_result$p.value, digits = 3))),
                    hjustvar = c(0) ,
                    vjustvar = c(1)) #<- adjust

    p[[paste0(celltype_tmp, "_", region_tmp)]] <- ggplot(data_to_plot_tmp, aes(x = meanEdit_sc, y = meanEdit_sn, color = celltype)) + geom_point(size=2) + scale_color_manual(values = cluster_col) + 
        #geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific Site Specific Editing (sc vs. sn) in ", celltype_tmp, " in ",region_tmp, " region (n = ", dim(data_to_plot_tmp)[1], ")")) + #coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16), legend.position = "none")  + xlab("Mean Editing Level Per Site for single cell Samples") + ylab("Mean Editing Level Per Site for single nucleus Samples")

        }
    }
}
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecificediting_foreachcelltype_foreachgenicregion_withouttrendline.pdf", width = 10, height = 10)
p
dev.off()

# run 2 version splitting the sites into alu repeat types
p <- list()
for (celltype_tmp in celltype_order) {
    print(celltype_tmp)
    for (alu_tmp in alu_order){
        sc_table_tmp <- edit_table_sc[sapply(strsplit(names(edit_table_sc), split = "_"), "[[", 1) == celltype_tmp]
        sn_table_tmp <- edit_table_sn[sapply(strsplit(names(edit_table_sn), split = "_"), "[[", 1) == celltype_tmp]

        sc_table_tmp <- do.call(rbind, sc_table_tmp)
        sn_table_tmp <- do.call(rbind, sn_table_tmp)

        sc_table_tmp <- sc_table_tmp[sc_table_tmp$Alu == alu_tmp, ]
        sn_table_tmp <- sn_table_tmp[sn_table_tmp$Alu == alu_tmp, ]

        inter_sites_tmp <- intersect(sc_table_tmp$Chrom_Pos, sn_table_tmp$Chrom_Pos)
        print(alu_tmp)
        print(length(inter_sites_tmp))

        if (length(inter_sites_tmp) > 3){
            sc_table_tmp <- sc_table_tmp[sc_table_tmp$Chrom_Pos %in% inter_sites_tmp,]
            sn_table_tmp <- sn_table_tmp[sn_table_tmp$Chrom_Pos %in% inter_sites_tmp,]

            library(dplyr)
            data_to_plot_sc <- sc_table_tmp %>%
                group_by(Chrom_Pos) %>%
                summarise(meanEditLevel = mean(EditLevel)) %>%
                distinct()
            data_to_plot_sn <- sn_table_tmp %>%
                group_by(Chrom_Pos) %>%
                summarise(meanEditLevel = mean(EditLevel)) %>%
                distinct()

            data_to_plot_tmp <- data.frame(site = inter_sites_tmp, meanEdit_sc = NA, meanEdit_sn = NA, celltype = NA)
            for (i in 1:nrow(data_to_plot_tmp)){
                data_to_plot_tmp[i, "meanEdit_sc"] <- data_to_plot_sc[data_to_plot_sc$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
                data_to_plot_tmp[i, "meanEdit_sn"] <- data_to_plot_sn[data_to_plot_sn$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
            }
            data_to_plot_tmp$celltype <- celltype_tmp

            cor_test_result <- cor.test(data_to_plot_tmp$meanEdit_sc, data_to_plot_tmp$meanEdit_sn)
            annotations <- data.frame(
                    xpos = c(-Inf),
                    ypos =  c(Inf),
                    annotateText = c(paste0(" R = ", format(cor_test_result$estimate, digits = 3), "\n ", "p = ", format(cor_test_result$p.value, digits = 3))),
                    hjustvar = c(0) ,
                    vjustvar = c(1)) #<- adjust

    p[[paste0(celltype_tmp, "_", alu_tmp)]] <- ggplot(data_to_plot_tmp, aes(x = meanEdit_sc, y = meanEdit_sn, color = celltype)) + geom_point(size=2) + scale_color_manual(values = cluster_col) + 
        #geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Cell Type Specific Site Specific Editing (sc vs. sn) in ", celltype_tmp, " in ",alu_tmp, " region (n = ", dim(data_to_plot_tmp)[1], ")")) + #coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16), legend.position = "none")  + xlab("Mean Editing Level Per Site for single cell Samples") + ylab("Mean Editing Level Per Site for single nucleus Samples")

        }
    }
}
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecificediting_foreachcelltype_foreachalutype_withouttrendline.pdf", width = 10, height = 10)
p
dev.off()

###### site specific version series 2 by each sample
## reference table
# dbsnp filter
REDIref <- read.table("/media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_GRCh38.txt", header = T, sep = "\t")
REDIref$Chrom_Pos <- paste0(REDIref$Region,"_",REDIref$Position)
head(REDIref)
dim(REDIref) # [1] 15681871       27
dim(REDIref[REDIref$dbsnp == "-",]) # [1] 13569494       27
SNP_sites <- REDIref[REDIref$dbsnp != "-","Chrom_Pos"]; head(SNP_sites)
ACALAP_sites <- REDIref[grepl("^AC|^AL|^AP", REDIref$Gene.wgEncodeGencodeBasicV34),"Chrom_Pos"]; head(ACALAP_sites); length(ACALAP_sites)
#length(REDIref[grep("AC[[:digit:]]|AL[[:digit:]]|AP[[:digit:]]", REDIref$Gene.wgEncodeGencodeBasicV34),"Chrom_Pos"])
MITO_sites <- REDIref[REDIref$Region == "GRCh38_chrM","Chrom_Pos"]; head(MITO_sites)


## read in sc site specific editing
input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_H1_v7/"
edit_table <- list()

all_edit_table <- list.files(input_dir)

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
}

length(edit_table) # 59
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sc <- edit_table
#edit_table_sc <- do.call(rbind, edit_table_sc)
#dim(edit_table_sc)
#write.table(edit_table_sc,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sc_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sc_noACALAP <- edit_table_sc[!grepl("^AC|^AL|^AP", edit_table_sc$Gene),]
#write.table(edit_table_sc_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sc_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

edit_table_sc <- edit_table

## read in sn site specific editing, keep wk10 only

input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/edit_output_TS_v4/"
edit_table <- list()

all_edit_table <- list.files(input_dir)
all_edit_table <- all_edit_table[sapply(strsplit(all_edit_table, split = "_"),"[[",1) == "wk10"]

for (i in all_edit_table){
    edit_table_tmp <- read.table(paste0(input_dir,i), header = 1)
    if (dim(edit_table_tmp)[1] == 0){
        print("empty file ignored")
    }
    else {
    edit_table_tmp$Age <- sapply(strsplit(i, split = "_"), "[[", 1)
    edit_table_tmp$Teratoma <- sapply(strsplit(i, split = "_"), "[[", 2)
    edit_table_tmp$Section <- sapply(strsplit(i, split = "_"), "[[", 3)
    edit_table_tmp$Chrom_Pos <- paste0(edit_table_tmp$Chrom, "_", edit_table_tmp$Position)
          
    edit_table[[sapply(strsplit(i, split = "_dedup.txt"), "[[", 1)]] <- edit_table_tmp
    }
}

length(edit_table) # 83
head(edit_table[[1]])

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$NonRefEdited >= 3,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$TotalEditedReads >= 5,]
    edit_table[[i]] <- edit_table[[i]][edit_table[[i]]$Chrom != "GRCh38_chrM",]
    edit_table[[i]] <- edit_table[[i]][!edit_table[[i]]$Chrom_Pos %in% SNP_sites,]
}
#edit_table_sn <- edit_table
#edit_table_sn <- do.call(rbind, edit_table_sn)
#dim(edit_table_sn)
#write.table(edit_table_sn,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP.txt", sep = ",", col.names = T, row.names = T, quote = F)

for (i in names(edit_table)){
    edit_table[[i]] <- edit_table[[i]][!grepl("^AC|^AL|^AP", edit_table[[i]]$Gene),]
}

#edit_table_sn_noACALAP <- edit_table_sn[!grepl("^AC|^AL|^AP", edit_table_sn$Gene),]
#write.table(edit_table_sn_noACALAP,"Figure1_WildtypeH1_R1/20240307_scVSsn/edit_table_sn_minEdit3_totEdit5_noM_noSNP_noACALAP.txt", sep = ",", col.names = T, row.names = T, quote = F)
sapply(edit_table, dim)[1,]

# drop tables that are empty after filtering
edit_table <- edit_table[sapply(edit_table,dim)[1,] != 0]

edit_table_sn <- edit_table

## for each of the cell type, for each of the three repeats/for each of the genic regions, locate the same sites that are showing up in the sc and sn, scatter plot
region_order <-  c("exonic","splicing","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","UTR5","UTR3","intronic","upstream","downstream","intergenic")
alu_order <- c("ALU","REP","NONREP")

# intersect(names(table(sapply(strsplit(names(edit_table_sc), split = "_"), "[[", 1))), names(table(sapply(strsplit(names(edit_table_sn), split = "_"), "[[", 1))))

# run 1 version NOT splitting the sites into repeats or genic locations first

    sc_table_tmp <- edit_table_sc
    sn_table_tmp <- edit_table_sn

    sc_table_tmp <- do.call(rbind, sc_table_tmp)
    sn_table_tmp <- do.call(rbind, sn_table_tmp)

    inter_sites_tmp <- intersect(sc_table_tmp$Chrom_Pos, sn_table_tmp$Chrom_Pos)
    print(length(inter_sites_tmp))

    sc_table_tmp <- sc_table_tmp[sc_table_tmp$Chrom_Pos %in% inter_sites_tmp,]
    sn_table_tmp <- sn_table_tmp[sn_table_tmp$Chrom_Pos %in% inter_sites_tmp,]

    library(dplyr)
    data_to_plot_sc <- sc_table_tmp %>%
        group_by(Chrom_Pos) %>%
        summarise(meanEditLevel = mean(EditLevel), n=n()) %>%
        distinct()
    data_to_plot_sn <- sn_table_tmp %>%
        group_by(Chrom_Pos) %>%
        summarise(meanEditLevel = mean(EditLevel), n=n()) %>%
        distinct()

    data_to_plot_tmp <- data.frame(site = inter_sites_tmp, meanEdit_sc = NA, meanEdit_sn = NA, n_sc = NA, n_sn = NA)
    for (i in 1:nrow(data_to_plot_tmp)){
        data_to_plot_tmp[i, "meanEdit_sc"] <- data_to_plot_sc[data_to_plot_sc$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
        data_to_plot_tmp[i, "meanEdit_sn"] <- data_to_plot_sn[data_to_plot_sn$Chrom_Pos == data_to_plot_tmp[i, "site"], "meanEditLevel"]
        data_to_plot_tmp[i, "n_sc"] <- data_to_plot_sc[data_to_plot_sc$Chrom_Pos == data_to_plot_tmp[i, "site"], "n"]
        data_to_plot_tmp[i, "n_sn"] <- data_to_plot_sn[data_to_plot_sn$Chrom_Pos == data_to_plot_tmp[i, "site"], "n"]
    }

data_to_plot <- data_to_plot_tmp
table(data_to_plot$n_sc) ; table(data_to_plot$n_sn)

head(data_to_plot)
dim(data_to_plot) # 21313
dim(data_to_plot[data_to_plot$n_sc >= 2 & data_to_plot$n_sn >= 2,])
data_to_plot <- data_to_plot[data_to_plot$n_sc >= 2 & data_to_plot$n_sn >= 2,] ; table(data_to_plot$n_sc) ; table(data_to_plot$n_sn)
data_to_plot <- data_to_plot[data_to_plot$n_sc >= 3 & data_to_plot$n_sn >= 3,] ; table(data_to_plot$n_sc) ; table(data_to_plot$n_sn)
data_to_plot <- data_to_plot[data_to_plot$n_sc >= 4 & data_to_plot$n_sn >= 4,] ; table(data_to_plot$n_sc) ; table(data_to_plot$n_sn)
data_to_plot <- data_to_plot[data_to_plot$n_sc >= 4 & data_to_plot$n_sn >= 5,] ; table(data_to_plot$n_sc) ; table(data_to_plot$n_sn)


cor_test_result <- cor.test(data_to_plot$meanEdit_sc, data_to_plot$meanEdit_sn)
annotations <- data.frame(
        xpos = c(-Inf),
        ypos =  c(Inf),
        annotateText = c(paste0(" R = ", format(cor_test_result$estimate, digits = 3), "\n ", "p = ", format(cor_test_result$p.value, digits = 3))),
        hjustvar = c(0) ,
        vjustvar = c(1)) #<- adjust

library(ggplot2)
# first plot all cell types on the same plot, different color
pdf("Figure1_WildtypeH1_R1/20240307_scVSsn/scVSsn_sitespecificediting_withouttrendline_teratomalevel_x5y3_atleast3samples_lightblue.pdf", width = 10, height = 10)
ggplot(data_to_plot, aes(x = meanEdit_sc, y = meanEdit_sn)) + geom_point(size=2, color = "lightblue") +
        #geom_smooth(linewidth=1.5, se = F, method = "lm", color = "grey", fullrange = T) + 
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = "black") + 
    theme_classic() + ggtitle(paste0("Teratoma Level Site Specific Editing (sc vs. sn)")) + #coord_cartesian(xlim = c(0.45,0.6), ylim = c(1, 1.7)) + 
    theme(axis.text=element_text(size=20,angle = 90, vjust = 0.5, hjust=1),axis.title=element_text(size=18,face="bold"), axis.title.y.right = element_text(size=18,face="bold"),
    legend.title = element_text(size=18), #change legend title font size
            legend.text = element_text(size=16))  + xlab("Mean Editing Level Per Site for single cell Samples") + ylab("Mean Editing Level Per Site for single nucleus Samples")
dev.off()

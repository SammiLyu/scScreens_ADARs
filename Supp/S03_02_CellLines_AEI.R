library(Seurat)

load("Supp1_CellLine_R1/WT3CL_YWsubset_intObj_20230418.Robj")

### obtain AEI through RNAEditingIndexer

input_dir <- "/media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_CL_v1/output_index_20230420/"
aei_table <- list()
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","HSC","Adipogenic-MSC-Fib",
                    "Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle",
                    "Kidney-Prog")

for (i in names(table(WT_3CL$annot_v1))){
  for (j in names(table(WT_3CL$teratoma))){
    file_tmp <- paste0(input_dir, i, "_", j, "_dedup/EditingIndex.csv")
    if (file.exists(file_tmp)) {
      print("file is there")
      aei_table_tmp <- read.table(file_tmp, header = 1, sep = ",")
      aei_table[[paste0(i, "_", j)]] <- aei_table_tmp
    }
  }
}

aei_table

# reformat
aei_table_df <- data.frame(matrix(NA, nrow = length(aei_table), ncol = ncol(aei_table[[1]])))
row.names(aei_table_df) <- names(aei_table)
colnames(aei_table_df) <- colnames(aei_table[[1]])
for (i in names(aei_table)){
  aei_table_df[i,] <- aei_table[i][[1]]
}

barplot(aei_table_df$A2GEditingIndex)

write.table(aei_table_df, 'Supp1_CellLine_R1/AEI_CL_bySample_byTer_20230424.txt', col.names = T, row.names = T, quote = FALSE, sep = ',')
aei_table_df <- read.table('Supp1_CellLine_R1/AEI_CL_bySample_byTer_20230424.txt', header = T, sep =',')


### obtain TPM from featureCounts
adar_tpm <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/fc_CL_v1/adar_exp_tpm.txt', header = T)
adar_tpm$Sample <- sapply(strsplit(adar_tpm$Sample, split = '20230419_'), "[[", 2)
adar_tpm
adar_tpm$CellType <- sapply(strsplit(adar_tpm$Sample, split = '_PG|_HU|_H9'), "[[", 1)
adar_tpm$CellType <- str_replace_all(string = adar_tpm$CellType, pattern = "_", replacement = "-")
adar_tpm$CL <- gsub("^.*_", "", adar_tpm$Sample)
write.table(adar_tpm, 'Supp1_CellLine_R1/adarTPM_bySample_byTer_20230507.txt', col.names = T, row.names = F, quote = FALSE, sep = ',')
adar_tpm <- read.table('Supp1_CellLine_R1/adarTPM_bySample_byTer_20230507.txt', header = T, sep =',')

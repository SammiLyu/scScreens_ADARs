library(WGCNA)
library(Seurat)


load(file = "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")
KO_sub

Idents(KO_sub) <- KO_sub$annot_v1

KO_sub <- ScaleData(KO_sub)

for (i in names(table(KO_sub$annot_v1))){
  celltype=i
  print(i)
  seur<-subset(KO_sub, idents = i)
  
  myPower=1
  
  savename=paste("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/WGCNA_Results/Genes_",celltype,".Robj",sep="")
  
  print("Run WGCNA!")
  
  net<-RunWGCNA(seur,savename,myPower)
  save(net,file=savename)	
}

net_list <- list()
for (i in names(table(KO_sub$annot_v1))){
  load(paste0("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/WGCNA_Results/Genes_",i,".Robj"))
  net_list[[i]] <- net
}
sapply(net_list, length)
sum(sapply(net_list, length))

Gene_tab <- data.frame(matrix(NA, nrow = sum(sapply(net_list, length)), ncol = 2, dimnames = list(NULL, c("CellType","Gene"))))
Gene_tab$CellType <- gsub("[0-9]*$","",names(unlist(net_list)))
Gene_tab$Gene <- unlist(net_list)
dim(Gene_tab)

write.table(Gene_tab,'/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/WGCNA_module_genes_20230308.csv', quote = F, sep = ",",row.names = F)

Reduce(intersect, net_list)
Reduce(rbind, net_list)

heat.mat <- t(apply(GetAssayData(KO_sub, assay = "RNA", slot = "scale.data")[Gene_tab$Gene,], 1, function(x) {
  tapply(x, KO_sub$annot_v1, mean)
}))
heat.mat[heat.mat > 5] <- 5

ggHeat(heat.mat, clustering = "none", x.lab.size = 11, y.lab.size = 11)

##
##Basic code used for WGCNA 
##
RunWGCNA<-function(seur,output,myPower)
{
  print("Filter Data!")
  dat=as.matrix(seur@assays$RNA@scale.data)
  dat<-dat[!(rowMeans(dat, na.rm = TRUE) < 0.25),] #get rid of lowly expressed genes
  dat<-dat[!(rowMeans(dat, na.rm = TRUE) > 9),] # get rid of highly expressed genes
  
  dat <- dat[!(rownames(dat) %in% grep(pattern = "^RP", x = rownames(dat), value = TRUE)),] ##Remove ribosomal genes
  
  dat=t(dat)
  
  #save(colnames(dat),file=paste(output,".genes.Robj",sep=""))
  
  return(colnames(dat))
  
  print("Remaining data:")
  print(head(colnames(dat)))
  print(head(rownames(dat)))
  
  print(dim(dat))
  
  print("Get modules!")
  
  
  networkType = "signed"
  minModuleSize = 7 
  mergeCutHeight = 0.15
  net = blockwiseModules(
    # Input data
    dat,
    
    # Data checking options
    checkMissingData = TRUE,
    
    # Options for splitting data into blocks
    blocks = NULL,
    maxBlockSize = 5000,   # xin jin change this from 5000 to 2000
    randomSeed = 59069,
    
    # load TOM from previously saved file?
    loadTOM = FALSE,
    
    # Network construction arguments: correlation options
    corType = "bicor", # more robust for non-normal? Song 2012
    maxPOutliers = 0.1,    # XinJin changed from dedault 0.9
    quickCor = 0,
    pearsonFallback = "individual",
    cosineCorrelation = FALSE,
    
    # Adjacency function options
    # this is where power gets input
    power = myPower,
    networkType = networkType,
    
    # Topological overlap options
    TOMType = "signed",
    TOMDenom = "min",
    
    # Saving or returning TOM
    getTOMs = NULL,
    saveTOMs = TRUE,
    saveTOMFileBase = "blockwiseTOM",
    
    # Basic tree cut options
    deepSplit = 3,   # xin jin changed from 2 to 3 to break down module size
    detectCutHeight = 0.995,
    minModuleSize = minModuleSize,
    
    # Advanced tree cut options
    maxCoreScatter = NULL, minGap = NULL,
    maxAbsCoreScatter = NULL, minAbsGap = NULL,
    minSplitHeight = NULL, minAbsSplitHeight = NULL,
    useBranchEigennodeDissim = FALSE,
    minBranchEigennodeDissim = mergeCutHeight,
    pamStage = TRUE, pamRespectsDendro = TRUE,
    
    # Gene reassignment, module trimming, and module "significance" criteria
    reassignThreshold = 1e-6,
    minCoreKME = 0.5,
    minCoreKMESize = minModuleSize/3,
    minKMEtoStay = 0.3,
    
    # Module merging options
    mergeCutHeight = 0.15,
    impute = TRUE,
    trapErrors = TRUE,
    
    # Output options
    numericLabels = FALSE,
    
    # Options controlling behaviour
    nThreads = 0,
    verbose = 3, 
    indent = 0)
  
  mods=NULL
  
  print("Save modules!")
  save(net,file=output)
  
  print("Done!")
  
  return(net)
  
}


if(!interactive())
{
  print("Load data and get arguments")
  load("../SeuratObjects/key.500genes.Robj")
  args = commandArgs(trailingOnly=TRUE)
  
  celltype=args[1]
  seur<-SubsetData(key,WhichCells(key,celltype))
  
  myPower=as.numeric(args[2])
  
  savename=paste("Results/Genes.",celltype,".Robj",sep="")
  
  print("Run WGCNA!")
  
  net<-RunWGCNA(seur,savename,myPower)
  save(net,file=savename)	
  print("Done!")
}
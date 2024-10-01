### https://github.com/klarman-cell-observatory/ivPerturbSeq/blob/master/CellComposition/AdamsMethod.R
### only modified to remove cell types of small cell number


library(Seurat)
library(dplyr)
library(tidyr)



### redoing cell type composition collapsing germ layers
load("/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")
KO_sub

elltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Hematopoietic",
                    "Adipogenic-MSC-Fib","Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle",
                    "Kidney-Prog")

celltype_order <- c("Smooth-Muscle","SCP","Retinal-Epi","Pericytes","Neurons","Neural-Progenitors","MyoFib","Muscle","MSC-Fib","Kidney-Prog",
                    "Hematopoietic","Gut-Epi","Cycling-MSC-Fib","Chondrogenic-MSC-Fib",
                    "Adipogenic-MSC-Fib")

tab <- CellComp_Poisson(KO_sub, celltype = "annot_v1",perturatbations = "guide",batch = "exp",ref = "NTC", cutoff=100, min_cells_celltype = 0)
tab$CellType <- factor(tab$CellType, levels = celltype_order)

pdf('Figure3_ADARKO_R1/KO_4ter_ivPerturb_celltypeComp_norm_NTC_minCells0_20230309.pdf')
ggplot(tab,aes(x=Gene,y=CellType,fill=Estimate,size=-log(padj,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

write.table(tab, "Figure3_ADARKO_R1/KO_4ter_ivPerturb_celltypeComp_norm_NTC_minCells0_20230309.csv", quote = F, row.names = F, sep = ",")

#tab <- CellComp_Poisson(KO_sub, celltype = "annot_v1",perturatbations = "guide",batch = "exp",ref = "NTC", cutoff=100, min_cells_celltype = 10)
#tab$CellType <- factor(tab$CellType, levels = celltype_order)

#pdf('Figure3_ADARKO_R1/KO_4ter_ivPerturb_celltypeComp_norm_NTC_minCells10_20230309.pdf')
#ggplot(tab,aes(x=Gene,y=CellType,fill=Estimate,size=-log(padj,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#dev.off()

#tab <- tab[tab$padj < 0.05,]
#pdf('Figure3_ADARKO_R1/KO_4ter_ivPerturb_celltypeComp_norm_NTC_sig_20230309.pdf')
#ggplot(tab,aes(x=Gene,y=CellType,fill=Estimate,size=-log(padj,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#dev.off()

# norm to AAVS1

tab <- CellComp_Poisson(KO_sub, celltype = "annot_v1",perturatbations = "guide",batch = "exp",ref = "AAVS1", cutoff=100, min_cells_celltype = 0)
tab$CellType <- factor(tab$CellType, levels = celltype_order)

pdf('Figure3_ADARKO_R1/KO_4ter_ivPerturb_celltypeComp_norm_AAVS1_minCells0_20230309.pdf')
ggplot(tab,aes(x=Gene,y=CellType,fill=Estimate,size=-log(padj,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

write.table(tab, "Figure3_ADARKO_R1/KO_4ter_ivPerturb_celltypeComp_norm_AAVS1_minCells0_20230309.csv", quote = F, row.names = F, sep = ",")


tab <- tab[tab$padj < 0.05,]
ggplot(tab,aes(x=Gene,y=CellType,fill=Estimate,size=-log(pval,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##Tests to see if differences in composition between perturbations. Fee in:
##seur, the Seurat object
##celltype, the column in meta.data containing celltype info
##perturbation, the column in meta.data containing celltype info
##batch, the column in meta.data containing batch info
##cutoff, the min number of cells per batch per guide to perform test
##min_cells_celltype, the min number of cells needed in AAVS1-KO cells for a cell type to be included
##
CellComp_Poisson<-function(seur,celltype="CellType",perturatbations="perturbation",batch="batch",ref = "ref", cutoff=50, min_cells_celltype = 10)
{
  
  print("Clean Data")
  meta=seur@meta.data[,c(celltype,perturatbations,batch)]
  colnames(meta)=c("CellType","Pert","Batch")
  
  meta <- meta %>% group_by(CellType,Pert,Batch) %>% summarise(Num=length(CellType)) %>% as.data.frame()
  
  meta2<-meta %>% group_by(Pert,Batch) %>% summarise(Tot=sum(Num)) %>% as.data.frame()
  meta=left_join(meta,meta2)
  
  meta3 <- meta %>% group_by(Pert, CellType) %>% summarise(Tot=sum(Num)) %>% as.data.frame() 
  cts_to_keep <- meta3[meta3$Pert == "AAVS1" & meta3$Tot > min_cells_celltype, "CellType"]
  
  ##so meta is dataframe of 5 columns: celltype, perturbation, batch, Number of total cells of that celltype/pert/batch, and number of total cells of that pert/batch
  
  meta=meta[meta[,"Tot"]>cutoff,]
  meta=meta[meta$CellType %in% cts_to_keep,]
  
  meta["Pert"]=relevel(factor(meta[,"Pert"]),ref=ref)
  
  lst=list()
  for(i in unique(meta[,"CellType"])){lst[[i]]=meta[meta[,"CellType"]==i,]}
  
  
  print("Fit model!")
  out<-lapply(lst,function(cur){
    celltype=cur[1,"CellType"]
    
    print(celltype)
    cur["logTot"]=log(cur[,"Tot"])
    fit<-glm(Num~offset(logTot)+Batch+Pert,data=cur,family="poisson")
    tab=summary(fit)
    tab=tab$coefficients
    tab=data.frame(tab)
    tab=tab[grep("Pert",rownames(tab)),]
    tab["Gene"]=sub("Pert","",rownames(tab))
    tab["CellType"]=celltype
    tab=tab[,c(5,6,4,1,2,3)]
    colnames(tab)[3]="pval"
    return(tab)
  })
  
  tab=do.call(rbind,out)
  
  rownames(tab)=NULL
  
  tab=tab[order(tab[,"pval"]),]
  
  tab["padj"]=p.adjust(tab[,"pval"],"fdr")
  
  print("Done!")
  
  return(tab)
  
}




ErrorBars_Poisson<-function(seur,celltype="CellType",perturatbations="perturbation",batch="batch",ref = "ref",cutoff=10)
{
  
  print("Clean Data")
  meta=seur@meta.data[,c(celltype,perturatbations,batch)]
  colnames(meta)=c("CellType","Pert","Batch")
  
  meta <- meta %>% group_by(CellType,Pert,Batch) %>% summarise(Num=length(CellType)) %>% as.data.frame()
  
  meta2<-meta %>% group_by(Pert,Batch) %>% summarise(Tot=sum(Num)) %>% as.data.frame()
  
  meta=left_join(meta,meta2)
  
  ##so meta is dataframe of 5 columns: celltype, perturbation, batch, Number of total cells of that celltype/pert/batch, and number of total cells of that pert/batch
  
  meta=meta[meta[,"Tot"]>cutoff,]
  
  
  meta["Pert"]=relevel(factor(meta[,"Pert"]),ref=ref)
  
  lst=list()
  for(i in unique(meta[,"CellType"])){for(j in unique(meta[,"Pert"])){lst[[paste(i,j,sep="_")]]=meta[meta[,"CellType"]==i & meta[,"Pert"]==j,]}}
  #ieta["Pert"]=relevel(factor(meta[,"Pert"]),ref=ref)
  
  
  print("Fit model!")
  out<-lapply(lst,function(cur){
    celltype=cur[1,"CellType"]
    
    #print(celltype)
    cur["logTot"]=log(cur[,"Tot"])
    fit<-glm(Num~offset(logTot)+1,data=cur,family="poisson")
    tab=summary(fit)
    tab=tab$coefficients
    tab=data.frame(tab)
    tab["Num"]=dim(cur)[1]
    #tab=tab[grep("Pert",rownames(tab)),]
    #tab["Gene"]=sub("Pert","",rownames(tab))
    #tab["CellType"]=celltype
    #tab=tab[,c(5,6,4,1,2,3)]
    #colnames(tab)[3]="pval"
    return(tab)
  })
  
  tab=do.call(rbind,out)
  
  tab=data.frame(tab)
  
  tab["Name"]=rownames(tab)
  tab["CellType"]=as.character(lapply(rownames(tab),function(x){strsplit(x,"_")[[1]][1]}))
  tab["Pert"]=as.character(lapply(rownames(tab),function(x){strsplit(x,"_")[[1]][2]}))
  
  return(tab)
  
  
  rownames(tab)=NULL
  
  tab=tab[order(tab[,"pval"]),]
  
  tab["padj"]=p.adjust(tab[,"pval"],"fdr")
  
  print("Done!")
  
  return(tab)
  
}

if(!interactive())
{
  print("Load Seurat")
  load("../SeuratObjects/key.500genes.Robj")
  tab=CellComp_Poisson(key,celltype="CellType",perturatbations="perturbation",batch="batch",cutoff=10)
  save(tab,file="CellComposition.Poisson.Robj")
  write.table(tab,file="CellComposition.Poisson.txt",sep="\t",quote=F,row.names=F)
  p=ggplot(tab,aes(x=Gene,y=CellType,fill=Estimate,size=-log(pval,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave("CellComposition.Poisson.pdf",p,width=14,height=5)
  
  
  
  #tab<-ErrorBars_Poisson(key)
  #save(tab,file="ForError.Bars.Robj")
  
}
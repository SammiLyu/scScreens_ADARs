source("/media/Scratch_SSD_Voyager/sammi/RNA_editing/scripts/PertLM.R")
source("/media/Scratch_SSD_Voyager/sammi/RNA_editing/scripts/WGCNA.Get.Eigen.R")
source("/media/Scratch_SSD_Voyager/sammi/RNA_editing/scripts/BatchSample.v2.R")

library(Seurat)

load(file = "/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_ADARB1-1rm_20230207.rda")

KO_sub$batch <- KO_sub$exp
KO_sub$perturbations <- KO_sub$guide

Idents(KO_sub) <- KO_sub$annot_v1

Gene_tab <- read.table('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/WGCNA_module_genes_20230308.csv', quote = "", sep = ",", row.names = NULL, header = T)

res <- TestWGCNA(KO_sub,Gene_tab,ref = "NTC", numRep=10000)
do.call(rbind, res)
pdf('Figure3_ADARKO_R1/KO_4ter_ivPerturb_transcriptDiff_norm_NTC_20230712.pdf')
drawResults(res, ref="NTC") + theme_classic() + RotatedAxis() + coord_flip()
dev.off()

pdf('Figure3_ADARKO_R1/KO_4ter_ivPerturb_transcriptDiff_norm_NTC_20230712_sideways.pdf')
drawResults(res, ref="NTC") + theme_classic() + RotatedAxis()
dev.off()

for (i in 1:length(res)){
  res[[i]][,"Cell_Type"] <- names(res)[i]
}

save(res, file = '/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/Transcript_level_effects_by_ivPerturbSeq_20230214.rda')
load('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/Transcript_level_effects_by_ivPerturbSeq_20230214.rda')

write.table(Reduce(rbind, res), file = "Figure3_ADARKO_R1/KO_4ter_ivPerturb_transcriptDiff_norm_NTC_20230309.csv", quote = F, row.names = T, sep = ",")


Idents(KO_sub) <- KO_sub$annot_v1
res2 <- TestWGCNA(KO_sub,Gene_tab,ref = "AAVS1", numRep=10000)
pdf('Figure3_ADARKO_R1/KO_4ter_ivPerturb_transcriptDiff_norm_AAVS1_20230712.pdf')
drawResults(res2, ref="AAVS1") + theme_classic() + RotatedAxis() + coord_flip()
dev.off()

pdf('Figure3_ADARKO_R1/KO_4ter_ivPerturb_transcriptDiff_norm_AAVS1_20230712_sideways.pdf')
drawResults(res2, ref="AAVS1") + theme_classic() + RotatedAxis()
dev.off()

for (i in 1:length(res2)){
  res2[[i]][,"Cell_Type"] <- names(res2)[i]
}

save(res2, file = '/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/Transcript_level_effects_by_ivPerturbSeq_normtoAAVS1_20230308.rda')
load('/media/Scratch_SSD_Voyager/sammi/RNA_editing/Figure3_ADARKO_R1/Transcript_level_effects_by_ivPerturbSeq_normtoAAVS1_20230308.rda')

write.table(Reduce(rbind, res2), file = "Figure3_ADARKO_R1/KO_4ter_ivPerturb_transcriptDiff_norm_AAVS1_20230309.csv", quote = F, row.names = T, sep = ",")


#drop kidney, SCP

TestWGCNA<-function(seur,tab,ref=ref,form=Cluster1~perturbations+batch+nGene,maxCells=0,seed=1,numRep=1000,getScore=F,minPert=10,plotIt=F)
{
  
  print("Pert only") # no need to run this subset because all cells in the feed in dataset is already filtered for perturbation
  print(length(seur))
  #seur=SubsetData(seur,names(seur@ident)[!is.na(seur@meta.data$perturbation)])
  #print(length(seur@ident))
  
  mods<-unique(tab[,1])[-c(6,14)] # remove kidney
  
  ret=list()
  
  scores<-list()
  
  for(i in mods)
  {
    print(i)
    
    print("Subset")
    tab_cur=tab[tab[,1]==i,]
    
    cur<-subset(seur, idents = tab_cur[1,"CellType"])

    
    print(length(colnames(cur)))
    
    print("Make list")
    
    genes<-as.character(tab_cur[,"Gene"])
    
    lst<-list()
    lst[["Genes"]]=genes
    
    print("Get Module")
    #cur<-AddModuleScore(cur,lst)
    cur<-getEigen(cur,lst[[1]])
    
    if(plotIt)
    {
      FeaturePlot(cur,"Cluster1")
    }
    
    if(maxCells>0)
    {
      set.seed(seed)
      print("Downsample!")
      lst<-table(cur@meta.data$guide)
      maxCells=2*median(as.numeric(lst))
      print(maxCells)
      cur<-batchSample(cur,total=maxCells,maxPerc=1.0)
    }
    
    print("Perform DE")
    cur@meta.data["nGene"]=scale(cur@meta.data[,"nFeature_RNA"])
    #cur@meta.data["XIST"]=scale(cur@assays$RNA@data["XIST",]) # determine sex
    
    mrk<-lm.pert(cur@meta.data[,c("perturbations","batch","nGene","Cluster1")],form,useBatch=T,dep_var = "Cluster1",ref = ref,numRep=numRep,minPert=minPert)
    
    
    ret[[i]]=mrk
    scores[[i]]=cur@meta.data
  }
  if(getScore)
  {
    return(scores)
  }
  return(ret)
  
}


drawResults<-function(out,meanCenter=F,centerBy="",ref="GFP")
{
  for(i in names(out)){out[[i]]["Module"]=i;out[[i]]=add_row(out[[i]],perturbation=ref,Module=i,Effect_Size=0, padj=0.99);if(meanCenter){out[[i]]["Effect_Size"]=out[[i]][,"Effect_Size"]-mean(out[[i]][,"Effect_Size"])}}
  
  toPlot=c()
  
  
  for(i in out){toPlot=rbind(toPlot,i)}
  
  #p=ggplot(toPlot,aes(x=Module,y=perturbation,fill=Effect_Size))+geom_tile()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradient2(low="blue",high="red",mid="white")
  p=ggplot(toPlot,aes(x=factor(Module, levels = intersect(celltype_order, Module)),y=factor(perturbation, levels = c("NTC","ADARB2","ADARB1","ADAR1","AAVS1")),fill=Effect_Size,size=-log(padj,10)))+geom_point(shape=21)+scale_fill_gradient2(low="blue",mid="white",high="red")+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(p)
  
}

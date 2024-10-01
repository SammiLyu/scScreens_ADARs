### EMD to measure KO screen effects - whether each guide has similar effects on changing transcriptome?
### code adopted from https://github.com/yanwu2014/teratoma-analysis-code/blob/fd1be4e669027236e20a7d85e1cccf4d9a0743d2/Figure4/03_disease_screen_cluster_enrich.R
library(perturbLM)
library(swne)

# one command from perturbLM is renamed
GenotypeClusterCounts <- function(group.list.1, group.list.2) {
  .Deprecated("GroupClusterCounts", package="perturbLM",
              msg = "Use GroupClusterCounts instead",
              old = as.character(sys.call(sys.parent()))[1L])
  df <- matrix(0, length(group.list.1), length(group.list.2))
  rownames(df) <- names(group.list.1)
  colnames(df) <- names(group.list.2)
  
  for(i in 1:length(group.list.1)) {
    group.1.cells <- group.list.1[[i]]
    for(j in 1:length(group.list.2)) {
      group.2.cells <- group.list.2[[j]]
      df[i,j] <- length(intersect(group.1.cells, group.2.cells))
    }
  }
  return(df)
}

load("Figure3_ADARKO_R1/KO_4ter_project_SO_annot_collapsed_20230201.rda")
KO$guide <- sapply(strsplit(KO$guide_id, split = "_|-"), "[[", 1) ; table(KO$guide)
celltype_order <- c("Neural-Progenitors","Neurons","Retinal-Epi","SCP","Gut-Epi","Hematopoietic",
                    "Adipogenic-MSC-Fib","Chondrogenic-MSC-Fib","Cycling-MSC-Fib","MSC-Fib","MyoFib","Muscle","Pericytes","Smooth-Muscle",
                    "Kidney-Prog")

## Filter away NAs and small clusters
idents <- factor(KO$annot_v1, levels = celltype_order)
idents <- droplevels(idents[idents %in% names(table(idents)[table(idents) > 200])])
table(idents); length(table(idents))
length(idents)

## Load guides
min.cells <- 100
meta_tmp <- split(as.data.frame(KO@meta.data[names(idents),]), f = as.data.frame(KO@meta.data[names(idents),])$guide_id)
guide.dict <- list()
for (i in names(meta_tmp)){
  guide.dict[[i]] <- row.names(meta_tmp[[i]])
}
guide.dict <- guide.dict[sapply(guide.dict, length) > min.cells]
# since current guide assignment is mutually exclusive - cells with KO guide do NOT have NTC, no need to remove any NTC cells that also got a true guide
ntc.guide.dict <- guide.dict[grepl("NTC", names(guide.dict))]
guide.dict[["NTC"]] <- unique(unlist(ntc.guide.dict, F, F))
guide.dict[["ADARB1_1"]] <- NULL # remove this outlier guide
sapply(guide.dict, length) ; sum(sapply(guide.dict, length)) # one additional NTC cb list in the list of guide.dict
# 13513+9030

## Make gene dictionary
#genes <- unique(sapply(names(guide.dict), ExtractField, field = 1, delim = "-"))
names(table(KO$guide))
gene.dict <- lapply(names(table(KO$guide)), function(g) {
  unique(unlist(guide.dict[grepl(g, names(guide.dict))]))
})
names(gene.dict) <- names(table(KO$guide))
gene.dict <- gene.dict[sapply(gene.dict, length) > min.cells]
sapply(gene.dict, length) ; sum(sapply(gene.dict, length)) # gene.dict already only has 1 NTC list
# 13513

## Compute cluster distances
umap.emb <- Embeddings(KO, "ref.umap")[names(idents),]
cl.emb <- apply(t(umap.emb), 1, function(x) tapply(x, idents, mean))
cl.dist <- as.matrix(dist(cl.emb))

## Compute guide level EMD
guide.emd <- computeEMD(guide.dict, idents, cl.dist)
guide.emd.ntc <- sort(guide.emd["NTC",])
ggBarplot(guide.emd.ntc, fill.color = "skyblue")

## Compute gene level EMD
gene.emd <- computeEMD(gene.dict, idents, cl.dist)
gene.emd.ntc <- sort(gene.emd["NTC",])
pdf('Figure3_ADARKO_R1/KO_4ter_EMD_genelevel_20230220.pdf')
ggBarplot(gene.emd.ntc, fill.color = "skyblue")
dev.off()
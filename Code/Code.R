
#单细胞分析----
source("SeuratV4CellAnnotateFunc.R")
library(Seurat)

count <- read.csv("F:\\单细胞联合孟德尔随机化探讨帕金森病脑血管周围巨噬细胞相关分子机制的研究\\数据\\scRNA-GSE157783\\IPDCO_hg_midbrain_UMI.tsv",
                  header = T,sep = "\t")
count$id <- rownames(count)
id <- read.csv("F:\\单细胞联合孟德尔随机化探讨帕金森病脑血管周围巨噬细胞相关分子机制的研究\\数据\\scRNA-GSE157783/IPDCO_hg_midbrain_genes.tsv",
               header = T,sep = "\t")
id$id <- rownames(id)

counts <- merge(id,count,by="id")
counts <- counts[,-c(1,3:8)]

library('org.Hs.eg.db')
library('clusterProfiler')
ENTREZID<- bitr(counts[,1], fromType = "ENSEMBL", toType=c("SYMBOL"),OrgDb = org.Hs.eg.db)

counts2 <- merge(ENTREZID,counts,by.x="ENSEMBL",by.y="gene")
counts2 <- counts2[,-1]

counts3 <- counts2[!duplicated(counts2[,1]), ]  
rownames(counts3) <- counts3[,1]
counts3 <- counts3[,-1]

Seurat_object <- CreateSeuratObject(
  counts = counts3,
  min.cells = 3,
  min.features = 200)

sample <- read.csv("F:\\单细胞联合孟德尔随机化探讨帕金森病脑血管周围巨噬细胞相关分子机制的研究\\数据\\scRNA-GSE157783\\IPDCO_hg_midbrain_cell.tsv",
                   header = T,sep = "\t",row.names = 1)
sample <- sample[,2]

Seurat_object <- AddMetaData(Seurat_object,
                             sample,
                             col.name = "patient" 
                              )  
Seurat_object$orig.ident <- Seurat_object$patient

dim(Seurat_object)

Seurat_object$group <- ifelse(substr(Seurat_object$orig.ident,2,2) %in% c(1:6),"Control","Disease")
table(Seurat_object$group,Seurat_object$orig.ident)


load("DatabaseImmuneCell.se.Rdata")
load("BlueprintEncode.Rdata")
load("MonacoImmune.Rdata")
load("HumanPrimaryCellAtlas.Rdata")
load("NovershternHematopoietic.Rdata")


scRNA <- scRNAAutoAnno(
  Path = "./",
  SeuratObject = Seurat_object,
  Multi = TRUE,
  ref = list(
    BlueprintEncode,
    DatabaseImmuneCell.se,
    NovershternHematopoietic,
    HumanPrimaryCellAtlas,
    MonacoImmune
  ),
  labels = list(
    BlueprintEncode$label.main,
    DatabaseImmuneCell.se$label.main,
    NovershternHematopoietic$label.main,
    HumanPrimaryCellAtlas$label.main,
    MonacoImmune$label.main
  ),
  Idents = "RNA_snn_res.0.2"
)



set.seed(2024)
library(Seurat)
library(ggplot2)
library(stringr)
library(scRNAtoolVis)
ggsave_fun <- function(filename,plot = last_plot(),plot_Device=c('.png','.pdf'),width=7,height=7,dpi=600, ...) {
  filenames <- paste0(filename,plot_Device)
  for (filename_tmp in filenames) {
    ggsave(filename = filename_tmp,plot=plot,width = width,height = height,dpi=dpi,...)
  }
}

data.merge.harmony <- readRDS("../.Cluster/data.merge.harmony.2000.rds")


Idents(data.merge.harmony) <- 'RNA_snn_res.0.2'
sub.Fibro.markers <- FindAllMarkers(data.merge.harmony,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sub.Fibro.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(sub.Fibro.markers,"Findall.markers.RNA_snn_res.0.2.csv")

table(data.merge.harmony$RNA_snn_res.0.2)

features <-c( 
  "MOBP","MOG","PLP1",      
  "GFAP","AQP4","RFX4",  
  "CSF1R","C1QB","CD163",   
  "PCDH15","VCAN","TNR", 
  "MYT1L","STMN2","SYT1",     
  "CLDN5","FLT1","VWF",   
  "RGS5","PDGFRB","MYO1B",    
  "CFAP299","FOXJ1","CFAP54"   
)

p <- DotPlot(data.merge.harmony,features = features,cols = "RdYlBu",group.by = "RNA_snn_res.0.2"
)+ scale_size_continuous(range = c(0, 10))+ theme(
  panel.border = element_rect(colour = 'black'),
  axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ),
  legend.position = "right"
);p

cluster.label <- data.merge.harmony$RNA_snn_res.0.2
cluster.label <- gsub("^0$", "Oligodendrocytes", cluster.label)
cluster.label <- gsub("^1$", "Oligodendrocytes", cluster.label)
cluster.label <- gsub("^2$", "Astrocytes", cluster.label)
cluster.label <- gsub("^3$", "Perivascular macrophages", cluster.label)
cluster.label <- gsub("^4$", "Oligodendrocyte precursor cells", cluster.label)
cluster.label <- gsub("^5$", "Neurons", cluster.label)
cluster.label <- gsub("^6$", "Endothelial cells", cluster.label)
cluster.label <- gsub("^7$", "Pericytes", cluster.label)
cluster.label <- gsub("^8$", "Oligodendrocytes", cluster.label)
cluster.label <- gsub("^9$", "Ependymal cells", cluster.label)
cluster.label <- gsub("^10$", "Neurons", cluster.label)
cluster.label <- gsub("^11$", "Oligodendrocytes", cluster.label)
cluster.label <- gsub("^12$", "Neurons", cluster.label)
cluster.label <- gsub("^13$", "Oligodendrocytes", cluster.label)




data.merge.harmony <- AddMetaData(data.merge.harmony, cluster.label, col.name = "cellType_1")
data.merge.harmony$cellType_1 <- factor(data.merge.harmony$cellType_1,c("Oligodendrocytes","Astrocytes",
                                                                        "Perivascular macrophages","Oligodendrocyte precursor cells",
                                                                        "Neurons","Endothelial cells",
                                                                        "Pericytes","Ependymal cells"))

p <- DimPlot(data.merge.harmony,group.by = "cellType_1",label = T,label.size = 3,reduction = "umap");p   
ggsave_fun(filename = paste0("3.CellAnnotate/2.Umap.plot"),plot = p,width = 9,height = 7)

p <- DotPlot(data.merge.harmony,features = features,cols = "RdYlBu",group.by = "cellType_1"
)+ scale_size_continuous(range = c(0, 10))+ theme(
  panel.border = element_rect(colour = 'black'),
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ),
  legend.position = "top",
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.8, "cm"),
);p
ggsave_fun(filename = paste0("3.CellAnnotate/3.Dotplot"),plot = p,width = 15,height = 7)



my_pal2 <- c("#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4","#FCBF6E", "#83AD00","#9ebcda",
             "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31","#F37B7D","#FF8E00")
cellRatio<-cellRatioPlot(object = data.merge.harmony,
                         sample.name = "group",    
                         celltype.name = "cellType_1",
                         flow.curve = 0.5,
                         fill.col = my_pal2)+ theme(
                           axis.text.x = element_text(
                             angle = 45,
                             hjust = 1,
                             vjust = 1
                           ));cellRatio

ggsave_fun(filename = paste0("3.CellAnnotate/4.CellAnnoted.cellType.group.ratio.plot"),plot = cellRatio,width = 7,height = 7)

saveRDS(data.merge.harmony,"CellAnnoted.cellType_1.rds")


#AUCell 评估不同通路在不同细胞上的活性----

library(msigdbr)
h.human <- msigdbr(species = "Homo sapiens", category = "H")
h.names <- unique(h.human$gs_name)
h.sets <- vector("list", length=length(h.names))
names(h.sets) <- h.names
for(i in names(h.sets)){
  h.sets[[i]] <- subset(h.human, gs_name == i, "gene_symbol") %>% unlist(.)
}

map_names = function(seur=NULL, names=NULL){
  
  names = as.list(names)
  genes = c()
  feats = c()
  
  
  data = seur@assays[['RNA']]@data
  meta = seur@meta.data
  
  
  if(!is.null(names)){
    genes = sapply(names, function(a){intersect(a, rownames(data))}, simplify=F)
    feats = sapply(names, function(a){intersect(a, colnames(meta))}, simplify=F)
  }
  
  
  genes = genes[lengths(genes) > 0]
  feats = feats[lengths(feats) > 0]
  
  
  if(length(genes) > 0){
    if(is.null(names(genes))){names(genes) = sapply(genes, paste, collapse='.')}
  }
  
  
  if(length(feats) > 0){
    if(is.null(names(feats))){names(feats) = sapply(feats, paste, collapse='.')}
  }
  
  return(list(genes=genes, feats=feats))
}

score_cells = function(seur=NULL, names=NULL, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL){
  
  data = seur@assays[['RNA']]@data
  meta = seur@meta.data
  if(!is.null(groups)){groups = setNames(groups, colnames(seur@assays[['RNA']]@data))}
  scores = NULL
  
  
  res = map_names(seur=seur, names=names)
  genes = res$genes
  feats = res$feats
  genes.use = unique(do.call(c, genes))
  
  
  if(!is.null(cells.use)){
    data = data[,cells.use,drop=F]
    meta = meta[cells.use,,drop=F]
    if(!is.null(groups)){groups = groups[cells.use]}
  }
  
  group_genes = function(x, method){
    
    if(nrow(x) == 1){return(x[1,,drop=F])}
    if(method == 'sum'){
      t(colSums(x))
    } else if(method == 'mean'){
      t(colMeans(x, na.rm=T))
    } else if(method == 'scale'){
      x = t(scale(t(x)))
      t(colMeans(x, na.rm=T))
    } else if(method == 'scale2'){
      x = t(scale(t(x), center=F))
      t(colMeans(x, na.rm=T))
    } else if(method == 'none'){
      x
    } else {
      stop('Error: invalid combine_genes method')
    }
  }
  
  group_cells = function(x, groups, method){
    
    if(is.null(groups)){return(x)}
    if(method %in% c('n', 'sum')){
      if(method == 'n'){x = x > 0}
      x = t(data.frame(aggregate(t(x), list(groups), sum, na.rm=T), row.names=1))
    } else {
      if(method == 'alpha'){x = x > 0}
      if(method == 'mu'){x[x == 0] = NA}
      x = t(data.frame(aggregate(t(x), list(groups), mean, na.rm=T), row.names=1))
    }
    x[is.na(x)] = 0
    x
  }
  
  names.use = unique(c(names(genes), names(feats)))
  
  
  name_map = sapply(names.use, function(a) c(genes[[a]], feats[[a]]), simplify=F)
  do.flat = all(lengths(name_map) == 1)
  if(do.flat == TRUE){
    genes[['flat']] = do.call(c, genes)
    feats[['flat']] = do.call(c, feats)
    names.iter = 'flat'
    combine_genes = 'none'
  } else {
    names.iter = names.use
  }
  
  backup = scores
  scores = lapply(names.iter, function(name){
    
    if(name %in% names(genes)){
      si = data[genes[[name]],,drop=F]
    } else {
      si = c()
    }
    if(name %in% names(feats)){
      if(is.null(si)){
        si = t(meta[,feats[[name]],drop=F])
      } else {
        si = rBind(si, t(meta[,feats[[name]],drop=F]))
      }
    }
    si = as.matrix(as.data.frame(si))
    si = group_genes(si, method=combine_genes)
    si = group_cells(si, groups=groups, method=group_stat)
    si = data.frame(t(si))
  })
  
  
  if(do.flat == TRUE){
    scores = scores[[1]][,make.names(name_map[names.use]),drop=F]
  } else {
    do.collapse = all(lapply(scores, ncol) == 1)
    if(do.collapse == TRUE){
      scores = as.data.frame(do.call(cbind, scores))
    }
  }
  
  names.use = make.names(names.use)
  names(scores) = names.use
  
  if(!is.null(backup)){
    if(is.data.frame(scores)){
      scores = cbind.data.frame(scores, backup)
    } else {
      scores = c(scores, backup)
    }
  }
  
  if(nrow(scores) == 0){scores = NULL}
  return(scores)
}
scRNA <- readRDS("../CellAnnoted.cellType_1.rds")
Idents(scRNA)<-scRNA$cellType_1
table(scRNA$cellType_1)


AUC_scores <- score_cells(subset(scRNA), names = h.sets)

library(corrplot)
library(pheatmap)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)

AUC_scores <- data.frame(t(AUC_scores),check.names = F)
identical(colnames(scRNA),colnames(AUC_scores))
scRNA[["HALLMARK"]] <- CreateAssayObject(AUC_scores,min.cells = 0,min.features = 0)
rownames(scRNA@assays$HALLMARK@data) <- gsub("-","_",rownames(scRNA@assays$HALLMARK@data))
pdf("1.HALLMARK_AverageHeatmap.pdf",width = 10,height = 12)
scRNAtoolVis::AverageHeatmap(scRNA,assays = "HALLMARK",markerGene = rownames(scRNA@assays$HALLMARK@data))
dev.off()



#CellCall----
library(cellcall)
library(ggplot2)
library(Seurat)

pbmc <- readRDS("../CellAnnoted.cellType_1.rds")
pbmc_cmuni <- CreateObject_fromSeurat(Seurat.object = pbmc,
                                      slot = "counts",
                                      cell_type = "cellType_1", data_source = "UMI",
                                      scale.factor = 10^4,
                                      Org = "Homo sapiens")

pbmc_cmuni <- readRDS("CellCall/pbmc_cmuni.rds")
table(pbmc$cellType_1)

(PathSet <- grep("^Perivascular macrophages", colnames(pbmc_cmuni@data$expr_l_r_log2_scale), value = T))

pathway.hyper.list <- lapply(PathSet, function(i){
  print(i)
  tmp <- getHyperPathway(data = pbmc_cmuni@data$expr_l_r_log2_scale, object = pbmc_cmuni, cella_cellb = i, Org="Homo sapiens")
  return(tmp)
})

myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=PathSet)
myPub.df$cc <- gsub("\\.", "\\_", myPub.df$cc)
p <- plotBubble(myPub.df)
ggsave(p, filename = paste0("1.Cell - pathway map.pdf"), width = 10, height = 8)

my_pal2 <- c("#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4","#FCBF6E")

cell_color <- data.frame(color=my_pal2, stringsAsFactors = FALSE)
rownames(cell_color) <- unique(pbmc$cellType_1)


pdf("2.Circle.pdf",width = 10)
ViewInterCircos(object = pbmc_cmuni@data$expr_l_r_log2_scale, font = 2, 
                cellColor = cell_color,
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = T)

dev.off()


#鉴定脑血管周围巨噬细胞亚型差异基因----

library(Seurat)
library(ggplot2)
library(stringr)
pbmcSub <- readRDS("../CellAnnoted.cellType_1.rds")

Idents(pbmcSub) <- pbmcSub$cellType_1

logFCfilter= 0
adjPvalFilter=0.05
pbmcSub.markers=FindAllMarkers(object = pbmcSub,
                               only.pos = FALSE,
                               min.pct = 0.1,  
                               logfc.threshold = logFCfilter)
sig.cellMarkers=pbmcSub.markers[(abs(as.numeric(as.vector(pbmcSub.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmcSub.markers$p_val_adj))<adjPvalFilter),]

write.csv(sig.cellMarkers, "pbmcSub.markers.csv", quote = F)

sig.cellMarkers <- read.csv("pbmcSub.markers.csv",row.names = 1,header = T)

library(ggrepel)
set.seed(1)
library(ggplot2)
library(scRNAtoolVis)
p <- markerVocalno(markers = sig.cellMarkers,
                   topn = 5,
                   log2FC = 1,
                   labelCol = ggsci::pal_npg()(9))

p
ggsave(filename = "1.Marker_gene_pointplot.pdf", plot = p, width = 12, height =5.5)




#孟德尔随机化分析----

library(TwoSampleMR)
library(MRPRESSO)
library(GagnonMR)
library(ieugwasr)
library(ggplot2)
library(stringr)
Gwasid <- "ebi-a-GCST90018894"
CutP1 <- 1e-5

CutP2 <- 5e-5

CandGen <- read.table("genes-1-570.txt",sep="\t",check.names=F,header=T)

CandGen <- unique(CandGen$gene);length(CandGen)
gene_human <- read.table("mart_export.txt",sep="\t",header=T,check.names=F)
outTab0 <- outTab1 <- outTab2 <- KeyGen <- NULL
SensiLoo <- 1;HeteroIVW <- 0
for (i in CandGen) {
  
  
  if(nrow(gene_human[gene_human$`Gene name` == i,]) != 0){
    
    dir.create(i)
    
    n2 <- n1 <- 1
    repeat{
      try({
        expd1 <- extract_instruments(outcomes = paste0('eqtl-a-',gene_human[gene_human$`Gene name` == i,]$`Gene stable ID`), p1 = CutP1,
                                     access_token = NULL, clump = T, p2 = CutP2, r2 = 0.001, kb = 10000)
      })
      if(!exists("expd1")){
        Sys.sleep(1)
        try({
          expd0 <- extract_instruments(outcomes = paste0('eqtl-a-',gene_human[gene_human$`Gene name` == i,]$`Gene stable ID`), p1 = 1,
                                       access_token = NULL, clump = T, p2 = CutP2, r2 = 0.001, kb = 10000)
        })
        if(exists("expd0")){
          if(is.null(expd0)){
            if(n1 == 3){
              break
            }
            n1 <- n1 + 1
          }else{
            if(!(min(expd0$pval.exposure) < CutP1)){
              break
            }else{
              Sys.sleep(1)
              if(n2 == 3){
                break
              }
              n2 <- n2 + 1
            }
          }
        }
      }else{
        if(is.null(expd1)){
          Sys.sleep(1)
          try({
            expd0 <- extract_instruments(outcomes = paste0('eqtl-a-',gene_human[gene_human$`Gene name` == i,]$`Gene stable ID`), p1 = 1,
                                         access_token = NULL, clump = T, p2 = CutP2, r2 = 0.001, kb = 10000)
          })
          if(exists("expd0")){
            if(is.null(expd0)){
              if(n1 == 3){
                break
              }
              n1 <- n1 + 1
            }else{
              if(!(min(expd0$pval.exposure) < CutP1)){
                break
              }else{
                Sys.sleep(1)
                if(n2 == 3){
                  break
                }
                n2 <- n2 + 1
              }
            }
          }
          Sys.sleep(1)
        }else{
          break
        }
      }
    }
    
    
    if(exists("expd1")){
      if(!is.null(expd1)){
        if(nrow(expd1) > 0){
          
          expd1$Symbol <- i
          outTab0 <- rbind(outTab0, expd1)
          n3 <- 1
          repeat{
            try({
              outd1 <- TwoSampleMR::extract_outcome_data(expd1$SNP, Gwasid, access_token = NULL)
            })
            if(exists("outd1")){
              if(is.null(outd1)){
                Sys.sleep(1)
                if(n3 == 3){
                  break
                }
                n3 <- n3 + 1
              }
              if(!is.data.frame(outd1)){
                Sys.sleep(1)
              }else{
                break
              }
            }else{
              Sys.sleep(1)
            }
          }
          
          
          if(!is.null(outd1)){
            
            dat1 <- TwoSampleMR::harmonise_data(expd1, outd1)
            outTab1 <- merge(outTab1, dat1, all = T)
            
            if(nrow(dat1) == 1){
              res <- TwoSampleMR::mr(dat1, method_list = c("mr_wald_ratio"))
              if(nrow(res) > 0){
                
                
                p <- my_mr_scatter_plot(res, dat1)
                pdf(paste0("./",i,"/PointPlot.pdf"),width = 8,height = 8)
                print(p)
                invisible(dev.off())
                
                
                res_OR <- generate_odds_ratios(res)
                write.table(res_OR,file = paste0("./",i,"/odds_ratios.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
              }else{
                repeat{
                  a <- file.rename(paste0("./",i),paste0("./mr_",i))
                  if(a){
                    break
                  }
                }
              }
            }else{
              res <- TwoSampleMR::mr(dat1, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
              if(nrow(res) == 0){
                res <- TwoSampleMR::mr(dat1, method_list = c("mr_wald_ratio"))
                if(nrow(res) > 0){
                  
                  
                  p <- my_mr_scatter_plot(res, dat1)
                  pdf(paste0("./",i,"/PointPlot.pdf"),width = 8,height = 8)
                  print(p)
                  invisible(dev.off())
                  
                  
                  res_OR <- generate_odds_ratios(res)
                  write.table(res_OR,file = paste0("./",i,"/odds_ratios.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                }else{
                  repeat{
                    a <- file.rename(paste0("./",i),paste0("./mr_",i))
                    if(a){
                      break
                    }
                  }
                }
              }else{
                p <- mr_scatter_plot(res, dat1)
                pdf(paste0("./",i,"/PointPlot.pdf"),width = 8,height = 8)
                print(p)
                invisible(dev.off())
                
                
                res_OR <- generate_odds_ratios(res)
                write.table(res_OR,file = paste0("./",i,"/odds_ratios.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                
                
                res_ht <- mr_heterogeneity(dat1)
                write.table(res_ht,file = paste0("./",i,"/heterogeneity.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                
                
                if(nrow(dat1) > 3){
                  res_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                                          SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat1,
                                          NbDistribution = 1000,  SignifThreshold = 0.05, seed = 123456)
                  b <- gsub("<","",res_presso$`MR-PRESSO results`$`Global Test`$Pvalue)
                  save(res_presso,file = paste0("./",i,"/presso_",b,"_.rda"))
                }
                
                
                res_pl <- mr_pleiotropy_test(dat1)
                write.table(res_pl,file = paste0("./",i,"/pleiotropy.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                
                
                res_loo <- mr_leaveoneout(dat1)
                p <- mr_leaveoneout_plot(res_loo)
                pdf(paste0("./",i,"/LeaveoneoutPlot.pdf"),width = 8,height = ((nrow(res_loo)/10) + 3))
                print(p)
                invisible(dev.off())
                write.table(res_loo,file = paste0("./",i,"/leaveoneout.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                
                
                res_single <- mr_singlesnp(dat1)
                p <- mr_forest_plot(res_single)
                pdf(paste0("./",i,"/ForestPlot.pdf"),width = 6.5,height = ((nrow(res_single)/10) + 3))
                print(p)
                invisible(dev.off())
                p <- mr_funnel_plot(res_single)
                pdf(paste0("./",i,"/FunnelPlot.pdf"),width = 8,height = 8)
                print(p)
                invisible(dev.off())
                p <- mr_density_plot(res_single, res)
                pdf(paste0("./",i,"/DensityPlot.pdf"),width = 8,height = 5.5)
                print(p)
                invisible(dev.off())
                write.table(res_single,file = paste0("./",i,"/singlesnp.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
              }
            }
            if(nrow(res) > 0){
              
              res$Symbol <- i
              outTab2 <- rbind(outTab2, res)
              if(min(res$pval) < 0.05){
                
                KeyGen <- c(KeyGen, i)
                
                if(SensiLoo == 1){
                  if(unique(res$nsnp) > 1){
                    if(length(setdiff(res_loo$SNP, "All")) == 0){
                      repeat{
                        a <- file.rename(paste0("./",i),paste0("./lo_",i))
                        if(a){
                          break
                        }
                      }
                      KeyGen <- setdiff(KeyGen, i)
                    }
                  }
                }
                
                if(HeteroIVW == 1){
                  if(unique(res$nsnp) > 1){
                    if(res_ht[res_ht$method == "Inverse variance weighted",]$Q_pval < 0.05){
                      repeat{
                        a <- file.rename(paste0("./",i),paste0("./ht_",i))
                        if(a){
                          break
                        }
                      }
                      KeyGen <- setdiff(KeyGen, i)
                    }
                  }
                }
              }else{
                repeat{
                  a <- file.rename(paste0("./",i),paste0("./pv_",i))
                  if(a){
                    break
                  }
                }
              }
            }
          }else{
            repeat{
              a <- file.rename(paste0("./",i),paste0("./ou_",i))
              if(a){
                break
              }
            }
          }
        }
      }else{
        repeat{
          a <- file.rename(paste0("./",i),paste0("./ex_",i))
          if(a){
            break
          }
        }
      }
    }else{
      repeat{
        a <- file.rename(paste0("./",i),paste0("./ex_",i))
        if(a){
          break
        }
      }
    }
    Sys.sleep(1)
    rm(list = intersect(ls(),c("expd0", "expd1", "outd1", "res", "dat1", "res_presso", "res_OR", "res_ht", "res_pl", "res_loo", "res_single", "p")))
  }
}
Sys.sleep(1)
write.table(outTab0,file = paste0("./0.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab1,file = paste0("./1.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab2,file = paste0("./2.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

KeyGen <- unique(KeyGen);length(KeyGen)
tmp <- outTab2[outTab2$Symbol %in% KeyGen & outTab2$method == "Inverse variance weighted",]
if(nrow(tmp[which(tmp$pval < 0.05),]) > 2){
  cat("IVW pval < 0.05")
  KeyGen <- tmp[which(tmp$pval < 0.05),]$Symbol
}
length(KeyGen)
save(KeyGen,file = "./KeyGen.rda")#load("./KeyGen.rda")
Sys.sleep(1)

outTab0 <- outTab0[,c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "Symbol")]
colnames(outTab0) <- str_split_fixed(colnames(outTab0),"\\.",2)[,1]
write.table(outTab0,file = paste0("./1.table.xls"),sep = "\t",row.names = F,col.names = T,quote = F)
Sys.sleep(1)



files <- list.dirs("./", full.names = F, recursive = F)
files <- files[!(str_split_fixed(files, "_", 2)[,1] %in% c("ex"))]
dir.create("PointPlot");dir.create("LeaveoneoutPlot");dir.create("ForestPlot");dir.create("FunnelPlot");dir.create("DensityPlot")
Sys.sleep(1)
for (i in files) {
  
  file.copy(paste0("./",i,"/PointPlot.pdf"), paste0("./PointPlot/",gsub("pv_|mr_","",i),".pdf"), overwrite = T)
}
Sys.sleep(1)
files <- KeyGen
for (i in files) {
  
  file.copy(paste0("./",i,"/odds_ratios.txt"), paste0("./PointPlot/",i,".txt"), overwrite = T)
  file.copy(paste0("./",i,"/PointPlot.pdf"), paste0("./PointPlot/",i,".pdf"), overwrite = T)
  file.copy(paste0("./",i,"/LeaveoneoutPlot.pdf"), paste0("./LeaveoneoutPlot/",i,".pdf"), overwrite = T)
  file.copy(paste0("./",i,"/ForestPlot.pdf"), paste0("./ForestPlot/",i,".pdf"), overwrite = T)
  file.copy(paste0("./",i,"/FunnelPlot.pdf"), paste0("./FunnelPlot/",i,".pdf"), overwrite = T)
  file.copy(paste0("./",i,"/DensityPlot.pdf"), paste0("./DensityPlot/",i,".pdf"), overwrite = T)
}




#共定位----
library(locuscomparer)
library(TwoSampleMR)
library(GagnonMR)
library(ieugwasr)
library(ggplot2)
library(coloc)
Gwasid <- "ebi-a-GCST90018894"
DirSMR <- "C:/Users/ROG/Desktop/共定位/smr-1.3.1-win-x86_64/"
DireQTL <- "C:/Users/ROG/Desktop/共定位/eQTLGen/"
GwasCaseNum <- 2638;GwasControlNum <- 477380

(CandGen <- c("ARHGAP22","CTSB","CYTH4","GPNMB","HLA-DRA","ITGAX","KYNU",
              "LCP1","LPAR6","LRRK2","OSBPL3","PELI1","PREX1","RCSD1",
              "RIN3","SAT1","TBC1D1","USP15"))


gene_human <- read.table("mart_export.txt",sep="\t",header=T,check.names=F)

outTab <- list()
for (i in CandGen) {
  
  if(nrow(gene_human[gene_human$`Gene name` == i,]) != 0){
    
    for (j in gene_human[gene_human$`Gene name` == i,]$`Gene stable ID`) {
      
      
      if(nrow(gene_human[gene_human$`Gene name` == i,]) > 1){
        dir.create(paste0(i,"_",j))
        SaveDir <- paste0(i,"_",j)
      }else{
        dir.create(i)
        SaveDir <- i
      }
      
      
      shell(paste0(DirSMR,"smr-1.3.1-win.exe --beqtl-summary ",
                   DireQTL,"eQTLGen --query 5e-2 --gene ",
                   j," --out ./",SaveDir,"/myquery"))
      Sys.sleep(1)
      if(file.exists(paste0(SaveDir,"/myquery.txt"))){
        
        eqtl <- read.table(paste0(SaveDir,"/myquery.txt"),sep="\t",header=T,check.names=F)
        eqtl$maf <- ifelse(eqtl$Freq > 0.5, (1-eqtl$Freq), eqtl$Freq)
        head(eqtl)
        write.table(eqtl,row.names = F,file=paste0(SaveDir,"/eqtl.txt"),sep="\t",quote=F)
        write.table(eqtl[,c("SNP", "p")],row.names = F,file=paste0(SaveDir,"/eqtl_in.txt"),sep="\t",quote=F)
        repeat{
          try({
            gwas <- ieugwasr::associations(variants=paste0(unique(eqtl$Chr), ":",min(eqtl$BP),"-",max(eqtl$BP)), id = Gwasid)
          })
          if(is.data.frame(gwas)){
            break
          }else{
            Sys.sleep(1)
          }
        }
        gwas <- as.data.frame(gwas)
        gwas <- gwas[gwas$rsid != "",]
        gwas$maf <- ifelse(gwas$eaf > 0.5, (1-gwas$eaf), gwas$eaf)
        gwas$varbeta <- gwas$se^2
        gwas <- gwas[order(gwas$rsid, gwas$maf, decreasing = T),]
        gwas <- gwas[!duplicated(gwas$rsid, fromLast = F),]
        head(gwas)
        write.table(gwas,row.names = F,file=paste0(SaveDir,"/gwas.txt"),sep="\t",quote=F)
        write.table(gwas[,c("rsid","p")],row.names = F,file=paste0(SaveDir,"/gwas_in.txt"),sep="\t",quote=F)
        
        
        p <- locuscompare(in_fn1=paste0("./",SaveDir,"/gwas_in.txt"), in_fn2=paste0("./",SaveDir,"/eqtl_in.txt"), 
                          title1="GWAS", title2="eQTL", marker_col1= "rsid", 
                          pval_col1="p", marker_col2="SNP", pval_col2="p")
        pdf(file=paste0("./",SaveDir,"/locuscompare.pdf"),width=7,height=6)
        print(p)
        dev.off()
        
        repeat{
          try({
            a <- ieugwasr::gwasinfo(id = Gwasid, access_token = NULL)
          })
          if(nrow(a) == 0){
            Sys.sleep(1)
          }else{
            break
          }
        }
        if(!is.element("ncase",colnames(a)) | !is.element("ncontrol",colnames(a))){
          if(!is.element("ncase",colnames(a))){
            a$ncase <- GwasCaseNum
          }
          if(!is.element("ncontrol",colnames(a))){
            a$ncontrol <- GwasControlNum
          }
        }
        colnames(eqtl)[1] <- "rsid"
        if(min(eqtl$p) <= 0){
          eqtl[which(eqtl$p <= 0),]$p <- sort(unique(eqtl$p[!(eqtl$p <= 0)]), decreasing = F)[1] / 10
        }
        if(min(gwas$p) <= 0){
          gwas[which(gwas$p > 1),]$p <- sort(unique(gwas$p[!(gwas$p <= 0)]), decreasing = F)[1] / 10
        }
        inputDat <- merge(eqtl, gwas, by = "rsid", all = F, suffixes = c("_eqtl", "_gwas"))
        
        if(length(inputDat$maf_gwas[!is.na(inputDat$maf_gwas)]) == 0){
          
          Res <- coloc.abf(dataset1 = list(pvalues = inputDat$p_eqtl, type = "quant", N = 31864, snp = inputDat$rsid, MAF=inputDat$maf_eqtl),
                           dataset2 = list(pvalues = inputDat$p_gwas, type = "cc", s=a$ncase/unique(gwas$n), N = unique(gwas$n), snp = inputDat$rsid, beta = inputDat$beta, varbeta = inputDat$varbeta))
        }else{
          Res <- coloc.abf(dataset1 = list(pvalues = inputDat$p_eqtl, type = "quant", N = 31864, snp = inputDat$rsid, MAF=inputDat$maf_eqtl),
                           dataset2 = list(pvalues = inputDat$p_gwas, type = "cc", s=a$ncase/unique(gwas$n), N = unique(gwas$n), snp = inputDat$rsid, MAF=inputDat$maf_gwas))
        }
        
        tmp <- Res$results
        write.table(tmp[order(tmp$SNP.PP.H4, decreasing = T),],row.names = F,col.names = T,file=paste0(SaveDir,"/Colocalization.",a$ncase,".",a$ncontrol,".txt"),sep="\t",quote=F)
        write.table(as.data.frame(Res$summary),row.names = T,col.names = F,file=paste0(SaveDir,"/Colocalization.summary.txt"),sep="\t",quote=F)
        write.table(as.data.frame(Res$priors),row.names = T,col.names = F,file=paste0(SaveDir,"/Colocalization.priors.txt"),sep="\t",quote=F)
        outTab[[SaveDir]] <- Res
      }else{
        unlink(SaveDir, recursive = T)
      }
    }
    
    if(nrow(gene_human[gene_human$`Gene name` == i,]) > 1){
      if(length(grep(paste0("^",i,"_"),list.dirs("./",full.names = F,recursive = F),value = T)) == 1){
        file.rename(grep(paste0("^",i,"_"),list.dirs("./",full.names = F,recursive = F),value = T), i)
      }
    }
  }
}
save(outTab,file = "./outTab.rda")#load("./outTab.rda")



CutNum <- 0.90
KeyGen <- NULL
ResTab <- NULL
for (i in names(outTab)) {
  
  a <- sort(outTab[[i]]$results$SNP.PP.H4, decreasing = T)
  if(a[1] > CutNum){
    
    SaveDir <- grep(unlist(strsplit(i,"_"))[1],list.dirs("./",full.names = F,recursive = F),value = T)
    if(length(SaveDir) == 1){
      file.copy(paste0("./",SaveDir,"/locuscompare.pdf"), paste0("./",SaveDir,".pdf"), overwrite = T)
    }else{
      for (j in SaveDir) {
        file.copy(paste0("./",j,"/locuscompare.pdf"), paste0("./",j,".pdf"), overwrite = T)
      }
    }
    tmp <- outTab[[i]]$results
    tmp <- tmp[which(tmp$SNP.PP.H4 > CutNum),]
    tmp$Symbol <- i
    ResTab <- rbind(ResTab, tmp)
    KeyGen <- c(KeyGen, unlist(strsplit(i,"_"))[1])
  }
}
length(KeyGen)
save(KeyGen,file = "./KeyGen.rda")
write.table(ResTab,row.names = F,col.names = T,file=paste0("./Colocalization.result.xls"),sep="\t",quote=F)



#反向MR----


library(TwoSampleMR)
library(MRPRESSO)
library(stringr)
library(dplyr)
ExpSet <- "ebi-a-GCST90018894"

OutSet <- c("eqtl-a-ENSG00000100599")  
CutP1 <- 1e-5
CutP2 <- 5e-5
LimitF <- 0;

ResDirName <- "Result"
dir.create(ResDirName)
for (i in ExpSet) {
  
  n <- 1
  repeat{
    try({
      expd1 <- TwoSampleMR::extract_instruments(i, p1 = CutP1, p2 = CutP2, access_token = NULL,
                                                clump = T, r2 = 0.001, kb = 10000)
    })
    if(exists("expd1")){
      if(is.null(expd1)){
        Sys.sleep(1)
        if(n == 3){
          break
        }
        n <- n + 1
      }else{
        if(nrow(expd1) == 0){
          Sys.sleep(1)
        }else{
          break
        }
      }
    }else{
      Sys.sleep(1)
    }
  }
  
  if(!is.null(expd1)){
    
    expd1$F.value <- (expd1$beta.exposure)^2 / (expd1$se.exposure)^2
    if(LimitF == 1){
      expd1 <- expd1[which(expd1$F.value > CutF),]
    }
    if(nrow(expd1) > 0){
      
      write.table(expd1,file = paste0("./",ResDirName,"/00.",i,".txt"),sep = "\t",row.names = F,col.names = T,quote = F)
      for (j in OutSet) {
        
        dir.create(paste0("./",ResDirName,"/",i,".",j))
        
        repeat{
          try({
            outd1 <- TwoSampleMR::extract_outcome_data(expd1$SNP, j, access_token = NULL)
          })
          if(exists("outd1")){
            if(!is.data.frame(outd1)){
              Sys.sleep(1)
            }else{
              break
            }
          }else{
            Sys.sleep(1)
          }
        }
        if(!is.null(outd1)){
          
          dat1 <- TwoSampleMR::harmonise_data(expd1, outd1)
          write.table(dat1,file = paste0("./",ResDirName,"/",i,".",j,"/01.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
          
          res <- TwoSampleMR::mr(dat1,
                                 method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
          if(nrow(res) > 0){
            
            p <- mr_scatter_plot(res, dat1)
            pdf(paste0("./",ResDirName,"/",i,".",j,"/1.PointPlot.pdf"),width = 8,height = 8)
            print(p)
            invisible(dev.off())
            write.table(res,file = paste0("./",ResDirName,"/",i,".",j,"/02.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
            
            
            res_OR <- generate_odds_ratios(res)
            write.table(res_OR,file = paste0("./",ResDirName,"/",i,".",j,"/03.table - odds_ratios.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
            
            
            res_ht <- mr_heterogeneity(dat1)
            write.table(res_ht,file = paste0("./",ResDirName,"/",i,".",j,"/04.table - heterogeneity.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
            
            
            res_pl <- mr_pleiotropy_test(dat1)
            write.table(res_pl,file = paste0("./",ResDirName,"/",i,".",j,"/05.table - pleiotropy.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
            
            
            res_loo <- mr_leaveoneout(dat1)
            p <- mr_leaveoneout_plot(res_loo)
            pdf(paste0("./",ResDirName,"/",i,".",j,"/5.LeaveoneoutPlot.pdf"),width = 8,height = ((nrow(res_loo)/10) + 3))
            print(p)
            invisible(dev.off())
            write.table(res_loo,file = paste0("./",ResDirName,"/",i,".",j,"/06.table - leaveoneout.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
            
            
            res_single <- mr_singlesnp(dat1)
            p <- mr_forest_plot(res_single)
            pdf(paste0("./",ResDirName,"/",i,".",j,"/3.ForestPlot.pdf"),width = 8,height = ((nrow(res_single)/10) + 3))
            print(p)
            invisible(dev.off())
            p <- mr_funnel_plot(res_single)
            pdf(paste0("./",ResDirName,"/",i,".",j,"/2.FunnelPlot.pdf"),width = 8,height = 8)
            print(p)
            invisible(dev.off())
            p <- mr_density_plot(res_single, res)
            pdf(paste0("./",ResDirName,"/",i,".",j,"/4.DensityPlot.pdf"),width = 8,height = 5.5)
            print(p)
            invisible(dev.off())
            write.table(res_single,file = paste0("./",ResDirName,"/",i,".",j,"/07.table - singlesnp.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
            
            Sys.sleep(1)
            rm(list = intersect(ls(), "expd1"))
          }
        }
      }
    }
  }
}


library(forestplot)

rt <- read.table(file = "OR.txt",sep = "\t",header = T)   
HR=rt[,10:12]
print(HR)

hr=sprintf("%.3f",HR[,"or"])
hrLow=sprintf("%.3f",HR[,"or_lci95"])
hrHigh=sprintf("%.3f",HR[,"or_uci95"])
pVal=rt[,"pval"]
method=rt[,"method"]
exposure=rt[,"exposure"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
clrs <- fpColors(box="red",line="darkblue", summary="royalblue", zero = "lightgray")

tabletext <- list(c("exposure", rt$exposure),
                  append("outcome", rt[,"outcome"]),
                  append("Method",method),
                  append("Nsnp", rt[,"nsnp"]),
                  append("abs(B)", sprintf("%.4f",abs(rt[,"b"]))),
                  append("OR(95%CI)",paste0(hr,"(",hrLow,"-",hrHigh,")")),
                  append("Pvalue", pVal))


p <- forestplot(tabletext,
                rbind(rep(NA, 3), HR),
                col=clrs, graph.pos = 4,
                graphwidth=unit(40, "mm"),
                xlog=T, zero = 1, lwd.zero = 3.5,
                lwd.ci=4, boxsize=0.10,  
                xlab="Odds Ratios",
                txt_gp=fpTxtGp(ticks=gpar(cex=1.1), xlab=gpar(cex=1.25)))
pdf(file=paste0("1.Forest plot.pdf"),width = 15,height = 8,onefile = FALSE)
print(p)
dev.off()



#拟时序分析----

library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(devtools)
library(ggsci)
library(reticulate)
library(SeuratDisk)
library(data.table)

scRNA=readRDS("../CellAnnoted.cellType_1.rds")
scRNA <- subset(scRNA,subset = cellType_1 %in% c("Perivascular macrophages"))
scRNA$cellType_1 <- factor(scRNA$cellType_1,unique(scRNA$cellType_1))
table(scRNA$cellType_1)
pbmcSub <- scRNA



logFCfilter=1    
adjPvalFilter=0.05 


sig.markers<-read.csv("pbmcSub.markers.csv",sep = ",",header = T,row.names = 1)



monocle.matrix=GetAssayData(object = pbmcSub, slot = "data", assay = "RNA")

monocle.sample=pbmcSub@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.markers=sig.markers


data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
Sys.sleep(1)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
# save(cds,file = "cds.rda")#load("cds.rda")

cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))

Sys.sleep(1)


rm(list = setdiff(ls(), c("pbmcSub", "cds", "adjPvalFilter", "logFCfilter")))
cds <- reduceDimension(cds,
                       
                       max_components = 2, reduction_method = 'DDRTree')

cds <- orderCells(cds)

Sys.sleep(1)

library(ggsci)

pdf(file="1.Trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()

pdf(file="2.Trajectory.celltype.pdf",width=7,height=6)

plot_cell_trajectory(cds, color_by = "group")
dev.off()

save(cds,file = "cds.rda")#load("../cds.rda")
Sys.sleep(1)

groups=subset(pData(cds),select='State')
pbmcSub=AddMetaData(object=pbmcSub, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
  pbmc.markers=FindMarkers(pbmcSub, ident.1 = i, min.pct = 0.25, group.by = 'group')
  sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
  sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
  write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
  if(nrow(sig.markers) > 0){
    sig.markers$cluster <- i
    geneList[[i]]=sig.markers
  }
}

unionGenes=do.call(rbind,geneList)
unionGenes <- unionGenes %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
unionGenes <- sort(unique(unionGenes$Gene))
write.table(file="05.monocleDiff.union.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)

pdf(file="3.TrajectoryHeatmap.branch.pdf",width=5,height=6)
Sys.sleep(1)

plot_genes_branched_heatmap(cds[unionGenes,], branch_point = 1, show_rownames = F)
dev.off()
load("cds.rda")
keygenes <- read.table("KeyGenes.txt")
keygenes <- keygenes$V1
pdf(file="4.Pseudotime.pdf",width=10,height=6)
cds_subset <- cds[keygenes,]
Sys.sleep(1)

p1 <-plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <-plot_genes_in_pseudotime(cds_subset, color_by = "cellType_1")
p3 <-plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
p1|p2|p3

dev.off()



#GSEA GSVA motif 关键基因与免疫细胞的相关性 关键基因的表达----
source("SeuratV4CellAnnotateFunc.R")
source("KeyGenePipeline.R")

mat <- read.table("GSE20292.txt",header = T,row.names = 1,check.names = F) %>% as.matrix()
mat <- log2(mat+1)

KeyGenePipeline(
  Path = "./",
  mat = mat,
  KeyGene = c("HLA-DRA","OSBPL3","RIN3"),
  species = NA,
  KeyCell = NA,
  logFCFilter = 1,
  adjPvalFilter = 0.05
)


scRNA <- readRDS("CellAnnoted.cellType_1.rds")
scRNA$group <- scRNA$group
scRNA <- scRNAAutoAnno(
  Path = "./",
  SeuratObject = scRNA,
  Multi = TRUE,
  KeyGene = c("HLA-DRA","OSBPL3","RIN3")
)









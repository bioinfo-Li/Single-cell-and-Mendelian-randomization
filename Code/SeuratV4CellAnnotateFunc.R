require(Seurat)
require(ggplot2)
require(harmony)
require(dplyr)
require(reshape2)
require(scRNAtoolVis)
require(ggpubr)
require(clusterProfiler)
options(dplyr.print_max = 1e9)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(123456)

ggsave_fun <- function(filename, plot = last_plot(), plot_Device = c(".pdf"), width = 7, height = 7, dpi = 600, ...) {
  filenames <- paste0(filename, plot_Device)
  for (filename_tmp in filenames) {
    ggsave(filename = filename_tmp, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
}

paramSweep_V4 <- function(seu, PCs = 1:10, sct = FALSE, num.cores = 1) {
  require(Seurat)
  require(fields)
  require(parallel)
  
  pK <- c(0.0005, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)

  
  min.cells <- round(nrow(seu@meta.data) / (1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]

  
  orig.commands <- seu@commands

  
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }

  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }

  
  
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
      FUN = parallel_paramSweep,
      n.real.cells,
      real.cells,
      pK,
      pN,
      data,
      orig.commands,
      PCs,
      sct, mc.cores = num.cores
    )
    stopCluster(cl)
  } else {
    output2 <- lapply(as.list(1:length(pN)),
      FUN = parallel_paramSweep,
      n.real.cells,
      real.cells,
      pK,
      pN,
      data,
      orig.commands,
      PCs,
      sct
    )
  }

  
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }

  
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

doubletFinder_V4 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat)
  require(fields)
  require(KernSmooth)

  
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
    return(seu)
  }

  if (reuse.pANN == FALSE) {
    
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells / (1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    
    orig.commands <- seu@commands

    
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
        normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
        scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
        margin = orig.commands$NormalizeData.RNA@params$margin
      )

      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
        selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
        loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
        clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
        mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
        dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
        num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
        binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
        nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
        mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
        dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff
      )

      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
        features = orig.commands$ScaleData.RNA$features,
        model.use = orig.commands$ScaleData.RNA$model.use,
        do.scale = orig.commands$ScaleData.RNA$do.scale,
        do.center = orig.commands$ScaleData.RNA$do.center,
        scale.max = orig.commands$ScaleData.RNA$scale.max,
        block.size = orig.commands$ScaleData.RNA$block.size,
        min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block
      )

      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
        features = orig.commands$ScaleData.RNA$features,
        npcs = length(PCs),
        rev.pca = orig.commands$RunPCA.RNA$rev.pca,
        weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
        verbose = FALSE
      )
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc() 
    }

    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)

      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }

    
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)

    
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells)) / k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbors_that_are_doublets <- neighbors[neighbors > n_real.cells]
          if (length(neighbors_that_are_doublets) > 0) {
            neighbor_types[i, ] <-
              table(doublet_types1[neighbors_that_are_doublets - n_real.cells]) +
              table(doublet_types2[neighbors_that_are_doublets - n_real.cells])
            neighbor_types[i, ] <- neighbor_types[i, ] / sum(neighbor_types[i, ])
          } else {
            neighbor_types[i, ] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) <- levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, ct]
      }
    }
    return(seu)
  }
}

human2mouse <- function(x) {
  x1 <- substr(x, 1, 1)
  x2 <- substr(x, 2, nchar(x))
  return(paste0(x1, tolower(x2)))
}

getScatterplot <- function(object, gene1, gene2, cor.method = "pearson", jitter.num = 0.15, pos = TRUE) {
  if (!gene1 %in% rownames(object)) {
    print("gene1 was not found")
    if (!gene2 %in% rownames(object)) {
      print("gene2 was not found")
    }
  } else {
    exp.mat <- GetAssayData(object = object, assay = "RNA") %>%
      .[c(gene1, gene2), ] %>%
      as.matrix() %>%
      t() %>%
      as.data.frame()
    if (pos) {
      if (nrow(exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]) > (nrow(exp.mat) * 0.01)) {
        exp.mat <- exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]
      } else {
        exp.mat <- exp.mat[which(exp.mat[, 1] > 0 | exp.mat[, 2] > 0), ]
      }
    }
    colnames(exp.mat) <- c("Var1", "Var2")
    plots <- ggplot(data = exp.mat, mapping = aes_string(x = "Var1", y = "Var2")) +
      geom_smooth(method = "lm", se = T, color = "red", size = 1) +
      stat_cor(method = cor.method) +
      labs(x = gene1, y = gene2) +
      geom_jitter(width = jitter.num, height = jitter.num, color = "black", size = 1, alpha = 1) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.text = element_text(colour = "black", size = 10),
        axis.text = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(size = 1, linetype = "solid", colour = "black"),
        panel.background = element_rect(fill = "white")
      )
    return(plots)
  }
}

capitalize_first <- function(string) {
  words <- strsplit(string, " ")[[1]]
  capitalized_words <- toupper(substring(words, 1, 1))
  rest_of_words <- substring(words, 2)
  capitalized_string <- paste0(capitalized_words, rest_of_words, collapse = " ")
  return(capitalized_string)
}


scRNAAutoAnno <- function(Path = ".", SeuratObject = NA, Multi = TRUE, ref, labels, Idents = "RNA_snn_res.0.2", DoubletFinder = TRUE, species = NA, KeyCell = NA, KeyGene = NA, logFCfilter = 0.585, adjPvalFilter = 0.05, assay = "RNA", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 20, nfeatures = 2000, npcs = 50) {
  cat("1.质量控制\n")
  cat("2.数据标准化\n")
  cat("3.自动化注释\n")
  cat("7.关键基因表达丰度\n")
  cat("\033[0;34m请选择要分析的内容，按序号填写，英文逗号分隔")
  AnalysisIndex <- scan(sep = ",", quiet = TRUE)

  AdvanceIndex <- length(which(AnalysisIndex %in% c(1:12) == FALSE)) 
  if (is.na(KeyCell) & AdvanceIndex > 0) {
    message("是否需要在分析过程中挑选关键细胞？yes/no")
    isKeyCell <- scan(sep = ",", what = "character", quiet = TRUE)
    if (isKeyCell == "yes") {
      cat("1.细胞含量变化最大\n")
      cat("2.Cellchat\n")
      cat("3.亚群贡献度\n")
      cat("\033[0;34m从哪一步分析中挑选关键细胞？\n")
      KeyCellIndex <- scan(sep = ",", quiet = TRUE)
    }
  }

  scRNA <- SeuratObject

  if (!is.na(species)) {
    mouse <- ifelse(species == "mouse", T, F)
  } else {
    mouse <- F
  }

  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  if (mouse) {
    s.genes <- human2mouse(s.genes)
    g2m.genes <- human2mouse(g2m.genes)
  }

  n <- 0
  KeyCells <- KeyCell
  for (step in AnalysisIndex) {
    n <- n + 1

    
    if (step == 1) {
      message("*************************质量控制*************************")

      dir <- paste0(Path, "/.QCFilter")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      grep("^[M,m][T,t]-", rownames(scRNA), value = T)
      grep("^R[P,p][SL,sl]", rownames(scRNA), value = T)
      scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, pattern = "^[M,m][T,t]-")
      scRNA[["percent.ribo"]] <- PercentageFeatureSet(scRNA, pattern = "^R[P,p][SL,sl]") 

      nomt <- ifelse(max(scRNA[["percent.mt"]]) == 0, T, F)
      if (nomt) {
        MT_genes <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")
        if (mouse) {
          MT_genes <- human2mouse(MT_genes)
        }
        MT_genes <- intersect(MT_genes, rownames(scRNA))
        scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, features = MT_genes)
      }

      mt_value <- stats::mad(scRNA$percent.mt)
      lower_mt <- median(scRNA$percent.mt) - 3 * mt_value
      upper_mt <- median(scRNA$percent.mt) + 3 * mt_value
      
      nF_value <- stats::mad(scRNA$nFeature_RNA)
      lower_nF <- median(scRNA$nFeature_RNA) - 3 * nF_value
      upper_nF <- median(scRNA$nFeature_RNA) + 3 * nF_value
      
      nC_value <- stats::mad(scRNA$nCount_RNA)
      lower_nC <- median(scRNA$nCount_RNA) - 3 * nC_value
      upper_nC <- median(scRNA$nCount_RNA) + 3 * nC_value
      
      nR_value <- stats::mad(scRNA$percent.ribo)
      lower_nR <- median(scRNA$percent.ribo) - 3 * nR_value
      upper_nR <- median(scRNA$percent.ribo) + 3 * nR_value

      scRNA <- subset(x = scRNA, subset = nFeature_RNA >= 200 & percent.mt <= upper_mt & nFeature_RNA <= upper_nF & nCount_RNA <= upper_nC & percent.ribo <= upper_nR)
      dim(scRNA)

      if (DoubletFinder) {
        require(DoubletFinder)
        scRNA_list <- SplitObject(scRNA, split.by = "orig.ident")
        newscRMA_list <- list()
        for (i in 1:length(scRNA_list)) {
          data <- scRNA_list[[i]]
          data <- NormalizeData(data)
          data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
          data <- ScaleData(data)
          data <- RunPCA(data)
          data <- RunUMAP(data, dims = 1:20)

          sweep.res.list <- paramSweep_V4(data, PCs = 1:20, sct = FALSE) 
          sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
          bcmvn <- find.pK(sweep.stats)
          p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))

          
          homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters) 
          Doubletrate <- ncol(data) * 8 * 1e-6
          nExp_poi <- round(Doubletrate * ncol(data))
          nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

          
          data <- doubletFinder_V4(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
          colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
          newscRMA_list[[i]] <- data
        }

        saveRDS(newscRMA_list, paste0(dir, "/doubletFinder_V4.rds"))

        plot_list <- list()
        for (i in 1:length(newscRMA_list)) {
          plot_list[[i]] <- DimPlot(newscRMA_list[[i]], group.by = "doublet_info") + ggtitle(newscRMA_list[[i]]$orig.ident[1]) + theme(plot.title = element_text(hjust = 0.5))
        }


        for (i in 1:length(plot_list)) {
          plot_list[[i]] <- plot_list[[i]] + NoAxes() + NoLegend()
        }

        p <- CombinePlots(plot_list, legend = "right")

        ggsave_fun(filename = paste0(dir, "/1.DoubletFinder"), plot = p, width = 12, height = 10)

        if (Multi) {
          scRNA <- merge(newscRMA_list[[1]], newscRMA_list[2:length(newscRMA_list)])
        } else {
          scRNA <- newscRMA_list[[1]]
        }

        scRNA <- subset(scRNA, subset = doublet_info == "Singlet")

        p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 4)
        
        ggsave_fun(filename = paste0(dir, "/2.VlnPlot"), plot = p, width = 20, height = 4)

        plot1 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
          NoLegend()
        plot2 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
          NoLegend()
        plot3 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo", raster = F) +
          NoLegend()

        plot1 + plot2 + plot3
        ggsave_fun(filename = paste0(dir, "/3.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
      } else {
        p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 4)
        
        ggsave_fun(filename = paste0(dir, "/1.VlnPlot"), plot = p, width = 20, height = 4)

        plot1 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
          NoLegend()
        plot2 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
          NoLegend()
        plot3 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo", raster = F) +
          NoLegend()

        plot1 + plot2 + plot3
        ggsave_fun(filename = paste0(dir, "/2.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
      }
    }

    
    if (step == 2) {
      message("*************************数据标准化*************************")

      dir <- paste0(Path, "/.Cluster")

      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      DefaultAssay(scRNA) <- assay

      scRNA <- NormalizeData(object = scRNA)
      scRNA <- CellCycleScoring(object = scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
      scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = nfeatures)
      top10 <- head(VariableFeatures(scRNA), 10)
      plot1 <- VariableFeaturePlot(scRNA)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
      plot1 + plot2
      ggsave_fun(filename = paste0(dir, "/1.VariableFeaturePlot"), plot = plot1 + plot2, width = 14)
      scRNA <- ScaleData(object = scRNA, vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
      scRNA <- RunPCA(scRNA, verbose = T, npcs = npcs, features = VariableFeatures(object = scRNA))
      p <- ElbowPlot(object = scRNA, ndims = npcs)
      ggsave_fun(filename = paste0(dir, "/2.ElbowPlot"), plot = p)

      if (Multi) {
        scRNA <- RunHarmony(object = scRNA, group.by.vars = "orig.ident", assay.use = assay, verbose = FALSE)
        scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:PC, verbose = T)
        
        scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "harmony", verbose = T) 
        scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
        p <- DimPlot(object = scRNA, reduction = "harmony", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/4.HarmonyPlot"), plot = p)
      } else {
        scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:PC, verbose = T)
        
        scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "pca", verbose = T)
        scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
      }

      pdf(paste0(dir, "/7.data.merge.harmony.pdf"))
      
      
      p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = "orig.ident", raster = F) + NoLegend()
      print(p)
      merge.res <- sapply(set.resolutions, function(x) {
        p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = paste0(assay, "_snn_res.", x), raster = F) + NoLegend()
        print(p)
      })
      dev.off()
      saveRDS(scRNA, paste0(dir, "/data.merge.harmony.2000.rds"))
    }

    
    if (step == 3) {
      message("*************************自动化注释*************************")

      dir <- paste0(Path, "/.CellAnnotate")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      
      

      Cluster.dir <- list.files(Path, pattern = "^.Cluster$", all.files = TRUE)

      p <- DimPlot(scRNA, group.by = Idents, label = T, reduction = "umap", raster = F) + NoLegend()
      print(p)
      ggsave_fun(filename = paste0(Path, "/", Cluster.dir, "/5.Umap.Cluster"), plot = p)

      CELL <- SingleR::SingleR(
        test = as.matrix(scRNA@assays$RNA@data),
        ref = ref,
        labels = labels,
        clusters = scRNA@meta.data[[Idents]]
      )
      
      cells <- CELL[, ncol(CELL):1]
      cells <- data.frame(cells)
      cell_id <- cells$pruned.labels
      names(cell_id) <- levels(scRNA@meta.data[[Idents]])
      Idents(scRNA) <- Idents
      
      if (anyNA(cell_id)) {
        scRNA <- RenameIdents(scRNA, cell_id)
        scRNA$cellType_1 <- Idents(scRNA)
        print(cell_id)
        p <- DimPlot(scRNA, group.by = "cellType_1", raster = F)
        print(p)
        message("某些细胞簇没注释到，是否强行为NA的细胞命名为unknow？yes/no")
        exit <- scan(what = "character", sep = ",", quiet = TRUE)
        
        if (exit == "yes") {
          cell_id[which(is.na(cell_id))] <- "unknow"
        }
        
        if (exit == "no") {
          message("某些细胞簇没注释到，具体请看右侧umap图，建议选择其他分辨率或参考数据集重新运行此函数，跳过1.质量控制和2.数据标准化，选择3重新注释。")
          return(scRNA)
        }
      }
      
      if (!anyNA(cell_id)) {
        cat("\033[0;34m由于使用了多个参考数据，参考数据细胞名字不同，请按顺序对下述细胞进行重命名，英文;分隔。\n例：将\033[31mCD8+ T cells,Neutrophils,CD8+ T-cells\033[0;34m更改为\033[31mCD8+ T cells,Neutrophils,CD8+ T cells\033[0;34m以统一细胞名字\n")
        p <- DimPlot(scRNA, group.by = Idents, raster = F,label = T)+NoLegend()
        print(p)
        print(paste0(paste0(names(cell_id),":",as.character(cell_id)), collapse = "; "))
        cell_id <- scan(what = "character", sep = ";", quiet = FALSE)
        cell_id <- gsub(".*:","",cell_id)
        names(cell_id) <- levels(scRNA@meta.data[[Idents]])
        
        scRNA <- RenameIdents(scRNA, cell_id)
        scRNA$cellType_1 <- Idents(scRNA)
        print(table(scRNA$cellType_1))
        p <- DimPlot(scRNA, group.by = "cellType_1", raster = F)
        print(p)
        saveRDS(scRNA, paste0(Path, "/result.rds"))

        Idents(scRNA) <- "cellType_1"
        scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

        write.csv(scRNA.markers, paste0(dir, "/Findall.markers.cellType_1.csv"))

        scRNA.markers %>%
          group_by(cluster) %>%
          top_n(n = 5, wt = avg_log2FC) -> features

        print((features))

        features <- features$gene
        cat("\033[0;34m是否要改写基因？yes/no\n")
        RenameGene <- scan(what = "character", sep = ",", quiet = TRUE)

        if (RenameGene == "yes") {
          print(paste0(as.character(features), collapse = ","))
          cat("\033[0;34m请改写\n")
          features <- scan(what = "character", sep = ",", quiet = FALSE)
        }

        
        my_pal2 <- c(
          "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A"
        )

        umap <- scRNA@reductions$umap@cell.embeddings %>%
          as.data.frame() %>%
          cbind(cellType = scRNA@meta.data$cellType_1)

        
        celltypepos <- umap %>%
          group_by(cellType) %>%
          summarise(
            umap_1 = median(UMAP_1),
            umap_2 = median(UMAP_2)
          )

        p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T) +
          ggrepel::geom_label_repel(aes(x = umap_1, y = umap_2, label = cellType, color = cellType),
            fontface = "bold",
            data = celltypepos,
            box.padding = 0.5, show.legend = FALSE
          )
        print(p)
        ggsave_fun(filename = paste0(dir, "/1.Umap1.plot"), plot = p, width = 10, height = 7)

      
        

        p <- DotPlot(scRNA, features = unique(features), cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 10)) + theme(
          panel.border = element_rect(colour = "black"),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          legend.position = "top",
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.8, "cm"),
        )
        print(p)
        if (length(unique(features)) >= 50) {
          ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 20, height = 9)
        }
        if (length(unique(features)) < 50) {
          ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 17, height = 9)
        }

        tryCatch(
          {
            cellRatio <- cellRatioPlot(
              object = scRNA,
              sample.name = "group",
              celltype.name = "cellType_1",
              flow.curve = 0.5,
              fill.col = my_pal2
            ) + theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1
              )
            )
            print(cellRatio)

            ggsave_fun(filename = paste0(dir, "/3.CellAnnoted.cellType.group.ratio.plot"), plot = cellRatio, width = 7, height = 7)

            data_p <- scRNA@meta.data[, c("group", "cellType_1")]

            p <- ggstatsplot::ggbarstats(
              data = data_p,
              x = group,
              y = cellType_1,
              palette = "Set3"
            ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            print(p)
            leg <- length(unique(data_p$cellType_1))
            ggsave_fun(filename = paste0(dir, "/4.Dysfunction"), plot = p, width = 2+1.5*leg, height = 7)

            EachSampleClusterDis <- lapply(unique(scRNA@meta.data$orig.ident), function(Sample) {
              seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$orig.ident %in% Sample, ]))
              seur@meta.data$cellType_1 %>%
                table() %>%
                data.frame() %>%
                magrittr::set_colnames(c("CellTypes", "Number")) %>%
                dplyr::mutate(Per = 100 * Number / sum(Number))
            })
            names(EachSampleClusterDis) <- unique(scRNA@meta.data$orig.ident)
            EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
            p1 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + 
              geom_col(position = "fill", width = 0.8) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Percentage of cell") +
              scale_fill_manual(values = my_pal2, name = "") +
              theme(
                axis.text = element_text(color = "black"),
                panel.background = element_blank(), 
                legend.title = element_blank(),
                axis.text.y = element_text(color = "black"),
                axis.line = element_line(color = "black"),
                axis.title.y = element_blank(), 
                legend.position = "right",
                
                
                legend.direction = "vertical", 
                legend.text = element_text(size = 10)
              ) +
              coord_flip()

            EachSampleClusterDis <- lapply(unique(scRNA@meta.data$group), function(Sample) {
              seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$group %in% Sample, ]))
              seur@meta.data$cellType_1 %>%
                table() %>%
                data.frame() %>%
                magrittr::set_colnames(c("CellTypes", "Number")) %>%
                dplyr::mutate(Per = 100 * Number / sum(Number))
            })

            my_pal1 <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#6387C5")

            names(EachSampleClusterDis) <- unique(scRNA@meta.data$group)
            EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
            p2 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + 
              geom_col(position = "fill", width = 0.8) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Percentage of cell") +
              scale_fill_manual(values = my_pal1, name = "") +
              theme(
                axis.text = element_text(color = "black", ),
                panel.background = element_blank(), panel.grid = element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_blank(), axis.title.x = element_blank(), 
                axis.ticks.y.left = element_blank(),
                axis.title.y = element_blank(), 
                axis.line.x = element_line(color = "black"),
                legend.position = "right",
                legend.direction = "vertical", 
                
                
                legend.text = element_text(size = 10)
              ) +
              coord_flip()

            
            cell_counts <- table(scRNA@meta.data$cellType_1)

            
            cell_counts_df <- data.frame(
              CellType = names(cell_counts),
              Count = as.numeric(cell_counts / 1000)
            )

            
            cell_counts_df <- cell_counts_df[order(cell_counts_df$Count, decreasing = T), ]

            p3 <- ggplot(cell_counts_df, aes(x = reorder(CellType, Count), y = Count)) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Number of cell(10^3)") +
              geom_bar(stat = "identity", fill = "#6387C5") +
              theme(
                axis.text = element_text(color = "black"),
                panel.background = element_blank(), panel.grid = element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_blank(),
                
                
                axis.ticks.y.left = element_blank(),
                
                axis.title.y = element_blank(), 
                axis.line.x = element_line(color = "black"),
                legend.text = element_text(size = 10)
              ) +
              coord_flip()
            print(p1 + p2 + p3)
            ggsave_fun(filename = paste0(dir, "/5.MergeRatio"), plot = p1 + p2 + p3, width = 15, height = 6)
            
            pB2_df <- table(scRNA@meta.data$cellType_1, scRNA@meta.data$group) %>% melt()
            colnames(pB2_df) <- c("Cluster", "Sample", "Number")
            pB2_df$Cluster <- factor(pB2_df$Cluster)

            pB2_df <- pB2_df %>%
              group_by(Sample) %>%
              mutate(Percentage = Number / sum(Number))

            pB2_df_tumor <- pB2_df %>%
              filter(Sample %in% c("Disease", "Tumor", "Cancer")) %>%
              select(Cluster, Tumor_Number = Number)

            pB2_df_control <- pB2_df %>%
              filter(Sample == c("Control", "Normal")) %>%
              select(Cluster, Control_Number = Number)

            pB2_df_combined <- merge(pB2_df_tumor, pB2_df_control, by = "Cluster")
            
            pB2_df_combined$rate <- (pB2_df_combined$Tumor_Number - pB2_df_combined$Control_Number) / pB2_df_combined$Control_Number
            KeyCell <- pB2_df_combined[abs(pB2_df_combined$rate) == max(abs(pB2_df_combined$rate)), ]$Cluster %>% as.character()
            print(paste0("变化最大的细胞：", KeyCell))
            write.csv(pB2_df_combined, paste0(dir, "/Cell.rate.csv"))
            
          },
          error = function(e) {
            message("没有分组信息，跳过此步骤\n", e)

            cellRatio <- cellRatioPlot(
              object = scRNA,
              sample.name = "orig.ident",
              celltype.name = "cellType_1",
              flow.curve = 0.5,
              fill.col = my_pal2
            ) + theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1
              )
            )
            print(cellRatio)

            ggsave_fun(filename = paste0(dir, "/3.Sample.ratio.plot"), plot = cellRatio, width = 7 + 0.25 * (length(unique(scRNA$orig.ident))), height = 7)
          }
        )
      }
    }

    
    if (step == 7) {
      message("*************************关键基因表达丰度*************************")

      dir <- paste0(Path, "/.KeyGene expression abundance")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }

      
      KeyGene <- intersect(rownames(scRNA), KeyGene)
      leg <- length(unique(KeyGene))
      
      p <- FeatureCornerAxes(
        object = scRNA, reduction = "umap",
        groupFacet = NULL,
        relLength = 0.5,
        relDist = 0.2,
        features = unique(KeyGene), nLayout = 3,
        aspect.ratio = 1, pSize = 0.1
      )
      if (leg <= 3) {
        ggsave_fun(filename = paste0(dir, "/1.Keygene.FeaturePlot"), plot = p, width = leg * 5, height = 5)
      }
      if (leg > 3) {
        ggsave_fun(filename = paste0(dir, "/1.Keygene.FeaturePlot"), plot = p, width = 15, height = ceiling(leg / 3) * 5)
      }

      p <- DotPlot(scRNA, group.by = "cellType_1", features = unique(KeyGene), cols = "RdYlBu", ) + scale_size_continuous(range = c(0, 10)) + theme(
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        legend.position = "right",
        text = element_text(size = 10)
      )
      if (leg <= 3) {
        ggsave_fun(filename = paste0(dir, "/2.Keygene.DotPlot"), plot = p, width = 5 + ceiling(leg / 3), height = 7)
      }
      if (leg > 3) {
        ggsave_fun(filename = paste0(dir, "/2.Keygene.DotPlot"), plot = p, width = 5 + 0.5 * leg, height = 7)
      }
      
      
      p <- VlnPlot(object = scRNA, features = unique(KeyGene))
      if (leg <= 3) {
        ggsave_fun(filename = paste0(dir, "/3.Keygene.Violin"), plot = p, width = leg * 5, height = 5)
      }
      if (leg > 3) {
        ggsave_fun(filename = paste0(dir, "/3.Keygene.Violin"), plot = p, width = 15, height = ceiling(leg / 3) * 5)
      }
    }

    

  }
  return(scRNA)
}

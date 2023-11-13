.libPaths("/home/weiyihu/package/share/4.0")
setwd("/data1/weiyihu/endometrium/data")
library(Seurat, lib.loc="/home/weiyihu/package/share/4.0")
library(dplyr, lib.loc="/home/weiyihu/package/share/4.0")
library(patchwork, lib.loc="/home/weiyihu/package/share/4.0")
library(Matrix, lib.loc="/home/weiyihu/package/share/4.0")
library(future)
library(future.apply)
library(plotly)

MAXCORES = 4
markerGene = list(T=c("CD3D", "CD8A", "CD4"),
                  Tref=c("FOXP3"),
                  NKT=c("FGFBP2"),
                  B=c("CD79A","CD20"),
                  NK=c("NKG7", "KLRC1", "KLRD1", "NCAM1"),
                  Neutrophils=c("CST3", "LYZ", "FCGR3B", "CSF3R"),
                  Monocytes=c("FCN1", "VCAN", "CD14", "CD16", "S100A12", "FCGR3A"),
                  Macrophage=c("CD14", "FOLR2", "CD163"),
                  DC=c("LILRA4", "CLEC9A", "CD1C"),
                  Mast=c("TPSAB1", "MS4A2", "KIT"),
                  ILC3=c("NCR2")
                  )
markerGeneSub = list(T=list(CTL=toupper(c("Pdcd1", "Nkg7", "Klrc1", "Klrd1", "Klrk1")),
                            CD4CD8=toupper(c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1")),
                            CD8=toupper(c("Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd8b1", "Gzma", "Gzmk", "Gzmm", "Prf1", "Tnfsf10")),
                            CD3=toupper(c("Cd3d", "Cd3e", "Cd3g")),
                            CD4=toupper(c("Cd3d", "Cd3e", "Cd3g", "Cd4")),
                            Tfh=toupper(c("Tcf7", "Sh2d1a", "Cd40lg", "Cd69", "Slamf6")),
                            NK=toupper(c("Klrb1c", "Klra3", "Fcer1g", "Gzma", "Gzmb", "Klra8 ,Klre1", "Klra7 ", "Klrk1", "Klrd1", "Klra4", "Serpinb6b", "Il2rb", "Ctsw", "Nkg7", "Xcl1", "Ccl5", "Ms4a4b", "Id2", "Ccl4")),
                            Th17=toupper(c("Il17a", "Il17f", "Il23r", "Blk", "Tcf12", "Gzma", "Gzmb", "Prf1", "Ccl3", "Ccl4", "Ccl5", "Ccl20")),
                            Th1=toupper(c("Cxcr3", "Ccl5", "Ms4a4b", "Ifng", "Tbx21", "Gzmk")),
                            Th2=toupper(c("Gata3", "Il13", "Il1rl1", "Plac8", "Igfbp7")),
                            Th22=toupper(c("Il13", "Il22", "Ccr4", "Cxcr6")),
                            Treg=toupper(c("Foxp3", "Ikzf2", "Il2ra", "Ctla4"))
                            ),
                     B=list(Memory_B=toupper(c("Ms4a1", "Cd27", "Cd40", "Cd80", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6")),
                            DZ_B=toupper(c("Birc5", "Ube2c", "Cdc20", "Stmn1")),
                            LZ_B=toupper(c("Aicda", "S1pr2", "Lipc", "Rgs1")),
                            Plasma_B=toupper(c("Ly6c1", "Ly6c2", "Xbp1", "Jchain")),
                            AntigenPresenting_B=toupper(c("Rsc22d3", "Ccr7", "Klf2", "Cd83", "Cd52", "Serpinb1a", "S1pr4", "Ms4a4c", "Plac8", "H2-K1", "H2-Q7", "Ly6a")),
                            IGSpecific_B=toupper(c("Pax5", "Rel", "Mcl1", "Mycbp2")))
                     )


# 一级函数
computeMeanExpressionInCluster = function(SeuratObject, geneSet, cluster, debug=FALSE){
  # 查询SeuratObject中指定geneSet的平均表达值
  SeuratObject = SeuratObject
  geneSet = geneSet
  cluster = cluster
  debug = debug
  
  info = t(SeuratObject@assays$integrated@scale.data)
  info = as.data.frame(info)
  geneSet = geneSet[geneSet %in% colnames(info)]
  if(length(geneSet)==0){return(-10000)}
  info = info[, geneSet, drop=FALSE]
  info = info[SeuratObject@active.ident==cluster, , drop=FALSE]
  info = colMeans(info)
  if(debug == FALSE){info = mean(info)}
  return(info)
}

# 一级函数
computeMarkersDifRes = function(seuratObject,
                                resolutionVector=seq(0.1, 2, 0.1)){
  # input:
  #   seuratObject|\t seuratObject|\t Seurat4
  #   resolutionVector|\t vector|\t =seq(0.1, 2, 0.1)
  # change:
  #   seuratObject在不同resolution下的FindAllMarkers()
  # output:
  #   markersDifRes|\t list|\t 一个list, key为resolution, value为FindAllMarkers()的结果
  seuratObject = seuratObject
  resolutionVector = resolutionVector
  
  markersDifRes = list()
  for(r in resolutionVector){
    r = as.character(r)
    Idents(seuratObject) = paste0("integrated_snn_res.", r)
    markersDifRes[[r]] = FindAllMarkers(seuratObject, only.pos=TRUE)
  }
  
  return(markersDifRes)
}

# 一级函数
computeCellTypeInResolution = function(markersDifRes, markerGene,
                                       resolutionVector=seq(0.1, 2, 0.1)){
  # input:
  #   markersDifRes|\t list|\t 一个list, key为resolution, value为在该resolution下的FindAllMarkers()结果
  #   markerGene|\t list|\t 一个list, key为细胞类型, value为该细胞类型的标记基因
  #   resolutionVector|\t vector|\t =seq(0.1, 2, 0.3)
  # change:
  #   获取scRNAMerged在不同resolution下的top10差异基因
  # output:
  #   markers|\t data.frame|\t 含有在不同resolution下, 不同cluster所表达的细胞类型标记基因
  markersDifRes = markersDifRes
  markerGene = markerGene
  resolutionVector = resolutionVector
  
  markers = list()
  for(resolution in resolutionVector){
    resolution = as.character(resolution)
    markers[[resolution]] = list(resolution=resolution)
    tempMarkers = markersDifRes[[resolution]]
    top10 = dplyr::top_n(dplyr::group_by(tempMarkers, cluster), n=10, wt=avg_log2FC)
    
    info = data.frame(row.names=c("resolution", "cluster", "cellType", "gene", "avg_log2FC", "pct.1", "pct.2"))
    for(i in seq(1, dim(top10)[1])){
      cluster = as.character(top10[[i, "cluster"]])
      gene = as.character(top10[[i, "gene"]])
      avg_log2FC = as.character(top10[[i, "avg_log2FC"]])
      pct.1 = as.character(top10[[i, "pct.1"]])
      pct.2 = as.character(top10[[i, "pct.2"]])
      for(cellType in names(markerGene)){if(gene %in% markerGene[[cellType]]){info[as.character(dim(info)[2]+1)] = c(resolution, cluster, cellType, gene, avg_log2FC, pct.1, pct.2)}}
    }
    markers[[as.character(resolution)]][["df"]] = info
  }
  dfTemp = data.frame(row.names=c("resolution", "cluster", "cellType", "gene", "avg_log2FC", "pct.1", "pct.2"))
  for(r in markers){
    df = r[["df"]]
    dfTemp = cbind(dfTemp, df)
  }
  dfTemp = t(dfTemp)
  markers[["df"]] = dfTemp
  
  return(markers[["df"]])
}

# 一级函数
makeListVlnPlot = function(seuratObject, markerGene){
  # input:
  #   seuratObject|\t seuratObject|\t Seurat4对象
  #   markerGene|\t list|\t 一个list, key为细胞类型, value为该细胞类型的标记基因
  # change:
  #   查看在不同resolution下, 标记基因在不同cluster的表达情况
  # output:
  #   listVlnPlot|\t list|\t 一个list, key为resolution, value为plot
  seuratObject = seuratObject
  markerGene = markerGene
  
  # 汇总在seuratObject中出现的所有标记基因
  gene = c()
  for(g in markerGene){
    g = g[g %in% toupper(rownames(seuratObject))]
    gene = append(gene, g)
  }
  # 绘制stackedVlnPlot
  listVlnPlot = list()
  for(i in seq(0.1, 2, 0.1)){
    i = as.character(i)
    print(i)
    Idents(seuratObject) = paste0("integrated_snn_res.", i)
    listVlnPlot[[i]] = VlnPlot(seuratObject, features=gene, pt.size=0, stack=T) + NoLegend()
  }
  
  return(listVlnPlot)
}

# 二级函数
getCellType = function(SeuratObject, cluster, markerGeneList, debug=FALSE){
  # 判断亚类的细胞类型
  SeuratObject = SeuratObject
  cluster = cluster
  markerGeneList = markerGeneList
  debug = debug
  
  result = c()
  progress = 0
  progressLen = length(markerGeneList)
  for(n in names(markerGeneList)){
    progress = progress + 1
    cat(progress, '/', progressLen, '\r')
    resultTemp = computeMeanExpressionInCluster(SeuratObject, geneSet=markerGeneList[[n]], cluster=cluster)
    result = c(result, resultTemp)
  }
  names(result) = names(markerGeneList)
  if(debug == FALSE){result = names(which.max(result))}
  return(result)
}


# 前置参数
if(TRUE){
  plan(multisession, workers=MAXCORES)
  options(future.globals.maxSize=10*1024*1024*1024,  # 10GB
          scipen=10)  
}

# 加载数据
if(TRUE){
  scRNAMerged = readRDS(file.path("share", "All_data_Merge.rds"))
  
}

# 合并前的处理
if(TRUE){
  # 分离样本数据
  scRNAList = SplitObject(scRNAMerged, split.by="orig.ident")
  saveRDS(scRNAList, "0003_scRNAList.rds")
  # 质控
  scRNAMerged = future_lapply(scRNAList, future.seed=TRUE, function(seuratObject){
    seuratObject <- subset(seuratObject, subset=nCount_RNA<100000 & nCount_RNA>1000)
  })
  # normalize and identify variable features for each dataset independently
  scRNAMerged = future_lapply(scRNAMerged, future.seed=TRUE, function(seuratObject){
    seuratObject = NormalizeData(seuratObject)
    seuratObject = FindVariableFeatures(seuratObject, selection.method="vst", nfeatures=2000)
    return(seuratObject)
  })
  # 保存数据
  saveRDS(scRNAMerged, "0003_scRNAMergedRaw.rds")
}



# harmony合并
if(FALSE){
  # 注释掉的部分在VScode中执行
  # scRNAMerged = readRDS("0003_scRNAMergedRaw.rds")
  # scRNAMerged = merge(scRNAMerged[[1]], y=c(scRNAMerged[-1]))
  # scRNAMerged <- NormalizeData(scRNAMerged) %>% 
  #                   FindVariableFeatures() %>% 
  #                   ScaleData() %>% 
  #                   RunPCA(verbose=FALSE)
  # scRNAMerged <- RunHarmony(scRNAMerged, group.by.vars="orig.ident", plot_convergence=TRUE)
  # saveRDS(scRNAMerged, "0003_scRNAMergedHarmony.rds")
}

# 准备确定resolution
if(TRUE){
  scRNAMerged = readRDS("0003_scRNAMergedHarmony.rds")
  ElbowPlot(scRNAMerged, ndims=25)
  pc.dim = 1:25
  scRNAMerged <- RunUMAP(scRNAMerged, reduction="harmony", dims=pc.dim)
  scRNAMerged <- FindNeighbors(scRNAMerged, reduction="harmony", dims=pc.dim)
  scRNAMerged = FindClusters(scRNAMerged, resolution=seq(0.1, 2, 0.1))
  DefaultAssay(scRNAMerged) <- "integrated"
  plot1 = DimPlot(scRNAMerged, reduction="umap", label=T) 
  plot2 = DimPlot(scRNAMerged, reduction="umap", group.by='orig.ident')
  pdf("0003_harmonyUMAP.pdf")
  print(plot1)
  print(plot2)
  dev.off()
  saveRDS(scRNAMerged, "0003_scRNAMergedClusters.rds")
}

# 确定resolution
if(FALSE){
  #scRNAMerged = readRDS("0003_scRNAMergedClusters.rds")
  #plot = clustree(scRNAMerged, prefix="integrated_snn_res.")
  #pdf("0003_clustree.pdf", width=14, height=14)
  #print(plot)
  #dev.off()
  resolution = 0.5 %>% as.character()
  
  Idents(scRNAMerged) = paste0("integrated_snn_res.", resolution)
  pdf("0003_UMAPOfClusters.pdf")
  DimPlot(scRNAMerged, reduction="umap", group.by=paste0("integrated_snn_res.", resolution), label=TRUE)
  dev.off()
}

# 确定细胞类型
if(TRUE){
  # 准备工作
  if(TRUE){
    # 确定不同resolution下的cluster间的差异基因
    markersDifRes = computeMarkersDifRes(seuratObject=scRNAMerged, resolutionVector=seq(0.1, 2, 0.3))
    # 确定不同resolution下cluster的差异基因中的标记基因
    dftemp = computeCellTypeInResolution(markersDifRes=markersDifRes,
                                         markerGene=markerGene,
                                         resolutionVector=seq(0.1, 2, 0.1))
    # 查看在不同resolution下, 标记基因在不同cluster的表达情况
    listVlnPlot = makeListVlnPlot(seuratObject=scRNAMerged, markerGene=markerGene)
  }
  # 判断细胞类型
  if(TRUE){
    # 结合clustree确定resolution值
    resolution = "0.5"
    Idents(scRNAMerged) = paste0("integrated_snn_res.", resolution)
    # 结合clustree判断cluster的可靠细胞类型, 及可疑细胞类型
    dfTemp = dfTemp
    # 结合标记基因在cluster中的平均表达量判断cluster的可疑细胞类型
    temp = list()
    progressLen = length(levels(Idents(scRNAMerged)))
    for(i in levels(Idents(scRNAMerged))){
      cat(i, '/', progressLen, '\n')
      temp[[as.character(i)]] = getCellType(scRNAMerged, cluster=as.numeric(i), markerGeneList=markerGene, debug=TRUE)
    }
    temp = do.call(data.frame, temp)
    temp = t(temp)
    rownames(temp) = gsub("X", "", rownames(temp))
    # 结合VlnPlot验证确定的细胞类型
    listVlnPlot[[resolution]] %>% print()
  }
  # 细胞类型注释
  if(TRUE){
    new.cluster.ids <- c("", "", "", "", "", "", "")
    names(new.cluster.ids) <- levels(scRNAMerged)
    scRNAMerged <- RenameIdents(scRNAMerged, new.cluster.ids)
    scRNAMerged@meta.data$cellType = Idents(scRNAMerged)
  }
  # 查看细胞类型注释的质量
  if(TRUE){
    # UMAP图 在指定的resolution下
    print(DimPlot(scRNAMerged, reduction="umap", group.by=paste0("integrated_snn_res.", resolution), pt.size=0.5, label=TRUE))
    # UMAP图 在分类的细胞类型下
    print(DimPlot(scRNAMerged, reduction="umap", pt.size=0.5, label=TRUE))
    # 热图 在cluster的top10下
    temp = subset(scRNAMerged, downsample=100)
    top10 = dplyr::top_n(dplyr::group_by(markersDifRes[[resolution]], cluster), n=10, wt=avg_log2FC)
    top10 = as.data.frame(top10)
    top10 = top10[, c("cluster", "gene", "avg_log2FC")]
    print(DoHeatmap(temp, features=top10$gene)+ NoLegend())
    # 热图 在cellType的标记基因下
    gene = unlist(markerGene) %>% unique()
    gene = gene[gene %in% rownames(scRNAMerged)]
    print(DoHeatmap(temp, features=gene)+ NoLegend())
    # 热图 在cluster的top10 + cellType的标记基因下
    levels(top10$cluster) = new.cluster.ids
    for(cellType in names(markerGene)){
      cellType2 = gsub("_.*", "", cellType)
      if(cellType2 %in% levels(top10$cluster)){
        gene = markerGene[[cellType]]
        gene = gene[gene %in% rownames(scRNAMerged)]
        gene = gene[!gene %in% top10$gene]
        for(i in gene){
          top10 = rbind(top10, c(cellType2, i, 0))
        }
      }
    }
    top10 = top10[order(match(top10$cluster, levels(Idents(scRNAMerged))), -as.numeric(top10$avg_log2FC)),]
    print(DoHeatmap(temp, features=top10$gene)+ NoLegend())
  }
}



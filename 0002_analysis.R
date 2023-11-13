.libPaths("/home/weiyihu/package/share/4.0")
setwd("/data1/weiyihu/endometrium/data")
library(Seurat, lib.loc="/home/weiyihu/package/share/4.0")
library(dplyr, lib.loc="/home/weiyihu/package/share/4.0")
library(patchwork, lib.loc="/home/weiyihu/package/share/4.0")
library(Matrix, lib.loc="/home/weiyihu/package/share/4.0")
library(future)
library(future.apply)
library(plotly)


MAXCORES = 6
markerGene = list(T = list(CTL = toupper(c("Pdcd1", "Nkg7", "Klrc1", "Klrd1", "Klrk1")),
                           CD4CD8 = toupper(c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1")),
                           CD8 = toupper(c("Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd8b1", "Gzma", "Gzmk", "Gzmm", "Prf1", "Tnfsf10")),
                           CD3 = toupper(c("Cd3d", "Cd3e", "Cd3g")),
                           CD4 = toupper(c("Cd3d", "Cd3e", "Cd3g", "Cd4")),
                           Tfh = toupper(c("Tcf7", "Sh2d1a", "Cd40lg", "Cd69", "Slamf6")),
                           NK = toupper(c("Klrb1c", "Klra3", "Fcer1g", "Gzma", "Gzmb", "Klra8 ,Klre1", "Klra7 ", "Klrk1", "Klrd1", "Klra4", "Serpinb6b", "Il2rb", "Ctsw", "Nkg7", "Xcl1", "Ccl5", "Ms4a4b", "Id2", "Ccl4")),
                           Th17 = toupper(c("Il17a", "Il17f", "Il23r", "Blk", "Tcf12", "Gzma", "Gzmb", "Prf1", "Ccl3", "Ccl4", "Ccl5", "Ccl20")),
                           Th1 = toupper(c("Cxcr3", "Ccl5", "Ms4a4b", "Ifng", "Tbx21", "Gzmk")),
                           Th2 = toupper(c("Gata3", "Il13", "Il1rl1", "Plac8", "Igfbp7")),
                           Th22 = toupper(c("Il13", "Il22", "Ccr4", "Cxcr6")),
                           Treg = toupper(c("Foxp3", "Ikzf2", "Il2ra", "Ctla4"))),
                  B = list(Memory_B=toupper(c("Ms4a1", "Cd27", "Cd40", "Cd80", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6")),
                           DZ_B=toupper(c("Birc5", "Ube2c", "Cdc20", "Stmn1")),
                           LZ_B=toupper(c("Aicda", "S1pr2", "Lipc", "Rgs1")),
                           Plasma_B=toupper(c("Ly6c1", "Ly6c2", "Xbp1", "Jchain")),
                           AntigenPresenting_B=toupper(c("Rsc22d3", "Ccr7", "Klf2", "Cd83", "Cd52", "Serpinb1a", "S1pr4", "Ms4a4c", "Plac8", "H2-K1", "H2-Q7", "Ly6a")),
                           IGSpecific_B=toupper(c("Pax5", "Rel", "Mcl1", "Mycbp2")))
                  )

# 一级函数
if(TRUE){
  ## remove the x-axis text and tick
  ## plot.margin to adjust the white space between each plot.
  ## ... pass any arguments to VlnPlot in Seurat
  modify_vlnplot<- function(obj, 
                            feature, 
                            pt.size = 0, 
                            plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                            ...) {
    p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
      xlab("") + ylab(feature) + ggtitle("") + 
      theme(legend.position = "none", 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.title.y = element_text(size = rel(1), angle = 0), 
            axis.text.y = element_text(size = rel(1)), 
            plot.margin = plot.margin ) 
    return(p)
  }
  
  ## extract the max value of the y axis
  extract_max<- function(p){
    ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
  }
  
  
  ## main function
  StackedVlnPlot<- function(obj, features,
                            pt.size = 0, 
                            plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                            ...) {
    
    plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
    
    # Add back x-axis title to bottom plot. patchwork is going to support this?
    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
      theme(axis.text.x=element_text(), axis.ticks.x = element_line())
    
    # change the y-axis tick to only max value 
    ymaxs<- purrr::map_dbl(plot_list, extract_max)
    plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                              scale_y_continuous(breaks = c(y)) + 
                              expand_limits(y = y))
    
    p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
    return(p)
  }
}

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
  scRNAMerged = readRDS(file.path("rds", "0001_scRNAMergedAfterAnnotation.rds"))
}


# main
# 对B细胞进行亚分类
if(TRUE){
  # 重新生成SeuratObject
  if(TRUE){
    data = scRNAMerged[, Idents(scRNAMerged) %in% c("B")]
    data = SplitObject(data, split.by="orig.ident")
    data = lapply(data, function(SeuratObject){SeuratObject@assays$RNA@counts})
    data = lapply(data, function(cm){CreateSeuratObject(cm)})
    for(n in names(data)){
      data[[n]]@meta.data[["sample"]] = n
    }
  }
  # 按study合并seuratObject
  if(TRUE){
    data2 = list(GSE164449=merge(data$GSE164449_Ctrl1, y=data$GSE164449_RPL1),
                 CRA002181=merge(data$CRA002181_Ctrl1, y=list(data$CRA002181_Ctrl2, data$CRA002181_Ctrl3,
                                                              data$CRA002181_RPL1, data$CRA002181_RPL2)),
                 HRA000237=merge(data$HRA000237_Ctrl1, y=list(data$HRA000237_Ctrl2, data$HRA000237_Ctrl3,
                                                              data$HRA000237_RPL1, data$HRA000237_RPL2,
                                                              data$HRA000237_RPL3)),
                 GSE130560=merge(data$GSE130560_Ctrl1, y=list(data$GSE130560_Ctrl2, data$GSE130560_Ctrl3,
                                                              data$GSE130560_RSA1, data$GSE130560_RSA2,
                                                              data$GSE130560_RSA3)))
  }
  # 质控
  if(TRUE){
    data2 = lapply(data2, function(seuratObject){
      seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern="^MT-")
      seuratObject[["percent.rb"]] <- PercentageFeatureSet(seuratObject, pattern="^RP[SL]")
      return(seuratObject)
    })
    data2 = lapply(data2, function(seuratObject){
      seuratObject <- subset(seuratObject, subset=nFeature_RNA>500 & nFeature_RNA<3000 & percent.mt<20)
      return(seuratObject)
    })
  }
  # 合并
  if(TRUE){
    data2 = lapply(data2, function(sample){
      sample = NormalizeData(sample)
      sample = FindVariableFeatures(sample, selection.method="vst", nfeatures=2000)
    })
    features = SelectIntegrationFeatures(object.list=data2)
    anchors = FindIntegrationAnchors(object.list=data2, anchor.features=features)
    data2 = IntegrateData(anchorset=anchors)
    DefaultAssay(data2) = "integrated"
  }
  if(TRUE){
    # 添加disease信息
    data2@meta.data[grepl("Ctrl", data2@meta.data$sample), "disease"] = "Healthy"
    data2@meta.data[grepl("RPL", data2@meta.data$sample), "disease"] = "RPL"
    data2@meta.data[grepl("RSA", data2@meta.data$sample), "disease"] = "RSA"
  }
  # 对合并后数据进行标准流程处理
  if(TRUE){
    data2 = ScaleData(data2)
    data2 = RunPCA(data2)
    ElbowPlot(data2, ndims=25)
    pc.dim = 1:5
    data2 <- RunUMAP(data2, reduction="pca", dims=pc.dim)
    data2 <- RunTSNE(data2, reduction="pca", dims=pc.dim)
    data2 <- FindNeighbors(data2, reduction="pca", dims=pc.dim)
    data2 <- FindClusters(data2, resolution=seq(0.1, 2, 0.1))
    gene = unlist(markerGene$B) %>% as.character()
    gene = gene[gene %in% rownames(data2)]
  }
  # 寻找不同resolution下的cluster间的差异基因
  if(TRUE){
    markersDifRes = list()
    for(r in seq(0.1, 2, 0.1)){
      r = as.character(r)
      Idents(data2) = paste0("integrated_snn_res.", r)
      markersDifRes[[r]] = FindAllMarkers(data2, only.pos=TRUE)
    }
  }
  # 寻找不同resolution下cluster的差异基因中的标记基因
  if(TRUE){
    markersB = list()
    for(resolution in seq(0.1, 2, 0.1)){
      resolution = as.character(resolution)
      markersB[[resolution]] = list(resolution=resolution)
      tempMarkers = markersDifRes[[resolution]]
      top10 = dplyr::top_n(dplyr::group_by(tempMarkers, cluster), n=10, wt=avg_log2FC)
      
      info = data.frame(row.names=c("resolution", "cluster", "cellType", "gene", "avg_log2FC", "pct.1", "pct.2"))
      for(i in seq(1, dim(top10)[1])){
        cluster = as.character(top10[[i, "cluster"]])
        gene = as.character(top10[[i, "gene"]])
        avg_log2FC = as.character(top10[[i, "avg_log2FC"]])
        pct.1 = as.character(top10[[i, "pct.1"]])
        pct.2 = as.character(top10[[i, "pct.2"]])
        for(cellType in names(markerGene[["B"]])){
          if(gene %in% markerGene[["B"]][[cellType]]){info[as.character(dim(info)[2]+1)] = c(resolution, cluster, cellType, gene, avg_log2FC, pct.1, pct.2)}
        }
      }
      markersB[[as.character(resolution)]][["df"]] = info
    }
    dfTemp = data.frame(row.names=c("resolution", "cluster", "cellType", "gene", "avg_log2FC", "pct.1", "pct.2"))
    for(r in markersB){
      df = r[["df"]]
      dfTemp = cbind(dfTemp, df)
    }
    dfTemp = t(dfTemp)
    markersB[["df"]] = dfTemp
  }
  # 统计在不同resolution下, 标记基因在不同cluster的表达情况
  if(TRUE){
    gene = c()
    for(g in markerGene[["B"]]){
      g = g[g %in% toupper(rownames(data2))]
      gene = append(gene, g)
    }
    p = list()
    for(i in seq(0.1, 2, 0.1)){
      i = as.character(i)
      print(i)
      Idents(data2) = paste0("integrated_snn_res.", i)
      #plot = StackedVlnPlot(temp, gene)
      plot = VlnPlot(data2, features=gene, pt.size=0, stack=T)+ NoLegend()
      p[[i]] = plot
    }
  }
  # 判断细胞类型
  if(TRUE){
    # 结合clustree确定resolution值
    resolution = "0.4"
    Idents(data2) = paste0("integrated_snn_res.", resolution)
    # 结合clustree判断cluster的可靠细胞类型, 及可疑细胞类型
    dfTemp = markersB[["df"]]
    # 结合标记基因在cluster中的平均表达量判断cluster的可疑细胞类型
    temp = list()
    progressLen = length(levels(Idents(data2)))
    for(i in levels(Idents(data2))){
      cat(i, '/', progressLen, '\n')
      temp[[as.character(i)]] = getCellType(data2, cluster=as.numeric(i), markerGeneList=markerGene[["B"]], debug=TRUE)
    }
    temp = do.call(data.frame, temp)
    temp = t(temp)
    rownames(temp) = gsub("X", "", rownames(temp))
    # 结合VlnPlot验证确定的细胞类型
    p[[resolution]] %>% print()
  }
  # 细胞类型注释
  if(TRUE){
    new.cluster.ids <- c("Antigen", "DZ", "Memory", "DZ", "LZ", "LZ", "Plasma")
    names(new.cluster.ids) <- levels(data2)
    data2 <- RenameIdents(data2, new.cluster.ids)
  }
  # 查看细胞类型注释的质量
  if(TRUE){
    # UMAP图 在指定的resolution下
    print(DimPlot(data2, reduction="umap", group.by=paste0("integrated_snn_res.", resolution), pt.size=0.5, label=TRUE))
    # UMAP图 在分类的细胞类型下
    print(DimPlot(data2, reduction="umap", pt.size=0.5, label=TRUE))
    # 热图 在cluster的top10下
    temp = subset(data2, downsample=100)
    top10 = dplyr::top_n(dplyr::group_by(markersDifRes[[resolution]], cluster), n=10, wt=avg_log2FC)
    top10 = as.data.frame(top10)
    top10 = top10[, c("cluster", "gene", "avg_log2FC")]
    print(DoHeatmap(temp, features=top10$gene)+ NoLegend())
    # 热图 在cellType的标记基因下
    gene = unlist(markerGene[["B"]]) %>% unique()
    gene = gene[gene %in% rownames(data2)]
    print(DoHeatmap(temp, features=gene)+ NoLegend())
    # 热图 在cluster的top10 + cellType的标记基因下
    levels(top10$cluster) = new.cluster.ids
    for(cellType in names(markerGene[["B"]])){
      cellType2 = gsub("_.*", "", cellType)
      if(cellType2 %in% levels(top10$cluster)){
        gene = markerGene[["B"]][[cellType]]
        gene = gene[gene %in% rownames(data2)]
        gene = gene[!gene %in% top10$gene]
        for(i in gene){
          top10 = rbind(top10, c(cellType2, i, 0))
        }
      }
    }
    top10 = top10[order(match(top10$cluster, levels(Idents(data2))), -as.numeric(top10$avg_log2FC)),]
    print(DoHeatmap(temp, features=top10$gene)+ NoLegend())
  }
  # 将B细胞亚型注释添加到scRNAMerged中
  if(TRUE){
    debug = scRNAMerged@active.ident %>% as.data.frame()
    colnames(debug) = c("cellType")
    debug["cellType"] = as.character(debug[["cellType"]])
    debug["barcode"] = rownames(debug)
    debug2 = data2@active.ident %>% as.data.frame()
    colnames(debug2) = c("cellType")
    debug2["barcode"] = rownames(debug2)
    debug2[["cellType"]] = as.character(debug2[["cellType"]]) %>% paste0(., "_B")
    debug3 = match(debug2[["barcode"]], debug[["barcode"]])
    debug[debug3, "cellType"] = debug2[["cellType"]]
    scRNAMerged@meta.data$cellType = debug$cellType
    Idents(scRNAMerged) = "cellType"
    DimPlot(scRNAMerged, reduction="umap", pt.size=0.5)
  }
}

# 统计数据
if(TRUE){
  # 统计在不同disease中不同类型细胞的占比
  if(TRUE){
    df <- data.frame(barcode=names(temp@active.ident),
                     cellType=temp@active.ident,
                     disease=temp@meta.data$disease)
    df = table(df$cellType, df$disease) %>% addmargins() %>% as.data.frame.array()
    for(c in colnames(df)){
      df[c] = df[c] / df["Sum", c]
    }
    df = df[-which(rownames(df)=="Sum"), -which(colnames(df)=="Sum")]
    df = t(df) %>% as.data.frame()
    df[["disease"]] = rownames(df)
  }
  # 根据统计数据绘制柱状图
  if(TRUE){
    fig <- plot_ly(df)
    for(c in colnames(df)[-which(colnames(df)=="disease")]){
      fig = add_trace(fig, x=df[["disease"]], y=df[[c]], type="bar", name=c)
    }
    fig = layout(fig, yaxis=list(title='Count'), barmode='stack')
  }
}




if(TRUE){
  # 按sample合并
  data3 = lapply(data, function(seuratObject){
    seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern="^MT-")
    seuratObject[["percent.rb"]] <- PercentageFeatureSet(seuratObject, pattern="^RP[SL]")
    return(seuratObject)
  })
  data3 = lapply(data3, function(seuratObject){
    seuratObject <- subset(seuratObject, subset=nFeature_RNA>500 & nFeature_RNA<3000 & percent.mt<20)
    return(seuratObject)
  })
  data3 = lapply(data3, function(seuratObject){
    seuratObject = NormalizeData(seuratObject)
    seuratObject = FindVariableFeatures(seuratObject, selection.method="vst", nfeatures=2000)
  })
  features = SelectIntegrationFeatures(object.list=data3)
  anchors = FindIntegrationAnchors(object.list=data3, anchor.features=features, dims=1:9)
  data3 = IntegrateData(anchorset=anchors)
  DefaultAssay(data3) = "integrated"
  data3 = ScaleData(data3, verbose=FALSE)
  data3 = RunPCA(data3, verbose=FALSE)
  ElbowPlot(data3, ndims=25)
  pc.dim = 1:5
  data3 <- RunUMAP(data3, reduction="pca", dims=pc.dim)
  data3 <- FindNeighbors(data3, reduction="pca", dims=pc.dim)
  data3 <- FindClusters(data3, resolution=seq(0.1, 2, 0.1))
}

# 手动比较
if(TRUE){
  Idents(data2) = "integrated_snn_res.0.4"
  debug = markersDifRes[["0.4"]]
  debug = debug[markersDifRes[["0.4"]]$cluster=="0",]
  debug = debug[debug$gene %in% gene,]
  p[["0.4"]]
  DimPlot(data2, reduction="umap", group.by="integrated_snn_res.0.4", pt.size=0.5, label=TRUE)
  FeaturePlot(data2, reduction="umap", features=markerGene[["B"]]$DZ_B, pt.size=0.5, min.cutoff=0, max.cutoff=2)
  getCellType(data2, cluster=1, markerGeneList=markerGene[["B"]], debug=TRUE)
}

# 绘制热图
if(TRUE){
  data2@meta.data$disease2 = gsub("RSA", "disease", data2@meta.data$disease)
  data2@meta.data$disease2 = gsub("RPL", "disease", data2@meta.data$disease2)
  m <- matrix(rnorm(9), nrow=3, ncol=3)
  m = df[, -which(colnames(df)=="disease")]
  m = t(m)
  fig <- plot_ly(
    x = colnames(m),
    y = rownames(m),
    z = m,
    type = "heatmap",
    colors = colorRamp(c("darkblue", "darkred"))
  )
  fig
}



# T细胞亚分类
if(TRUE){
  # 重新生成SeuratObject
  if(TRUE){
    data = scRNAMerged[, Idents(scRNAMerged) %in% c("T")]
    data = SplitObject(data, split.by="orig.ident")
    data = lapply(data, function(SeuratObject){SeuratObject@assays$RNA@counts})
    data = lapply(data, function(cm){CreateSeuratObject(cm)})
    for(n in names(data)){
      data[[n]]@meta.data[["sample"]] = n
    }
  }
  # 质控
  if(TRUE){
    data2 = lapply(data2, function(seuratObject){
      seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern="^MT-")
      seuratObject[["percent.rb"]] <- PercentageFeatureSet(seuratObject, pattern="^RP[SL]")
      return(seuratObject)
    })
    data2 = lapply(data2, function(seuratObject){
      seuratObject <- subset(seuratObject, subset=nFeature_RNA>500 & nFeature_RNA<3000 & percent.mt<20)
      return(seuratObject)
    })
  }
  # 合并
  if(TRUE){
    data2 = lapply(data2, function(sample){
      sample = NormalizeData(sample)
      sample = FindVariableFeatures(sample, selection.method="vst", nfeatures=2000)
    })
    features = SelectIntegrationFeatures(object.list=data2)
    anchors = FindIntegrationAnchors(object.list=data2, anchor.features=features)
    
    data2 = IntegrateData(anchorset=anchors)
    DefaultAssay(data2) = "integrated"
  }
  if(TRUE){
    # 添加disease信息
    data2@meta.data[grepl("Ctrl", data2@meta.data$sample), "disease"] = "Healthy"
    data2@meta.data[grepl("RPL", data2@meta.data$sample), "disease"] = "RPL"
    data2@meta.data[grepl("RSA", data2@meta.data$sample), "disease"] = "RSA"
  }
  # 对合并后数据进行标准流程处理
  if(TRUE){
    data2 = ScaleData(data2)
    data2 = RunPCA(data2)
    ElbowPlot(data2, ndims=25)
    pc.dim = 1:15
    data2 <- RunUMAP(data2, reduction="pca", dims=pc.dim)
    data2 <- RunTSNE(data2, reduction="pca", dims=pc.dim)
    data2 <- FindNeighbors(data2, reduction="pca", dims=pc.dim)
    data2 <- FindClusters(data2, resolution=seq(0.1, 2, 0.1))
    gene = unlist(markerGene$T) %>% as.character()
    gene = gene[gene %in% rownames(data2)]
    saveRDS(data2, "0002_tempTdata2.rds")
  }
  # 寻找不同resolution下的cluster间的差异基因
  if(TRUE){
    markersDifRes = list()
    for(r in seq(0.1, 2, 0.1)){
      r = as.character(r)
      Idents(data2) = paste0("integrated_snn_res.", r)
      markersDifRes[[r]] = FindAllMarkers(data2, only.pos=TRUE)
    }
    saveRDS(markersDifRes, "0002_tempMarkersDifResT.rds")
  }
  # 寻找不同resolution下cluster的差异基因中的标记基因
  if(TRUE){
    markersT = list()
    for(resolution in seq(0.1, 2, 0.1)){
      resolution = as.character(resolution)
      markersT[[resolution]] = list(resolution=resolution)
      tempMarkers = markersDifRes[[resolution]]
      top10 = dplyr::top_n(dplyr::group_by(tempMarkers, cluster), n=10, wt=avg_log2FC)
      
      info = data.frame(row.names=c("resolution", "cluster", "cellType", "gene", "avg_log2FC", "pct.1", "pct.2"))
      for(i in seq(1, dim(top10)[1])){
        cluster = as.character(top10[[i, "cluster"]])
        gene = as.character(top10[[i, "gene"]])
        avg_log2FC = as.character(top10[[i, "avg_log2FC"]])
        pct.1 = as.character(top10[[i, "pct.1"]])
        pct.2 = as.character(top10[[i, "pct.2"]])
        for(cellType in names(markerGene[["T"]])){
          if(gene %in% markerGene[["T"]][[cellType]]){info[as.character(dim(info)[2]+1)] = c(resolution, cluster, cellType, gene, avg_log2FC, pct.1, pct.2)}
        }
      }
      markersT[[as.character(resolution)]][["df"]] = info
    }
    dfTemp = data.frame(row.names=c("resolution", "cluster", "cellType", "gene", "avg_log2FC", "pct.1", "pct.2"))
    for(r in markersT){
      df = r[["df"]]
      dfTemp = cbind(dfTemp, df)
    }
    dfTemp = t(dfTemp)
    markersT[["df"]] = dfTemp
  }
  # 统计在不同resolution下, 标记基因在不同cluster的表达情况
  if(TRUE){
    gene = c()
    for(g in markerGene[["T"]]){
      g = g[g %in% toupper(rownames(data2))]
      gene = append(gene, g)
    }
    gene = unique(gene)
    p = list()
    for(i in seq(0.1, 2, 0.1)){
      i = as.character(i)
      print(i)
      Idents(data2) = paste0("integrated_snn_res.", i)
      #plot = StackedVlnPlot(temp, gene)
      plot = VlnPlot(data2, features=gene, pt.size=0, stack=T)+ NoLegend()
      p[[i]] = plot
    }
  }
  # 判断细胞类型
  if(TRUE){
    # 结合clustree确定resolution值
    resolution = "0.3"
    Idents(data2) = paste0("integrated_snn_res.", resolution)
    # 结合clustree判断cluster的可靠细胞类型, 及可疑细胞类型
    dfTemp = markersT[["df"]]
    # 结合标记基因在cluster中的平均表达量判断cluster的可疑细胞类型
    temp = list()
    progressLen = length(levels(Idents(data2)))
    for(i in levels(Idents(data2))){
      cat(i, '/', progressLen, '\n')
      temp[[i]] = getCellType(data2, cluster=as.numeric(i), markerGeneList=markerGene[["T"]], debug=TRUE)
    }
    temp = do.call(data.frame, temp)
    temp = t(temp)
    rownames(temp) = gsub("X", "", rownames(temp))
    # 结合VlnPlot验证确定的细胞类型
    p[[resolution]] %>% print()
  }
  # 细胞类型注释
  if(TRUE){
    new.cluster.ids <- c("Th1", "CD4", "NKT", "Th22", "Th17",
                         "Th1", "CTL", "Th2", "Treg", "CD8",
                         "CD8", "CD8")
    names(new.cluster.ids) <- levels(data2)
    data2 <- RenameIdents(data2, new.cluster.ids)
  }
  # 查看细胞类型注释的质量
  if(TRUE){
    # UMAP图 在指定的resolution下
    print(DimPlot(data2, reduction="umap", group.by=paste0("integrated_snn_res.", resolution), pt.size=0.5, label=TRUE))
    # UMAP图 在分类的细胞类型下
    print(DimPlot(data2, reduction="umap", pt.size=0.5, label=TRUE))
    # 热图 在cluster的top10下
    temp = subset(data2, downsample=100)
    top10 = dplyr::top_n(dplyr::group_by(markersDifRes[[resolution]], cluster), n=10, wt=avg_log2FC)
    top10 = as.data.frame(top10)
    top10 = top10[, c("cluster", "gene", "avg_log2FC")]
    print(DoHeatmap(temp, features=top10$gene)+ NoLegend())
    # 热图 在cellType的标记基因下
    gene = unlist(markerGene[["T"]]) %>% unique()
    gene = gene[gene %in% rownames(data2)]
    print(DoHeatmap(temp, features=gene)+ NoLegend())
    # 热图 在cluster的top10 + cellType的标记基因下
    levels(top10$cluster) = new.cluster.ids
    for(cellType in names(markerGene[["T"]])){
      cellType2 = gsub("_.*", "", cellType)
      if(cellType2 %in% levels(top10$cluster)){
        gene = markerGene[["T"]][[cellType]]
        gene = gene[gene %in% rownames(data2)]
        gene = gene[!gene %in% top10$gene]
        for(i in gene){
          top10 = rbind(top10, c(cellType2, i, 0))
        }
      }
    }
    top10 = top10[order(match(top10$cluster, levels(Idents(data2))), -as.numeric(top10$avg_log2FC)),]
    print(DoHeatmap(temp, features=top10$gene)+ NoLegend())
  }
  # 将T细胞亚型注释添加到scRNAMerged中
  if(TRUE){
    debug = scRNAMerged@active.ident %>% as.data.frame()
    colnames(debug) = c("cellType")
    debug["cellType"] = as.character(debug[["cellType"]])
    debug["barcode"] = rownames(debug)
    debug2 = data2@active.ident %>% as.data.frame()
    colnames(debug2) = c("cellType")
    debug2["barcode"] = rownames(debug2)
    debug2[["cellType"]] = as.character(debug2[["cellType"]]) %>% paste0(., "_T")
    debug3 = match(debug2[["barcode"]], debug[["barcode"]])
    debug[debug3, "cellType"] = debug2[["cellType"]]
    scRNAMerged@meta.data$cellType = debug$cellType
    Idents(scRNAMerged) = "cellType"
    DimPlot(scRNAMerged, reduction="umap", pt.size=0.5)
  }
}

saveRDS(scRNAMerged, "0002_scRNAMergedAfterCellType.rds")



# 最终的比例图
if(TRUE){
  # 统计在不同disease中不同类型细胞的占比
  if(TRUE){
    df <- data.frame(barcode=names(scRNAMerged@active.ident),
                     cellType=scRNAMerged@active.ident,
                     disease=scRNAMerged@meta.data$disease)
    df = table(df$cellType, df$disease) %>% addmargins() %>% as.data.frame.array()
    for(c in colnames(df)){
      df[c] = df[c] / df["Sum", c]
    }
    df = df[-which(rownames(df)=="Sum"), -which(colnames(df)=="Sum")]
    df = t(df) %>% as.data.frame()
    #df[["disease"]] = rownames(df)
  }
  df = df*10000
  df = log(df, base=2)
}


if(TRUE){
  df = as.matrix(df)
  fig <- plot_ly(
    x = rownames(df),
    y = colnames(df),
    z = df,
    type = "heatmap"
  )
  fig
}
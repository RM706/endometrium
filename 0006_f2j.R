libPaths(c("/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/",
             "/home/weiyihu/App/conda/envs/RStudio/lib/R/library"
             #"/home/weiyihu/package/share/4.0"
))
setwd("/data1/weiyihu/endometrium/data")
library("plotly", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("future", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("future.apply", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("tidyverse")
.libPaths(c("/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/"))
library("Seurat", lib.loc="/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/")
.libPaths(c("/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/",
            "/home/weiyihu/App/conda/envs/RStudio/lib/R/library"))
library("qpdf", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("esc", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("meta", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("jsonlite", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")


# 前置参数
if(TRUE){
  MAXCORES = 4
  inputFileName = "0004_scRNAMerged.rds"
  inputTFileName = "/data1/shumeng/5.NK_cell/6.TLOs/2.project_V2/2.Tcell/rPCA/Lympht.rename.last.rds"
  inputBFileName = "/data1/shumeng/5.NK_cell/6.TLOs/2.project_V2/3.Bcell/rPCA/Bcells.renames.rds"
  outputFileName = "0004_plot.pdf"
  
}

# 初始化
if(TRUE){
  source("/data1/weiyihu/endometrium/data/share/paraments.R")
  source("/data1/weiyihu/endometrium/data/share/function.R")
  plan(multisession, workers=MAXCORES)
  options(future.globals.maxSize=21*1024*1024*1024,  # 10GB
          scipen=10)  
}

# 加载数据
if(TRUE){
  scRNAMerged = readRDS(inputFileName)  # 未标注细胞类型
  metaInfo = readRDS("0004_metaData.rds")
}

# 绘制F2j基因集表达
if(TRUE){
  # 汇总数据
  dataHeatmap = list()
  # 计算基因表达
  if(TRUE){
    # 从geneSet中获取gene
    geneUsed = c(geneSet$stromal$Angiogenesis,
                 geneSet$stromal$AdhesionMolecule,
                 geneSet$stromal$LymphocyteRecruitment)
    geneUsed = unique(geneUsed)
    geneUsed = geneUsed[geneUsed %in% rownames(scRNAMerged)]
    # 获取Stromal细胞, 先按group进行分组, 再按sample分组, 去除Ctrl, RPL, TE, 添加Ctrl(Endometrium)
    barcodeStromalList = list()
    for(group in unique(scRNAMerged@meta.data[["group"]])){
      df = scRNAMerged@meta.data[(scRNAMerged@meta.data[["group"]]==group) &
                                   (scRNAMerged@meta.data[["cellType"]] %in% cellTypeClassified[["Stromal_cells"]]),
                                 "sample",
                                 drop=FALSE
      ]
      sampleGroup = unique(df[["sample"]])
      temp = list()
      for(sample in sampleGroup){
        temp[[sample]] = rownames(df[df[["sample"]]==sample, , drop=FALSE])
      }
      barcodeStromalList[[group]] = temp
    }
    barcodeStromalList = barcodeStromalList[-which(names(barcodeStromalList) %in% c("Ctrl", "TE", "RPL"))]
    df = scRNAMerged@meta.data[(scRNAMerged@meta.data[["groupSub"]]=="Ctrl(Endometrium)") &
                                 (scRNAMerged@meta.data[["cellType"]] %in% cellTypeClassified[["Stromal_cells"]]),
                               "sample",
                               drop=FALSE]
    temp = list()
    for(sample in unique(df[["sample"]])){
      temp[[sample]] = rownames(df[df[["sample"]]==sample, , drop=FALSE])
    }
    barcodeStromalList[["Ctrl(Endometrium)"]] = temp
    # 计算在每个group的每个样本中, geneUsed的表达
    barcodeAll = lapply(barcodeStromalList, function(x){unlist(x)}) %>% unlist()
    dfData = scRNAMerged@assays$RNA@counts[, barcodeAll]
    numGroup = 0
    numGroupLen = length(barcodeStromalList)
    geneExpression = lapply(barcodeStromalList, function(group){
      numGroup <<- numGroup +1
      numSample = 0
      numSampleLen = length(group)
      lapply(group, function(barcode){
        numSample <<- numSample +1
        cat(numGroup, '/', numGroupLen, '\t', numSample, '/', numSampleLen, '\t', '\r')
        df = dfData[, barcode]
        return(rowSums(df[geneUsed, ]) / sum(df) *10^6)
      })
    })
  }
  # 整理格式
  if(TRUE){
    # 将每个sample的基因表达情况转化为矩阵
    dfList = lapply(geneExpression, function(group){
      df = data.frame(row.names=names(group[[1]]))
      for(sample in names(group)){
        df[[sample]] = group[[sample]]
      }
      return(df)
    })
    
    df = data.frame(row.names=rownames(dfList[[1]]))
    for(dfGroup in dfList){
      df = cbind(df, dfGroup)
    }
    geneDel = apply(df, 1, function(c){sum(ifelse(c==0,1,0)) == length(c)})
    geneDel = names(geneDel)[geneDel]
    df = df[!rownames(df) %in% geneDel, ]
    dataHeatmap$raw = df  # 保存数据
  }
  # 统计数据
  if(TRUE){
    # 对gene在每个group的所有sample中取均值
    dfDebug = data.frame(row.names=rownames(df))
    for(group in names(sampleInGroup[[groupBy]])){
      sample = sampleInGroup[[groupBy]][[group]]
      sample = sample[sample %in% colnames(df)]
      dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
    }
    dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
    dataHeatmap$rawMean = dfDebug  # 保存数据
    
    ctrl = "Ctrl(Endometrium)"
    cutoffValue = NULL
    dfDebug[["sum"]] = apply(dfDebug, 1, function(x){sum(x)})
    dfDebug = apply(dfDebug, 2, function(x){x/dfDebug$sum})
    dfDebug = dfDebug[, -which(colnames(dfDebug)=="sum")]
    dfDebug = as.data.frame(dfDebug)
    dfDebug[, "orderCtrl"] = 1
    col = colnames(dfDebug)[-which(colnames(dfDebug) %in% c(ctrl, "orderCtrl"))]
    ctrl = apply(dfDebug[, ctrl, drop=FALSE], 1, function(x){mean(x)})
    for(i in col){
      score = ifelse(ctrl <= dfDebug[i], 0, 1)
      dfDebug[, "orderCtrl"] = dfDebug[, "orderCtrl"] + score
    }
    
    cutoffValue = ifelse(is.null(cutoffValue), floor(dim(dfDebug)[2]/2), cutoffValue)
    dfDebug = dfDebug[dfDebug[["orderCtrl"]]<=cutoffValue, ]
    
    df = df[rownames(df)[rownames(df) %in% rownames(dfDebug)], ]
    ref <- unlist(geneSetDebug)
    val <- rownames(dfDebug)
    val_sorted <- val[order(match(val, ref))]
    df = df[val_sorted, ]
    dataHeatmap$filter = df
    
    dfDebug = data.frame(row.names=rownames(df))
    for(group in names(sampleInGroup[[groupBy]])){
      sample = sampleInGroup[[groupBy]][[group]]
      sample = sample[sample %in% colnames(df)]
      dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
    }
    dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
    dataHeatmap$filterMean = dfDebug  # 保存数据
  }
  
  saveRDS(dataHeatmap, "0006_dataHeatmap.rds")
}






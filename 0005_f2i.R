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

# 绘制F2i基因集表达
if(TRUE){
  # 常规
  if(FALSE){
    dataHeatmap = list()
    groupBy = "group"
    # 处理免疫细胞
    if(TRUE){
      # 指定基因集
      geneUsed = c(unlist(geneSet$immune))
      #geneUsed = c(geneSet$LymphocyteRecruitment,geneSet$Chemokines,geneSet$Other)
      names(geneUsed) = NULL
      geneUsed = unique(geneUsed)
      geneUsed = geneUsed[geneUsed %in% rownames(scRNAMerged)]
      # 提取免疫细胞的barcode并按group分组
      barcodeImmuneList = list()
      for(group in unique(scRNAMerged@meta.data[[groupBy]])){
        barcodeSampleList = list()
        barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[[groupBy]]==group) & (scRNAMerged@meta.data$cellType %in% cellClustersClassified$immune_cells)] 
        barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
        for(sample in unique(barcodeSample[["sample"]])){
          barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
        }
        barcodeImmuneList[[group]] = barcodeSampleList
      }
      barcodeImmuneList = barcodeImmuneList[-which(names(barcodeImmuneList) %in% c("Ctrl", "TE", "RPL"))]
      barcodeSampleList = list()
      barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[["groupSub"]]=="Ctrl(Endometrium)") & (scRNAMerged@meta.data$cellType %in% cellClustersClassified$immune_cells)] 
      barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
      for(sample in unique(barcodeSample[["sample"]])){
        barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
      }
      barcodeImmuneList[["Ctrl(Endometrium)"]] = barcodeSampleList
      
      # 计算在所有免疫细胞中指定基因集的相对表达值
      barcodeAll = unlist(lapply(barcodeImmuneList, function(group){unlist(group)}))
      dfData = scRNAMerged@assays$RNA@counts[, barcodeAll]
      numGroup = 0
      numGroupLen = length(barcodeImmuneList)
      expressionList_1 = lapply(barcodeImmuneList, function(group){
        numGroup <<- numGroup +1
        numSample = 0
        numSampleLen = length(group)
        
        sampleExpression = lapply(group, function(barcode){
          numSample <<- numSample +1
          cat(numGroup, '/', numGroupLen, '\t', numSample, '/', numSampleLen, '\t', '\r')
          df = dfData[, barcode]
          expressionTemp = apply(df[geneUsed, ], 1, function(x){sum(x)})
          expressionTemp = (expressionTemp / sum(df)) *10^6
        })
        
        return(sampleExpression)
      })
    }
    # 处理GC细胞
    if(TRUE){
      # 提取geneSet中proliferation的gene
      #geneUsed = c(geneSet$Proliferation)
      geneUsed = c(unlist(geneSet$GC))
      names(geneUsed) = NULL
      geneUsed = unique(geneUsed)
      geneUsed = geneUsed[geneUsed %in% rownames(scRNAMerged)]
      # 提取GC细胞的barcode并按group分组
      barcodeGCList = list()
      for(group in unique(scRNAMerged@meta.data[[groupBy]])){
        barcodeSampleList = list()
        barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[[groupBy]]==group) & (scRNAMerged@meta.data$cellTypeSub %in% c("DZ B cells", "LZ B cells"))] 
        barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
        for(sample in unique(barcodeSample[["sample"]])){
          barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
        }
        barcodeGCList[[group]] = barcodeSampleList
      }
      barcodeGCList = barcodeGCList[-which(names(barcodeGCList) %in% c("Ctrl", "TE", "RPL"))]
      
      barcodeSampleList = list()
      barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[["groupSub"]]=="Ctrl(Endometrium)") & (scRNAMerged@meta.data$cellTypeSub %in% c("DZ B cells", "LZ B cells"))] 
      barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
      for(sample in unique(barcodeSample[["sample"]])){
        barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
      }
      barcodeGCList[["Ctrl(Endometrium)"]] = barcodeSampleList
      
      # 计算在所有GC细胞中指定基因集的相对表达值
      barcodeAll = unlist(lapply(barcodeGCList, function(group){unlist(group)}))
      dfData = scRNAMerged@assays$RNA@counts[, barcodeAll]
      numGroup = 0
      numGroupLen = length(barcodeGCList)
      expressionList_2 = lapply(barcodeGCList, function(group){
        numGroup <<- numGroup +1
        numSample = 0
        numSampleLen = length(group)
        
        sampleExpression = lapply(group, function(barcode){
          numSample <<- numSample +1
          cat(numGroup, '/', numGroupLen, '\t', numSample, '/', numSampleLen, '\t', '\r')
          df = dfData[, barcode, drop=FALSE]
          expressionTemp = apply(df[geneUsed, , drop=FALSE], 1, function(x){sum(x)})
          expressionTemp = (expressionTemp / sum(df)) *10^6
        })
        
        return(sampleExpression)
      })
    }
    # 合并免疫细胞与GC细胞的基因表达
    if(TRUE){
      dfList1 = lapply(expressionList_1, function(group){
        df = data.frame(row.names=names(group[[1]]))
        for(sample in names(group)){
          df[[sample]] = group[[sample]]
        }
        return(df)
      })
      dfList2 = lapply(expressionList_2, function(group){
        if(length(group)==0){return(data.frame())}
        df = data.frame(row.names=names(group[[1]]))
        for(sample in names(group)){
          df[[sample]] = group[[sample]]
        }
        return(df)
      })
      #dataHeatmapSample[[groupBy]][["immune"]] = dfList1
      #dataHeatmapSample[[groupBy]][["GC"]] = dfList2
      
      # 合并dfList1/dfList2
      dfList1 = lapply(dfList1, function(df){
        df = as.data.frame(t(df))
        df[, "sample"] = rownames(df)
        return(df)
      })
      dfList2 = lapply(dfList2, function(df){
        df = as.data.frame(t(df))
        df[, "sample"] = rownames(df)
        return(df)
      })
      dfList = list()
      for(group in names(dfList1)){
        df1 = dfList1[[group]]
        df2 = dfList2[[group]]
        df = dplyr::full_join(df1, df2, by="sample")
        df[is.na(df)] = 0
        rownames(df) = df[["sample"]]
        df = df[, -which(colnames(df)=="sample")]
        df = t(df)
        dfList[[group]] = df
      }
      geneUsed = lapply(dfList, function(df){rownames(df)}) %>% unlist() %>% unique()
      dfList = lapply(dfList, function(df){
        gene = setdiff(geneUsed, rownames(df))
        df = as.data.frame(t(df))
        df[gene] = 0
        df = as.data.frame(t(df))
      })
    }
    #dataHeatmapSample[[groupBy]][["raw"]] = dfList 
    # 整理格式
    if(TRUE){
      df = data.frame(row.names=rownames(dfList[[1]]))
      for(dfGroup in dfList){
        df = cbind(df, dfGroup)
      }
      
      geneDel = apply(df, 1, function(c){sum(ifelse(c==0,1,0)) == length(c)})  # 去掉在所有group中表达为0的gene
      geneDel = names(geneDel)[geneDel]  # 去掉在所有group中表达为0的gene
      df = df[!rownames(df) %in% geneDel, ]
      dataHeatmap$raw = df  # 保存数据
    }
    # 根据gene在Ctrl中的表达对gene进行筛选
    if(TRUE){
      #ctrl = ifelse(groupBy=="group", list(c("Ctrl")), list(c("Ctrl(Endometrium)", "Ctrl(Decidual)")))[[1]]
      ctrl = "Ctrl(Endometrium)"
      cutoffValue = NULL
      
      # 对gene在每个group的所有sample中取均值
      dfDebug = data.frame(row.names=rownames(df))
      for(group in names(sampleInGroup[["group2"]])){
        sample = sampleInGroup[["group2"]][[group]]
        sample = sample[sample %in% colnames(df)]
        dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
      }
      dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
      dataHeatmap$rawMean = dfDebug  # 保存数据
      
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
      ref <- unlist(c(geneSet$immune, geneSet$GC))
      val <- rownames(dfDebug)
      val_sorted <- val[order(match(val, ref))]
      df = df[val_sorted, ]
      dataHeatmap$filter = df  # 保存数据
      
      dfDebug = data.frame(row.names=rownames(df))
      for(group in names(sampleInGroup[["group2"]])){
        sample = sampleInGroup[["group2"]][[group]]
        sample = sample[sample %in% colnames(df)]
        dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
      }
      dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
      dataHeatmap$filterMean = dfDebug  # 保存数据
    }
    
    saveRDS(dataHeatmap, "0005_dataHeatmap.rds")
  }
  # counts +1
  if(TRUE){
    dataHeatmap = list()
    groupBy = "group"
    # 处理免疫细胞
    if(TRUE){
      # 指定基因集
      geneUsed = c(unlist(geneSet$immune))
      #geneUsed = c(geneSet$LymphocyteRecruitment,geneSet$Chemokines,geneSet$Other)
      names(geneUsed) = NULL
      geneUsed = unique(geneUsed)
      geneUsed = geneUsed[geneUsed %in% rownames(scRNAMerged)]
      # 提取免疫细胞的barcode并按group分组
      barcodeImmuneList = list()
      for(group in unique(scRNAMerged@meta.data[[groupBy]])){
        barcodeSampleList = list()
        barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[[groupBy]]==group) & (scRNAMerged@meta.data$cellType %in% cellClustersClassified$immune_cells)] 
        barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
        for(sample in unique(barcodeSample[["sample"]])){
          barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
        }
        barcodeImmuneList[[group]] = barcodeSampleList
      }
      barcodeImmuneList = barcodeImmuneList[-which(names(barcodeImmuneList) %in% c("Ctrl", "TE", "RPL"))]
      barcodeSampleList = list()
      barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[["groupSub"]]=="Ctrl(Endometrium)") & (scRNAMerged@meta.data$cellType %in% cellClustersClassified$immune_cells)] 
      barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
      for(sample in unique(barcodeSample[["sample"]])){
        barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
      }
      barcodeImmuneList[["Ctrl(Endometrium)"]] = barcodeSampleList
      
      # 计算在所有免疫细胞中指定基因集的相对表达值
      barcodeAll = unlist(lapply(barcodeImmuneList, function(group){unlist(group)}))
      dfData = scRNAMerged@assays$RNA@counts[, barcodeAll]
      numGroup = 0
      numGroupLen = length(barcodeImmuneList)
      expressionList_1 = lapply(barcodeImmuneList, function(group){
        numGroup <<- numGroup +1
        numSample = 0
        numSampleLen = length(group)
        
        sampleExpression = lapply(group, function(barcode){
          numSample <<- numSample +1
          cat(numGroup, '/', numGroupLen, '\t', numSample, '/', numSampleLen, '\t', '\r')
          df = dfData[, barcode]
          expressionTemp = apply(df[geneUsed, , drop=FALSE], 1, function(x){sum(x)})
          expressionSample = rowSums(df)
          expressionTemp = (expressionTemp / sum(expressionSample)) *10^6
          expressionTemp = expressionTemp +1
          return(expressionTemp)
        })
        
        return(sampleExpression)
      })
    }
    # 处理GC细胞
    if(TRUE){
      # 提取geneSet中proliferation的gene
      #geneUsed = c(geneSet$Proliferation)
      geneUsed = c(unlist(geneSet$GC))
      names(geneUsed) = NULL
      geneUsed = unique(geneUsed)
      geneUsed = geneUsed[geneUsed %in% rownames(scRNAMerged)]
      # 提取GC细胞的barcode并按group分组
      barcodeGCList = list()
      for(group in unique(scRNAMerged@meta.data[[groupBy]])){
        barcodeSampleList = list()
        barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[[groupBy]]==group) & (scRNAMerged@meta.data$cellTypeSub %in% c("DZ B cells", "LZ B cells"))] 
        barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
        for(sample in unique(barcodeSample[["sample"]])){
          barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
        }
        barcodeGCList[[group]] = barcodeSampleList
      }
      barcodeGCList = barcodeGCList[-which(names(barcodeGCList) %in% c("Ctrl", "TE", "RPL"))]
      
      barcodeSampleList = list()
      barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[["groupSub"]]=="Ctrl(Endometrium)") & (scRNAMerged@meta.data$cellTypeSub %in% c("DZ B cells", "LZ B cells"))] 
      barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
      for(sample in unique(barcodeSample[["sample"]])){
        barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
      }
      barcodeGCList[["Ctrl(Endometrium)"]] = barcodeSampleList
      
      # 计算在所有GC细胞中指定基因集的相对表达值
      barcodeAll = unlist(lapply(barcodeGCList, function(group){unlist(group)}))
      dfData = scRNAMerged@assays$RNA@counts[, barcodeAll]
      numGroup = 0
      numGroupLen = length(barcodeGCList)
      expressionList_2 = lapply(barcodeGCList, function(group){
        numGroup <<- numGroup +1
        numSample = 0
        numSampleLen = length(group)
        
        sampleExpression = lapply(group, function(barcode){
          numSample <<- numSample +1
          cat(numGroup, '/', numGroupLen, '\t', numSample, '/', numSampleLen, '\t', '\r')
          df = dfData[, barcode, drop=FALSE]
          
          expressionTemp = apply(df[geneUsed, , drop=FALSE], 1, function(x){sum(x)})
          expressionSample = rowSums(df)
          expressionTemp = (expressionTemp / sum(expressionSample)) *10^6
          expressionTemp = expressionTemp +1
          return(expressionTemp)
        })
        
        return(sampleExpression)
      })
    }
    # 合并免疫细胞与GC细胞的基因表达
    if(TRUE){
      dfList1 = lapply(expressionList_1, function(group){
        df = data.frame(row.names=names(group[[1]]))
        for(sample in names(group)){
          df[[sample]] = group[[sample]]
        }
        return(df)
      })
      dfList2 = lapply(expressionList_2, function(group){
        if(length(group)==0){return(data.frame())}
        df = data.frame(row.names=names(group[[1]]))
        for(sample in names(group)){
          df[[sample]] = group[[sample]]
        }
        return(df)
      })
      
      # 合并dfList1/dfList2
      dfList1 = lapply(dfList1, function(df){
        df = as.data.frame(t(df))
        df[, "sample"] = rownames(df)
        return(df)
      })
      dfList2 = lapply(dfList2, function(df){
        df = as.data.frame(t(df))
        df[, "sample"] = rownames(df)
        return(df)
      })
      dfList = list()
      for(group in names(dfList1)){
        df1 = dfList1[[group]]
        df2 = dfList2[[group]]
        df = dplyr::full_join(df1, df2, by="sample")
        df[is.na(df)] = 1  # 注意这个值, 代表了在相应样本中未检测到相应基因时的值(指定样本不存在GC细胞时, proliferation的gene的表达值)
        rownames(df) = df[["sample"]]
        df = df[, -which(colnames(df)=="sample")]
        df = t(df)
        dfList[[group]] = df
      }
      geneUsed = lapply(dfList, function(df){rownames(df)}) %>% unlist() %>% unique()
      dfList = lapply(dfList, function(df){
        gene = setdiff(geneUsed, rownames(df))
        df = as.data.frame(t(df))
        df[gene] = 1  # 注意这个值, 代表了在相应样本中未检测到相应基因时的值(指定group不存在GC细胞时, proliferation的gene的表达值)
        df = as.data.frame(t(df))
      })
    }
    # 整理格式
    if(TRUE){
      df = data.frame(row.names=rownames(dfList[[1]]))
      for(dfGroup in dfList){
        df = cbind(df, dfGroup)
      }
      
      geneDel = apply(df, 1, function(c){sum(ifelse(c==0,1,0)) == length(c)})  # 去掉在所有group中表达为0的gene
      geneDel = names(geneDel)[geneDel]  # 去掉在所有group中表达为0的gene
      df = df[!rownames(df) %in% geneDel, ]
      dataHeatmap$raw = df  # 保存数据
    }
    # 根据gene在Ctrl中的表达对gene进行筛选
    if(TRUE){
      #ctrl = ifelse(groupBy=="group", list(c("Ctrl")), list(c("Ctrl(Endometrium)", "Ctrl(Decidual)")))[[1]]
      ctrl = "Ctrl(Endometrium)"
      cutoffValue = NULL
      
      # 对gene在每个group的所有sample中取均值
      dfDebug = data.frame(row.names=rownames(df))
      for(group in names(sampleInGroup[["group2"]])){
        sample = sampleInGroup[["group2"]][[group]]
        sample = sample[sample %in% colnames(df)]
        dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
      }
      dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
      dataHeatmap$rawMean = dfDebug  # 保存数据
      
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
      ref <- unlist(c(geneSet$immune, geneSet$GC))
      val <- rownames(dfDebug)
      val_sorted <- val[order(match(val, ref))]
      df = df[val_sorted, ]
      dataHeatmap$filter = df  # 保存数据
      
      dfDebug = data.frame(row.names=rownames(df))
      for(group in names(sampleInGroup[["group2"]])){
        sample = sampleInGroup[["group2"]][[group]]
        sample = sample[sample %in% colnames(df)]
        dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
      }
      dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
      dataHeatmap$filterMean = dfDebug  # 保存数据
    }
    
    saveRDS(dataHeatmap, "0005_dataHeatmap.rds")
  }
  # counts +1 且 只处理免疫细胞
  if(FALSE){
    dataHeatmap = list()
    groupBy = "group"
    # 处理免疫细胞
    if(TRUE){
      # 指定基因集
      geneUsed = c(unlist(geneSet$immune))
      #geneUsed = c(geneSet$LymphocyteRecruitment,geneSet$Chemokines,geneSet$Other)
      names(geneUsed) = NULL
      geneUsed = unique(geneUsed)
      geneUsed = geneUsed[geneUsed %in% rownames(scRNAMerged)]
      # 提取免疫细胞的barcode并按group分组
      barcodeImmuneList = list()
      for(group in unique(scRNAMerged@meta.data[[groupBy]])){
        barcodeSampleList = list()
        barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[[groupBy]]==group) & (scRNAMerged@meta.data$cellType %in% cellClustersClassified$immune_cells)] 
        barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
        for(sample in unique(barcodeSample[["sample"]])){
          barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
        }
        barcodeImmuneList[[group]] = barcodeSampleList
      }
      barcodeImmuneList = barcodeImmuneList[-which(names(barcodeImmuneList) %in% c("Ctrl", "TE", "RPL"))]
      barcodeSampleList = list()
      barcode = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data[["groupSub"]]=="Ctrl(Endometrium)") & (scRNAMerged@meta.data$cellType %in% cellClustersClassified$immune_cells)] 
      barcodeSample = scRNAMerged@meta.data[barcode, "sample", drop=FALSE]
      for(sample in unique(barcodeSample[["sample"]])){
        barcodeSampleList[[sample]] = rownames(barcodeSample)[barcodeSample[["sample"]]==sample]
      }
      barcodeImmuneList[["Ctrl(Endometrium)"]] = barcodeSampleList
      
      # 计算在所有免疫细胞中指定基因集的相对表达值
      barcodeAll = unlist(lapply(barcodeImmuneList, function(group){unlist(group)}))
      dfData = scRNAMerged@assays$RNA@counts[, barcodeAll]
      numGroup = 0
      numGroupLen = length(barcodeImmuneList)
      expressionList_1 = lapply(barcodeImmuneList, function(group){
        numGroup <<- numGroup +1
        numSample = 0
        numSampleLen = length(group)
        
        sampleExpression = lapply(group, function(barcode){
          numSample <<- numSample +1
          cat(numGroup, '/', numGroupLen, '\t', numSample, '/', numSampleLen, '\t', '\r')
          df = dfData[, barcode]
          expressionTemp = apply(df[geneUsed, , drop=FALSE], 1, function(x){sum(x)})
          expressionTemp = expressionTemp +1
          expressionSample = rowSums(df)
          expressionSample = expressionSample +1
          expressionTemp = (expressionTemp / sum(expressionSample)) *10^6
          return(expressionTemp)
        })
        
        return(sampleExpression)
      })
    }
    # 
    if(TRUE){
      dfList1 = lapply(expressionList_1, function(group){
        df = data.frame(row.names=names(group[[1]]))
        for(sample in names(group)){
          df[[sample]] = group[[sample]]
        }
        return(df)
      })
    }
    #dataHeatmapSample[[groupBy]][["raw"]] = dfList 
    # 整理格式
    if(TRUE){
      df = data.frame(row.names=rownames(dfList1[[1]]))
      for(dfGroup in dfList1){
        df = cbind(df, dfGroup)
      }
      
      geneDel = apply(df, 1, function(c){sum(ifelse(c==0,1,0)) == length(c)})  # 去掉在所有group中表达为0的gene
      geneDel = names(geneDel)[geneDel]  # 去掉在所有group中表达为0的gene
      df = df[!rownames(df) %in% geneDel, ]
      dataHeatmap$raw = df  # 保存数据
    }
    # 根据gene在Ctrl中的表达对gene进行筛选
    if(TRUE){
      #ctrl = ifelse(groupBy=="group", list(c("Ctrl")), list(c("Ctrl(Endometrium)", "Ctrl(Decidual)")))[[1]]
      ctrl = "Ctrl(Endometrium)"
      cutoffValue = NULL
      
      # 对gene在每个group的所有sample中取均值
      dfDebug = data.frame(row.names=rownames(df))
      for(group in names(sampleInGroup[["group2"]])){
        sample = sampleInGroup[["group2"]][[group]]
        sample = sample[sample %in% colnames(df)]
        dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
      }
      dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
      dataHeatmap$rawMean = dfDebug  # 保存数据
      
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
      ref <- unlist(c(geneSet$immune, geneSet$GC))
      val <- rownames(dfDebug)
      val_sorted <- val[order(match(val, ref))]
      df = df[val_sorted, ]
      dataHeatmap$filter = df  # 保存数据
      
      dfDebug = data.frame(row.names=rownames(df))
      for(group in names(sampleInGroup[["group2"]])){
        sample = sampleInGroup[["group2"]][[group]]
        sample = sample[sample %in% colnames(df)]
        dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
      }
      dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
      dataHeatmap$filterMean = dfDebug  # 保存数据
    }
    
    saveRDS(dataHeatmap, "0005_dataHeatmap.rds")
  }
}






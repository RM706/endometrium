.libPaths(c("/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/",
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

# 一级函数
readyForForest = function(df, group1, group2){
  # input:
  #   df|\t  data.frame|\t  长数据, 第一列为geneSet, 第二列为组名group, 第三列为基因表达值value
  #   group1|\t  str|\t  第一组的组名
  #   group2|\t  str|\t  第二组的组名
  df = df
  group1 = group1
  group2 = group2
  
  
  result = data.frame(geneSet=NULL, n1=NULL, n2=NULL, t=NULL, es=NULL, se=NULL, p=NULL)
  for (gene in unique(df[["geneSet"]])){
    valueGroup1 = df[df[["geneSet"]]==gene & df[["group"]]==group1, "value"]
    valueGroup2 = df[df[["geneSet"]]==gene & df[["group"]]==group2, "value"]
    n1 = length(valueGroup1)
    n2 = length(valueGroup2)
    
    t = t.test(valueGroup1, valueGroup2) 
    t.stat = unname(t[["statistic"]])
    t.p = t$p.value
    
    es.obj = esc_t(t=t.stat,
                   grp1n=n1,
                   grp2n=n2,
                   es.type="g")
    es = es.obj[["es"]]
    se = es.obj[["se"]]
    tmp = data.frame(geneSet=gene, n1=n1, n2=n2, t=t.stat, es=es, se=se, p=t.p)
    result = rbind(result, tmp, stringsAsFactors=F)
  }
  
  return(result)
}

# 一级函数
readyForMeta = function(seuratObject, geneSet, sampleBy, groupBy, filterMarker=TRUE){
  # input:
  #   seuratObject
  #   geneSet
  seuratObject = seuratObject
  geneSet = geneSet
  sampleBy = sampleBy
  groupBy = groupBy
  filterMarker = filterMarker
  
  if(filterMarker){geneSet = geneSet[geneSet %in% rownames(seuratObject)]}
  data = seuratObject@assays[["rPCA"]]@data[geneSet, ] %>% as.data.frame() %>% t() %>% as.data.frame()
  data[["barcode"]] = rownames(data)
  
  dataMeta = seuratObject@meta.data[, c(sampleBy, groupBy), drop=FALSE]
  colnames(dataMeta) = c("sample", "group")
  dataMeta[["barcode"]] = rownames(dataMeta)
  
  data = left_join(data, dataMeta, by="barcode")
  data = gather(data, key="geneSet", value="value", -c("barcode", "group", "sample"))
  data = split(data, data[["geneSet"]])
  
  return(data)
}

# 一级函数
metaAnalysis = function(data, group1, group2){
  data = data
  group1 = group1
  group2 = group2
  
  statistic = lapply(data, function(dfGene){
    expressionGroup1 = dfGene[dfGene[["group"]]==group1, "value"]
    expressionGroup2 = dfGene[dfGene[["group"]]==group2, "value"]
    len1 = length(expressionGroup1) %>% as.numeric()
    len2 = length(expressionGroup2) %>% as.numeric()
    meanGroup1 = mean(expressionGroup1)
    meanGroup2 = mean(expressionGroup2)
    sdGroup1 = sd(expressionGroup1)
    sdGroup2 = sd(expressionGroup2)
    se <- sqrt((sdGroup1^2) / length(expressionGroup1) + (sdGroup2^2) / length(expressionGroup2))
    log2FC = log(meanGroup1/meanGroup2, base=2)
    # test
    t = t.test(expressionGroup1, expressionGroup2)
    t.stat = t$statistic %>% unname
    t.p = t$p.value
    es.obj = esc_t(t=t.stat, grp1n=len1, grp2n=len2, es.type="d")  # g
    es = es.obj$es
    se = es.obj$se
    return(list(se=se,
                log2FC=log2FC,
                p=t$p.value))
  })
  
  se = c()
  log2FC = c()
  for(i in statistic){
    se = append(se, i[["se"]])
    log2FC = append(log2FC, i[["log2FC"]])
  }
  
  meta_result <- metagen(
    TE = log2FC,  # 效应大小（log2FC）
    seTE = se,  # 标准误差
    studlab = names(statistic),  # 研究标识符
    comb.fixed = FALSE,  # 随机效应模型
    comb.random = TRUE,   # 允许随机效应
    # test
    data=test,
    prediction=F
  )
  
  return(meta_result)
}

# 一级函数
computeLog2FC = function(sampleGroup1, sampleGroup2, pCutoff=0.05){
  # input:
  #   sampleGroup1|\t  list|\t       key为sample名, value为double类型的基因表达值
  #   sampleGroup2|\t  list|\t       key为sample名, value为double类型的基因表达值
  #   pCutoff|\t       double|\t     p值阈值
  sampleGroup1 = sampleGroup1
  sampleGroup2 = sampleGroup2
  pCutoff = pCutoff
  
  sampleGroup1 = Filter(function(x){!is.null(x)}, sampleGroup1)
  sampleGroup2 = Filter(function(x){!is.null(x)}, sampleGroup2)
  geneGroup1 = lapply(sampleGroup1, function(geneExpression){names(geneExpression)}) %>%
    unlist() %>%
    unique()
  geneGroup2 = lapply(sampleGroup2, function(geneExpression){names(geneExpression)}) %>%
    unlist() %>%
    unique()
  geneList = list()
  for(gene in unique(c(geneGroup1, geneGroup2))){
    geneList[[gene]] = gene
  }
  geneList = future_lapply(geneList, function(gene){
    expressionGroup1 = lapply(sampleGroup1, function(geneExpression){
      if(gene %in% names(geneExpression)){return(geneExpression[gene])}else{return(0)}
    }) %>% unlist()
    expressionGroup2 = lapply(sampleGroup2, function(geneExpression){
      if(gene %in% names(geneExpression)){return(geneExpression[gene])}else{return(0)}
    }) %>% unlist()
    wilcoxResult = wilcox.test(expressionGroup1, expressionGroup2)
    meanExpressionGroup1 = mean(expressionGroup1)
    meanExpressionGroup2 = mean(expressionGroup2)
    log2FC = log2(meanExpressionGroup1/meanExpressionGroup2)
    
    result = list(wilcoxResult=wilcoxResult, log2FC=log2FC)
    return(result)
  })
  df = lapply(geneList, function(d){
    result = c(p=d$wilcoxResult$p.value,
               w=as.numeric(d$wilcoxResult$statistic),
               log2FC=d$log2FC)
    return(result)
  }) %>% as.data.frame() %>% t() %>% as.data.frame()
  df = df[df['p']<=pCutoff, ]
  df = df[order(df['p']), ]
  
  return(df)
}

# 二级函数
compareLog2FC = function(seuratObject, sampleList, group1, group2, pCutoff=0.05){
  seuratObject=seuratObject
  sampleList = sampleList
  group1 = group1
  group2 = group2
  pCutoff = pCutoff
  
  sampleGroup1 = seuratObject@meta.data[seuratObject@meta.data[["group"]]==group1, "sample.x"] %>%
    unique() %>%
    sampleList[.]
  sampleGroup2 = seuratObject@meta.data[seuratObject@meta.data[["group"]]==group2, "sample.x"] %>%
    unique() %>%
    sampleList[.]
  df = computeLog2FC(sampleGroup1=sampleGroup1, sampleGroup2=sampleGroup2, pCutoff=pCutoff)
  df = df[, c("p", "log2FC")]
  df[["group1"]] = group1
  df[["group2"]] = group2
  
  return(df)
}

# 初始化
if(TRUE){
  source("/data1/weiyihu/endometrium/data/share/paraments.R")
  source("/data1/weiyihu/endometrium/data/share/function.R")
  plan(multisession, workers=MAXCORES)
  options(future.globals.maxSize=21*1024*1024*1024,  # 10GB
          scipen=10)  
}



# main
if(TRUE){
  # 加载数据
  if(TRUE){
    scRNAMerged = readRDS(inputFileName)  # 未标注细胞类型
  }
  # meta.data中添加cellType信息
  #scRNAMerged@meta.data[["cellType"]] = scRNAMerged@active.ident
  #seuratObject = subset(scRNAMerged, cellType!="Proliferative cells")
  #seuratObject@meta.data$cellType = as.character(seuratObject@meta.data$cellType)

  # 统计在不同类型疾病中的不同类型细胞的比例
  if(TRUE){
    # 准备数据
    if(TRUE){
      df = seuratObject@meta.data
      #df = df[df[["type"]]!="CD45+", ]
      df = table(df[["cellType"]], df[["group"]]) %>%
        as.data.frame() %>%
        spread(., key="Var2", value="Freq")
      rownames(df) = df[["Var1"]]
      df = df[,-1]
      df = df[c("Stromal","Unciliated Epithelial","NK","Endothelial","Macrophages",
                "CD8+ T","CD4+ T","Monocytes","Ciliated Epithelial","B cells","ILC3","DC","Mast"), ]
      
      dfCellType = df %>% mutate_at(.vars=vars(everything()), list(function(x){x/sum(x)}))
      dfCellType = dfCellType %>% select(order(dfCellType["Stromal", ], decreasing=TRUE))
      #dfCellType = dfCellType %>% select("Ctrl", everything())
      dfCellType = dfCellType %>% select("Ctrl", "AS", "EMs","EC", "RPL")
      dfCellType = dfCellType * 100
    }
    # 绘图
    if(TRUE){
      plotCellType <- plot_ly(x=colnames(dfCellType),
                              y=as.numeric(dfCellType[1, ]),
                              type="bar",
                              name=rownames(dfCellType)[1],
                              marker=list(color=cellTypeColors[["base"]][[rownames(dfCellType)[1]]]),
                              width=400, height=500
                              )
      for(i in seq(2, dim(dfCellType)[1])){
        plotCellType = plotCellType %>% add_trace(y=as.numeric(dfCellType[i, ]),
                                                  name=rownames(dfCellType)[i],
                                                  marker=list(color=cellTypeColors[["base"]][[rownames(dfCellType)[i]]]))
      }
      plotCellType <- plotCellType %>% layout(barmode='stack',
                                              xaxis=list(categoryorder="trace",
                                                         tickangle="270",
                                                         showline=TRUE,
                                                         showticklabels=TRUE,
                                                         ticklen=5),
                                              yaxis=list(title=list("text"='Percentage(%)',
                                                                    "font"=list("color"="black",
                                                                                "family"="Arial",
                                                                                "size"=12)),
                                                         showline=TRUE,
                                                         showticklabels=TRUE,
                                                         tickmode="linear",
                                                         dtick=25,
                                                         ticklen=5),
                                              font=list(family="Arial",
                                                        color="black",
                                                        size=12),
                                              paper_bgcolor="white",
                                              plot_bgcolor="white")
    }
  }
  # 统计在不同类型疾病中的细胞类群的比例
  if(TRUE){
    # 准备数据
    if(TRUE){
      dfTemp = t(df) %>%
        as.data.frame()
      dfClusters = data.frame(row.names=rownames(dfTemp))
      for(i in colnames(dfTemp)){
        for(cellClusters in names(cellTypeClassified)){
          if(i %in% cellTypeClassified[[cellClusters]]){
            if(cellClusters %in% colnames(dfClusters)){
              dfClusters[[cellClusters]] = dfClusters[cellClusters] + dfTemp[[i]]
            }else{
              dfClusters[[cellClusters]] = dfTemp[[i]]
            }
          }
        }
      }
      dfClusters = t(dfClusters) %>% as.data.frame()
      dfClusters = dfClusters %>% select(colnames(dfCellType))
      dfClusters = dfClusters %>% mutate_at(.vars=vars(everything()), .funs=list(function(x){x/sum(x)}))
      dfClusters = dfClusters * 100
    }
    # 绘图
    if(TRUE){
      plotClusters = plot_ly(x=colnames(dfClusters),
                             y=as.numeric(dfClusters[1, ]),
                             type="bar",
                             name=rownames(dfClusters)[1],
                             markers=list(),
                             width=400, height=500
      )
      for(i in seq(2, dim(dfClusters)[1])){
        plotClusters = plotClusters %>% add_trace(y=as.numeric(dfClusters[i, ]),
                                                  name=rownames(dfClusters)[i],
                                                  marker=list())
      }
      plotClusters <- plotClusters %>% layout(barmode='stack',
                                              xaxis=list(categoryorder="trace",
                                                         tickangle="270",
                                                         showline=TRUE,
                                                         showticklabels=TRUE,
                                                         ticklen=5),
                                              yaxis=list(title='Percent(%)',
                                                         showline=TRUE,
                                                         showticklabels=TRUE,
                                                         tickmode="linear",
                                                         dtick=25,
                                                         ticklen=5),
                                              font=list(family="Arial",
                                                        color="black",
                                                        size=12),
                                              paper_bgcolor="white",
                                              plot_bgcolor="white")
    }
  }
  # 统计在不同类型疾病中免疫细胞与基质细胞的比例
  if(TRUE){
    # 准备数据
    if(TRUE){
      dfTemp = t(df) %>% as.data.frame()
      df3 = data.frame(row.names=rownames(dfTemp))
      for(i in colnames(dfTemp)){
        for(cellClusters in names(cellClustersClassified)){
          if(i %in% cellClustersClassified[[cellClusters]]){
            if(cellClusters %in% colnames(df3)){
              df3[[cellClusters]] = df3[cellClusters] + dfTemp[[i]]
            }else{
              df3[[cellClusters]] = dfTemp[[i]]
            }
          }
        }
      }
      df3 = t(df3) %>% as.data.frame()
      df3 = df3 %>% select(colnames(dfCellType))
      df3 = df3 %>% mutate_at(.vars=vars(everything()), .funs=list(function(x){x/sum(x)}))
      df3 = df3 * 100
    }
    # 绘图
    if(TRUE){
      plotLine = plot_ly(x=colnames(df3),
                         y=as.numeric(df3[1, ]),
                         name=gsub("_", " ", rownames(df3)[1]),
                         type="scatter",
                         mode="lines+markers",
                         width=600, height=400)
      plotLine = plotLine %>% add_trace(y=as.numeric(df3[2, ]),
                                        yaxis="y2",
                                        name=gsub("_", " ", rownames(df3)[2]),
                                        mode="lines+markers")
      plotLine <- plotLine %>% layout(margin=list('l'=50, 'r'=70, 't'=60, 'b'=50, 'pad'=10),
                                      xaxis=list(categoryorder="trace",
                                                 tickangle="270",
                                                 showline=TRUE,
                                                 showticklabels=TRUE,
                                                 ticklen=5,
                                                 showgrid=FALSE),
                                      yaxis=list(title="Percentage of<br>immune cells (%)",
                                                 showline=TRUE,
                                                 showticklabels=TRUE,
                                                 range=c(floor(min(df3[1, ])), ceiling(max(df3[1, ]))),
                                                 tickmode="linear",
                                                 dtick=10,
                                                 ticklen=5,
                                                 showgrid=TRUE),
                                      yaxis2=list(overlaying="y",
                                                  side="right",
                                                  title="Percentage of<br>stromal cells (%)",
                                                  range=c(floor(min(df3[2, ])), ceiling(max(df3[2, ]))),
                                                  showline=TRUE,
                                                  showticklabels=TRUE,
                                                  tickmode="linear",
                                                  dtick=5,
                                                  ticklen=5,
                                                  showgrid=TRUE
                                      ),
                                      legend=list(x=0.3, y=1.2,
                                                  orientation='h'),
                                      font=list(family="Arial",
                                                color="black",
                                                size=12),
                                      paper_bgcolor="white",
                                      plot_bgcolor="white"
      )
    }
  }
  # 绘制不同类型疾病中的细胞类型UMAP图
  if(TRUE){
    group_UMAP <- factor(seuratObject@meta.data$group, c("Ctrl", "TE", "AS", "EMs", "EC", "RPL"))
    seuratObject@meta.data$group_UMAP = group_UMAP
    pdf("temp4.pdf", width=15, height=10)
    DimPlot(seuratObject, reduction="umap", split.by="group_UMAP", group.by="cellType", ncol=3,
            cols=getCellTypeColor(cellType=seuratObject@meta.data[["cellType"]], colorList=cellTypeColors[["base"]]))
    dev.off()
  }
  # 统计在不同类型疾病中Tfh, Memory B, GC B, Plasma B细胞相对于健康组的细胞占比的log2FC
  if(TRUE){
    # 准备数据
    if(TRUE){
      df = seuratObject@meta.data
      df = table(df$cellTypeSub, df$group)
      df = as.data.frame(df)
      df = spread(df, key="Var2", value="Freq")
      rownames(df) = df[, 1]
      df = df[, -1]
      df["sum", ] = apply(df, 2, function(c){sum(c)})
      #temp = apply(df, 1, function(r){(r / (df["sum", ]-r))})
      temp = apply(df, 1, function(r){(r / (df["sum", ]))})
      df = df[0, ]
      for(i in names(temp)){
        tempDf = temp[[i]]
        rownames(tempDf) = i
        df = rbind(df, tempDf)
      }
      #for(i in colnames(df)){
      #  if(i=="Ctrl"){next}
      #  df[i] = log((df[i] / df["Ctrl"]), base=2)
      #}
      #df = df[-which(rownames(df)=="sum"), -which(colnames(df)=="Ctrl")]
      df = df[-which(rownames(df)=="sum"), ]
      
      # 筛选细胞类型
      df = df[c("Tfh", "Memory B cells", "DZ B cells", "LZ B cells", "Plasma cells"), ]
      rownames(df) = c("Tfh", "Memory B", "DZ B", "LZ B", "Plasma B")
      
      #df = df * 100
      #df = as.matrix(df)
      #df[is.infinite(df)] = min(df[!is.infinite(df)])  # 改变inf类型的值为df中的最小值
      #df = as.data.frame(df)
    }
    # 绘图
    if(TRUE){
      colorBarRange = abs(min(df)) + abs(max(df))
      colorBar0Value = abs(min(df)) / colorBarRange
      plotCellProportion <- plot_ly(x=colnames(df),
                                    y=rownames(df),
                                    z=as.matrix(df),
                                    type="heatmap",
                                    colorscale=data.frame(c(0, colorBar0Value, 1),
                                                          c("darkblue", "white", "darkred")),
                                    width=500, height=400)
      plotCellProportion = plotCellProportion# %>%
        #add_annotations(text=" ",
        #                x=-1, y=0.5,
        #                xref="x", yref="paper",
        #                textangle=270,
        #                font=list(family="Arial", size=12, color="black"),
        #                showarrow=FALSE)
      plotCellProportion = plotCellProportion %>%
        layout(margin=list('l'=20, 'r'=20, 't'=20, 'b'=20),
               font=list(family="Arial",
                         color="black",
                         size=12),
               xaxis=list(side="top",
                          tickangle=270,
                          ticklen=0,
                          showgrid=FALSE,
                          zeroline=0),
               yaxis=list(side="left",
                          title=" ",
                          ticklen=0,
                          showgrid=FALSE,
                          zeroline=0),
               legend=list(x=0, y=0)
        )
    }
  }
  # 保存图片
  if(TRUE){
    reticulate::py_run_string("import sys")
    save_image(plotCellType, file="temp1.pdf")
    save_image(plotClusters, file="temp2.pdf")
    save_image(plotLine, file="temp3.pdf")
    save_image(plotCellProportion, file="temp5.pdf")
    
    
    qpdf::pdf_combine(c("temp1.pdf", "temp2.pdf", "temp3.pdf", "temp4.pdf",
                        "temp5.pdf"),
                      output=outputFileName)
    file.remove(c("temp1.pdf", "temp2.pdf", "temp3.pdf", "temp4.pdf",
                  "temp5.pdf"))
  }
  
  
  # 已废弃, 新code存在于0005_f2i.R中
  # 准备数据, 绘制在不同group中的基因表达热图
  if(FALSE){
    # 常规
    if(FALSE){
      #dataHeatmapSample = list()
      groupBy = "groupSub"
      #dataHeatmapSample[[groupBy]] = list()
      # 处理免疫细胞
      if(TRUE){
        # 提取geneSet中非proliferation的gene
        geneUsed = c("CCR7", "SELL",
                     "ACKR3","CCL11","CCL16","CCL27","CCR10","CCR4","CCR5","CCR6","CD6",
                     "CXCL2","CXCL6","CXCL8","FBLN7","IL1R1","IL22","MADCAM1","MSMP","PF4V1",
                     "SDC1","STAT1","TFF2","TRAF6","CXCL13", "CCL19"
                     )
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
        geneUsed = c("IGF1R", "GNB4", "GNB1", "MAP2K1")
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
        dfList = dfList[-which(names(dfList) %in% c("RPL", "TE"))]
        df = data.frame(row.names=rownames(dfList[[1]]))
        for(dfGroup in dfList){
          df = cbind(df, dfGroup)
        }
        
        geneDel = apply(df, 1, function(c){sum(ifelse(c==0,1,0)) == length(c)})  # 去掉在所有group中表达为0的gene
        geneDel = names(geneDel)[geneDel]  # 去掉在所有group中表达为0的gene
        df = df[!rownames(df) %in% geneDel, ]
      }
      # 根据gene在Ctrl中的表达对gene进行筛选
      if(TRUE){
        #ctrl = ifelse(groupBy=="group", list(c("Ctrl")), list(c("Ctrl(Endometrium)", "Ctrl(Decidual)")))[[1]]
        ctrl = ifelse(groupBy=="group", list(c("Ctrl")), list(c("Ctrl(Endometrium)")))[[1]]
        cutoffValue = NULL
        
        # 对gene在每个group的所有sample中取均值
        dfDebug = data.frame(row.names=rownames(df))
        for(group in names(sampleInGroup[[groupBy]])){
          sample = sampleInGroup[[groupBy]][[group]]
          sample = sample[sample %in% colnames(df)]
          dfDebug[group] = apply(df[, sample], 1, function(x){mean(x)})
        }
        dfDebug = dfDebug[, colSums(is.na(dfDebug)) < nrow(dfDebug)]  # 删除全部为NA的列
        
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
      }
    }
    # 读取已保存数据
    if(FALSE){
      dataHeatmapSample = readRDS("0004_dataHeatmapSample.rds")
    }
}

#########################
#########################
#########################
#########################

# debug
if(TRUE){
  temp = subset(scRNAMerged, downsample=100)
  temp = scRNAMerged[,sample(dim(scRNAMerged)[2], 10000, replace=FALSE)]
  
  geneSet = c("CCR7", "CXCR5", "SELL", "CCR5", "Itgb2",
              "CD6", "CXCR13", "VCAN", "CXCR5", "CCL19",
              "SDC2", "ITGA6", "CCR6", "SDC3", "CXCR3",
              "CXCR6", "CXCR4", "MAP2K1", "GNB1", "EGFR",
              "GNG12", "EPHA2", "GNB4", "RAPGEF5", "TIAM1",
              "IGF1R")
  sampleBy = "sample.x"
  groupBy = "healthy"
  seuratObject = scRNAMerged
  filterMarker = TRUE
  group1 = "RPL"
  group2 = "Ctrl"
  
  data = readyForMeta(seuratObject=scRNAMerged,
                      geneSet=geneSet,
                      sampleBy=sampleBy,
                      groupBy=groupBy)
  meta_result = metaAnalysis(data=data,
                             group1=group1,
                             group2=group2)
  summary(meta_result)
  forest(meta_result)
  
  
  
}

if(TRUE){
  total = sample(seq(1,1000000), 10000, replace=TRUE)
  totalMean = mean(total)
  totalSd = sd(total)
  
  d = 0
  for(i in seq(1,10000,1)){
    sample1 = sample(total, 100, replace=FALSE)
    sampleMean = mean(sample1)
    sampleSd = sd(sample1)
    sampleSe = sampleSd/sqrt(100)
    if((sampleMean-1.96*sampleSe)<=totalMean & totalMean<=(sampleMean+1.96*sampleSe)){
      d = d+1
    }
  }
  
}

if(TRUE){
  barcodeTarget = rownames(scRNAMerged@meta.data)[scRNAMerged@meta.data[["cellType"]] %in% cellClustersClassified[["immune_cells"]]]
  temp = scRNAMerged[, barcodeTarget]
  seuratObject = temp
  group = unique(seuratObject@meta.data[["group"]])
  group = group[-which(group=="Ctrl")]
  markers = list()
  timeStart = as.numeric(Sys.time())
  for(i in seq(1, length(group))){
    if(i!=1){
      timeCurrent = as.numeric(Sys.time())
      timeUsed = timeCurrent - timeStart
      timeMean = timeUsed / (i-1)
      timeRemained = timeMean * (length(group)-i+1)
      cat(paste0(i, '/', length(group), " -- ", floor(timeRemained/60)," min ", round(timeRemained%%60), " s"), '\n')
    }
    groupTemp = group[i]
    markers[[groupTemp]] = FindMarkers(seuratObject,
                                       group.by="group",
                                       ident.1=groupTemp,
                                       ident.2="Ctrl",
                                       min.pct=0,
                                       logfc.threshold=0)
  }
  saveRDS(markers, "0004_markersImmune.rds")
  for(df in markers){
    print(df[geneSet, "avg_log2FC"])
  }
  
  
}

if(TRUE){
  seuratObject = scRNAMerged
  geneSet = list(lymphocyte_migration=c('ABL1', 'ABL2', 'ADAM10', 'ADAM17', 'ADAM8', 'ADTRP', 'AIF1', 'AIRE', 'AKT1', 'APOD', 'APP', 'ARTN', 'ASCL2', 'CCL1', 'CCL11', 'CCL13', 'CCL14', 'CCL15', 'CCL16', 'CCL17', 'CCL18', 'CCL19', 'CCL2', 'CCL20', 'CCL21', 'CCL22', 'CCL23', 'CCL24', 'CCL25', 'CCL26', 'CCL3', 'CCL3L1', 'CCL3L3', 'CCL4', 'CCL5', 'CCL7', 'CCL8', 'CCR2', 'CCR6', 'CCR7', 'CD200', 'CD200R1', 'CD99', 'CD99L2', 'CH25H', 'CKLF', 'CORO1A', 'CRK', 'CRKL', 'CRTAM', 'CX3CL1', 'CXCL10', 'CXCL11', 'CXCL12', 'CXCL13', 'CXCL16', 'CYP7B1', 'DEFA1', 'DEFA1B', 'DOCK8', 'ECM1', 'EXT1', 'F11R', 'FADD', 'FUT4', 'FUT7', 'GAS6', 'GATA3', 'GBA1', 'GCSAM', 'GCSAML', 'GPR15', 'GPR15LG', 'GPR183', 'HSD3B7', 'ICAM1', 'IL27RA', 'ITGA4', 'ITGAL', 'ITGB3', 'ITGB7', 'JAM2', 'KLRC4-KLRK1', 'KLRK1', 'LRCH1', 'MADCAM1', 'MIA3', 'MSMP', 'MSN', 'MYO1G', 'NEDD9', 'OXSR1', 'PADI2', 'PIK3CD', 'PIK3CG', 'PLEC', 'PTK2B', 'PYCARD', 'RET', 'RHOA', 'RIPK3', 'RIPOR2', 'S100A7', 'S1PR1', 'SAA1', 'SELENOK', 'SLC12A2', 'SLC8B1', 'SPN', 'SPNS2', 'STK10', 'STK39', 'TBX21', 'TMEM102', 'TNFRSF14', 'TNFSF14', 'WASL', 'WNK1', 'WNT5A', 'XCL1', 'XCL2', 'XG', 'ZAP70')
  )
  geneSet = lapply(geneSet, function(gene){gene[gene %in% rownames(seuratObject)]})
  group = seuratObject@meta.data[["group"]] %>% unique()
  group = group[-which(group=="Ctrl")]
  markers = list()
  for(i in group){
    print(i)
    markers[[i]] = FindMarkers(seuratObject,
                               features=geneSet,
                               group.by="group",
                               ident.1=i,
                               ident.2="Ctrl",
                               min.pct=0,
                               logfc.threshold=0)
  }
  for(i in markers){
    print(i[, "avg_log2FC"])
  }
}

if(TRUE){
  seuratObject = scRNAMerged
  geneSet = lapply(geneSet, function(gene){gene[gene %in% rownames(seuratObject)]})
  geneSetTemp = unlist(geneSet) %>% unique()
  markers = FindAllMarkers(seuratObject, features=geneSetTemp, min.pct=0, logfc.threshold=0, only.pos=FALSE)
  #saveRDS(markers, "0004_markersBetweenCellTypeInAllCell.rds")
  #markers = readRDS("0004_markersBetweenCellTypeInAllCell.rds")
  # 将gene结果归类到geneSet中
  markersGeneSet = lapply(geneSet, function(gene){markers[markers[["gene"]] %in% gene, ]})
  
  
  cellType = "CD8+ T"
  gene = "CXCR4"
  tempBarcode = colnames(scRNAMerged)[scRNAMerged@meta.data$cellType==cellType]
  temp = scRNAMerged[, tempBarcode]
  group = scRNAMerged@meta.data[["group"]] %>% unique()
  test = list()
  for(i in group){
    test[[i]] = FindMarkers(temp,
                            features=gene,
                            group.by="group",
                            ident.1=i,
                            ident.2="Ctrl",
                            min.pct=0,
                            logfc.threshold=0,
                            only.pos=FALSE
    )
  }
  unlist(test)
}



if(TRUE){
  
  temp = scRNAMerged[, barcodeCellType]
  #temp@meta.data[["group_cellType"]] = paste0(temp@meta.data$group, '_', temp@meta.data$cellType)
  #Idents(temp) = temp@meta.data$group_cellType
  group = list()
  for(i in unique(scRNAMerged@meta.data$cellType)){
    print(i)
    barcodeCellType = rownames(scRNAMerged@meta.data)[scRNAMerged@meta.data$cellType==i]
    temp = scRNAMerged[, barcodeCellType]
    group[[i]] = FindMarkers(temp, assay="RNA", features=geneSet$lymphocyte_migration,
                             group.by="group", ident.1="Endometrioid", ident.2="Ctrl",
                             min.pct=0.25, logfc.threshold=0.25, only.pos=TRUE)
  }
  debug = FindMarkers(temp, assay="rPCA", features=geneSet$lymphocyte_migration,
                      group.by="group", ident.1="RPL", ident.2="Ctrl",
                      min.pct=0, logfc.threshold=0, only.pos=FALSE)
}

if(TRUE){
  scRNAMerged@meta.data[["cellTypeCluster"]] = scRNAMerged@meta.data[["cellType"]]
  for(cellType in levels(scRNAMerged@meta.data[["cellTypeCluster"]])){
    if(cellType %in% cellClustersClassified[["immune_cells"]]){
      levels(scRNAMerged@meta.data[["cellTypeCluster"]])[levels(scRNAMerged@meta.data[["cellTypeCluster"]])==cellType] = "immune_cells"
    }else{
      levels(scRNAMerged@meta.data[["cellTypeCluster"]])[levels(scRNAMerged@meta.data[["cellTypeCluster"]])==cellType] = "other_cells"
    }
  }
  seuratObject = subset(scRNAMerged, cellTypeCluster=="immune_cells")
  # 获取每个样本中每个基因的counts占该样本中所有基因的counts的比例
  sampleList = list()
  for(sample in unique(seuratObject@meta.data[["sample.x"]])){
    temp = list()
    for(cellType in unique(seuratObject@meta.data[["cellType"]])){
      barcode = rownames(seuratObject@meta.data)[seuratObject@meta.data[["sample.x"]]==sample & seuratObject@meta.data[["cellType"]]==cellType]
      if(length(barcode)){temp[[cellType]] = barcode}
    }
    sampleList[[sample]] = temp
  }
  #temp = seuratObject@assays$RNA@counts
  Sys.time()
  sampleList = lapply(sampleList, function(barcodeList){
    barcode = unlist(barcodeList) %>% unique()
    temp2 = seuratObject@assays$RNA@counts[, barcode]
    rateList = lapply(barcodeList, function(barcode){
      data = temp2[, barcode, drop=FALSE]
      totalSum = sum(data)
      # 筛选出在该样本中10%及以上的细胞都表达的gene
      cutoffValue = dim(data)[2] * 0.1  # 10%的细胞
      gene = data >= 1
      gene = rowSums(gene)
      gene = gene[gene>=cutoffValue]
      gene = names(gene)
      # 计算gene的counts占该样本中所有基因的counts的比例
      data = data[gene, , drop=FALSE]
      geneSum = rowSums(data)
      geneRate = (geneSum / totalSum) * 10^6
      return(geneRate)
    })
    return(rateList)
  })
  Sys.time()
  
  group1 = "RPL"
  group2 = "Ctrl"
  cellType = "Bcells"
  pCutoff = 0.05
  
  sampleGroup1 = seuratObject@meta.data[seuratObject@meta.data[["group"]]==group1, "sample.x"] %>% unique()
  sampleGroup2 = seuratObject@meta.data[seuratObject@meta.data[["group"]]==group2, "sample.x"] %>% unique()
  sampleGroup1 = lapply(sampleList[sampleGroup1], function(cellTypeList){return(cellTypeList[[cellType]])})
  sampleGroup2 = lapply(sampleList[sampleGroup2], function(cellTypeList){return(cellTypeList[[cellType]])})
  sampleGroup1 = Filter(function(x){!is.null(x)}, sampleGroup1)
  sampleGroup2 = Filter(function(x){!is.null(x)}, sampleGroup2)
  geneGroup1 = lapply(sampleGroup1, function(geneExpression){names(geneExpression)}) %>%
                 unlist() %>%
                 unique()
  geneGroup2 = lapply(sampleGroup2, function(geneExpression){names(geneExpression)}) %>%
                 unlist() %>%
                 unique()
  geneList = list()
  for(gene in unique(c(geneGroup1, geneGroup2))){
    geneList[[gene]] = gene
  }
  geneList = future_lapply(geneList, function(gene){
    expressionGroup1 = lapply(sampleGroup1, function(geneExpression){
      if(gene %in% names(geneExpression)){return(geneExpression[gene])}else{return(0)}
    }) %>% unlist()
    expressionGroup2 = lapply(sampleGroup2, function(geneExpression){
      if(gene %in% names(geneExpression)){return(geneExpression[gene])}else{return(0)}
    }) %>% unlist()
    wilcoxResult = wilcox.test(expressionGroup1, expressionGroup2)
    meanExpressionGroup1 = mean(expressionGroup1)
    meanExpressionGroup2 = mean(expressionGroup2)
    log2FC = log2(meanExpressionGroup1/meanExpressionGroup2)
    
    result = list(wilcoxResult=wilcoxResult, log2FC=log2FC)
    return(result)
  })
  df = lapply(geneList, function(d){
    result = c(p=d$wilcoxResult$p.value,
                  w=as.numeric(d$wilcoxResult$statistic),
                  log2FC=d$log2FC)
    return(result)
  }) %>% as.data.frame() %>% t() %>% as.data.frame()
  df = df[df['p']<=pCutoff, ]
  df = df[order(df['p']), ]
  
  
  
  
}

# 免疫细胞
if(TRUE){
  scRNAMerged@meta.data[["cellTypeCluster"]] = scRNAMerged@meta.data[["cellType"]]
  for(cellType in levels(scRNAMerged@meta.data[["cellTypeCluster"]])){
    if(cellType %in% cellClustersClassified[["immune_cells"]]){
      levels(scRNAMerged@meta.data[["cellTypeCluster"]])[levels(scRNAMerged@meta.data[["cellTypeCluster"]])==cellType] = "immune_cells"
    }else{
      levels(scRNAMerged@meta.data[["cellTypeCluster"]])[levels(scRNAMerged@meta.data[["cellTypeCluster"]])==cellType] = "other_cells"
    }
  }
  seuratObject = subset(scRNAMerged, cellTypeCluster=="immune_cells")
  sampleList = list()
  for(sample in unique(seuratObject@meta.data[["sample.x"]])){
    barcode = rownames(seuratObject@meta.data)[seuratObject@meta.data[["sample.x"]]==sample]
    if(length(barcode)){
      sampleList[[sample]] = barcode
    }
  }
  sampleList = lapply(sampleList, function(barcode){
    data = seuratObject@assays$RNA@counts[, barcode]
    totalSum = sum(data)
    # 筛选出在该样本中10%及以上的细胞都表达的gene
    if(TRUE){
      cutoffValue = dim(data)[2] * 0.1  # 10%的细胞
      gene = data >= 1
      gene = rowSums(gene)
      gene = gene[gene>=cutoffValue]
      gene = names(gene)
      data = data[gene, , drop=FALSE]
    }
    # 计算gene的counts占该样本中所有基因的counts的比例
    geneSum = rowSums(data)
    geneRate = (geneSum / totalSum) * 10^6

    return(geneRate)
  })

  df = data.frame("log2FC"=c(), "group1"=c(), "gene"=c())
  for(group in unique(seuratObject@meta.data[["group"]])){
    cat(group, ", ")
    if(group=="Ctrl"){next}
    dfTemp = compareLog2FC(seuratObject=seuratObject,
                           sampleList=sampleList,
                           group1=group,
                           group2="Ctrl",
                           pCutoff=0.05)
    dfTemp = dfTemp[, c("log2FC", "group1")]
    dfTemp[["gene"]] = rownames(dfTemp)
    df = rbind(df, dfTemp)
  }
  df = spread(df, key="group1", value="log2FC")
  rownames(df) = df[["gene"]]
  df = df[, c(colnames(df)[-which(colnames(df)=="gene")])]
  saveRDS(df, "0004_tempExpressionDifGroup.rds")
  
  gene = c(geneSet$refLymphocyteRecruitment, geneSet$refChemokines)
  gene = geneSet$cell_chemotaxis
  gene = gene[gene %in% rownames(df)]
}

# 绘图
if(TRUE){
  dfTemp = df[gene, ]
  #dfTemp = dfTemp[rowSums(is.na(dfTemp))<=3, ]  # 筛选出缺失值少于指定个数的行
  dfTemp = as.matrix(dfTemp)
  dfTemp[is.na(dfTemp)] = 0
  dfTemp[is.infinite(dfTemp) & dfTemp>0] = max(dfTemp[!is.infinite(dfTemp)])
  dfTemp[is.infinite(dfTemp) & dfTemp<0] = min(dfTemp[!is.infinite(dfTemp)])
  
  dfTemp = as.data.frame(dfTemp)
  dfTemp = dfTemp[order(5*dfTemp$`Asherman’s Syndrome`+4*dfTemp$Endometrioid+3*dfTemp$endometriosis+2*dfTemp$RPL+dfTemp$Thin), ]
  
  colorBarRange = abs(min(dfTemp)) + abs(max(dfTemp))
  colorBar0Rate = abs(min(dfTemp)) / colorBarRange
  plotGeneProportion <- plot_ly(x=colnames(dfTemp),
                                y=rownames(dfTemp),
                                z=as.matrix(dfTemp),
                                type="heatmap",
                                width=500, height=500,
                                colorscale=data.frame(c(0, colorBar0Rate, 1),
                                                      c("rgb(49,54,149)", "lightgrey", "rgb(165,0,38)")),
                                colorbar=list(title="log2FC", titleside="left",
                                              x=-0.2, xref="container",
                                              y=0.5, yref="container",
                                              ticks="inside")
                                )
  plotGeneProportion = plotGeneProportion %>%
    add_annotations(text="immunes cells (including Macrophage, NK,<br>monocyte, T cells, B cells, DC, mast cells, and ILC3)",
                    x=-1, y=0.5,
                    xref="x", yref="paper",
                    textangle=270,
                    font=list(family="Arial", size=12, color="black"),
                    showarrow=FALSE)
  plotGeneProportion = plotGeneProportion %>%
    layout(margin=list('l'=50, 'r'=50, 't'=50, 'b'=50),
           font=list(family="Arial",
                     color="black",
                     size=12),
           xaxis=list(side="top",
                      tickangle=270,
                      ticklen=0,
                      showgrid=FALSE,
                      zeroline=0),
           yaxis=list(side="right",
                      title="",
                      ticklen=0,
                      showgrid=FALSE,
                      zeroline=0),
           legend=list(x=0, y=0)
    )
  plotGeneProportion
}

if(FALSE){
  # del
  col = "severe"
  debug$temp = ifelse(debug$disease==col, debug$disease, NA)
  debug$show = ifelse(is.na(debug$temp), debug$show, debug$temp)
  debug[debug$disease=="potential_endometrial_disorders", "show"] = "ctrl"
  tempMetaData = seuratObject@meta.data
  seuratObject@meta.data = debug
  pdf("temp.pdf", height=15, width=15)
  DimPlot(seuratObject, split.by="show", group.by="cellType", ncol=3, cols=c("#3de1ad","#83ccd2","#ff4777","#a98175","#1685a9","#426666","#e4c6d0",
                                                                             "#afdd22","#c97586","#4c221b","#884898","#9d2933","#424c50","#ffec47"))
  dev.off()
}

# 将seuratB/seuratT的cellTypeSub添加到scRNAMerged中
if(TRUE){
  # 将seuratB/seuratT的cellTypeSub添加到scRNAMerged中
  seuratB@meta.data$cellTypeSub = as.character(seuratB@active.ident)
  metaInfo = seuratB@meta.data
  metaInfo$barcode = rownames(metaInfo)
  metaInfo = metaInfo[, c("barcode", "cellTypeSub")]
  
  metaOld = scRNAMerged@meta.data
  metaOld$barcode = rownames(metaOld)
  metaOld = metaOld[, -which(colnames(metaOld)=="cellTypeSub")]
  
  metaNew = dplyr::left_join(metaOld, metaInfo, by="barcode")
  
  seuratT@meta.data$cellTypeSubTemp = as.character(seuratT@active.ident)
  metaInfo = seuratT@meta.data
  metaInfo$barcode = rownames(metaInfo)
  metaInfo = metaInfo[, c("barcode", "cellTypeSubTemp")]
  
  metaNew = dplyr::left_join(metaNew, metaInfo, by="barcode")
  if(FALSE %in% is.na(metaNew[!is.na(metaNew$cellTypeSubTemp), "cellTypeSub"])){warning("[Error]Conflict existed in cellTypeSub and cellTypeSubTemp")}
  metaNew$cellTypeSub = ifelse(is.na(metaNew$cellTypeSub), metaNew$cellTypeSubTemp, metaNew$cellTypeSub)
  
  rownames(metaNew) = metaNew$barcode
  metaNew = metaNew[, -which(colnames(metaNew)=="cellTypeSubTemp")]
  metaNew = metaNew[, -which(colnames(metaNew)=="barcode")]
}

if(TRUE){
  # 绘制F2h
  # 数据处理
  if(TRUE){
    geneSetRef = lapply(geneSetRef, function(geneSet){geneSet[geneSet %in% rownames(scRNAMerged)]})
    
    # 获取不同groupSub的barcode
    #scRNAMerged = AddModuleScore(scRNAMerged, features=geneSetRef, name=names(geneSetRef))
    barcodeList = list()
    for(g in unique(scRNAMerged@meta.data[["group"]])){
      barcodeList[[g]] = rownames(scRNAMerged@meta.data)[scRNAMerged@meta.data[["group"]]==g]
    }
    barcodeList[["Ctrl"]] = rownames(scRNAMerged@meta.data)[(scRNAMerged@meta.data$group=="Ctrl")]
    barcodeList = barcodeList[c("EMs(Endometrioma)", "EC(EC)", "AS(Severe)")]
    # 筛选出免疫细胞
    barcodeListImmune = lapply(barcodeList, function(barcode){
      barcode[scRNAMerged@meta.data[barcode, "cellType"] %in% cellClustersClassified[["immune_cells"]]]
    })
    
    # 根据barcode获取基因集打分
    geneSetRefScore = list()
    for(g in names(geneSetRef)){
      geneSetRefScore[[g]] = lapply(barcodeListImmune, function(barcode){scRNAMerged@meta.data[barcode, g]})
    }
    # 对基因集打分取中位数
    geneSetRefScore = lapply(geneSetRefScore, function(geneSet){
      lapply(geneSet, function(group){median(group)})
    })
    
    # 获取不同group中淋巴细胞的占比
    lymphocyteRate = lapply(barcodeList, function(barcode){
      lymphocyteNum = ifelse(scRNAMerged@meta.data[barcode, "cellType"] %in% cellTypeClassified[["Lymphocyte"]], TRUE, FALSE) %>% sum()
      rate = lymphocyteNum / length(barcode)
    }
    ) %>% unlist()
    lymphocyteRate = lymphocyteRate[order(names(lymphocyteRate))]
    
    correlationList = lapply(geneSetRefScore, function(geneScore){
      geneScore = geneScore %>% unlist()
      geneScore = geneScore[order(names(geneScore))]
      cor.test(geneScore, lymphocyteRate, method="pearson")
    })
  }
  # 绘制图片
  if(TRUE){
    colorsTemp = list("AS(Severe)"="blue", "EC(EC)"="red", "EMs(Endometrioma)"="green")
    plot = list()
    for(geneSetName in names(geneSetRefScore)){
      geneScore = geneSetRefScore[[geneSetName]] %>% unlist()
      geneScore = geneScore[order(names(geneScore))]
      lymphocyteRate = lymphocyteRate
      plot[[geneSetName]] = plot_ly()
      for(disease in names(lymphocyteRate)){
        plot[[geneSetName]] = plot[[geneSetName]] %>% add_trace(x=lymphocyteRate[disease],
                                                                y=geneScore[disease],
                                                                type="scatter",
                                                                mode="markers",
                                                                name=disease,
                                                                color=colorsTemp[[disease]]
                                                                )
      }
      
      plot[[geneSetName]] = plot[[geneSetName]] %>% layout(font=list(family="Arial",
                                                                     color="black",
                                                                     size=12),
                                                           xaxis=list(categoryorder="trace",
                                                                      tickangle="270",
                                                                      showline=TRUE,
                                                                      showticklabels=TRUE,
                                                                      ticklen=5,
                                                                      zeroline=FALSE,
                                                                      showgrid=FALSE),
                                                           yaxis=list(title='Percent(%)',
                                                                      showline=TRUE,
                                                                      showticklabels=TRUE,
                                                                      ticklen=5,
                                                                      zeroline=FALSE,
                                                                      showgrid=FALSE),
                                                           paper_bgcolor="white",
                                                           plot_bgcolor="white")
    }

    fig = subplot(plot[[1]], plot[[2]], plot[[3]],
                    nrows=1)
    fig = fig %>% layout()
  }
}



# test
if(FALSE){
  debug <- data.frame(A = c(1, 2, 3),
                     B = c(5, 10, 15))
  
  scaled_data <- scale(debug, center=FALSE, scale=FALSE)
  
}








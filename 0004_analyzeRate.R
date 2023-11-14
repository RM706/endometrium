.libPaths(c("/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/",
            "/home/weiyihu/App/conda/envs/RStudio/lib/R/library"
            #"/home/weiyihu/package/share/4.0"
            ))
setwd("/data1/weiyihu/endometrium/data")
library("plotly", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("future", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("future.apply", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
library("tidyverse")
library("Seurat", lib.loc="/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/")
library("qpdf", lib.loc="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")

# 前置参数
if(TRUE){
  MAXCORES = 4
  inputFileName = "/data1/shumeng/5.NK_cell/6.TLOs/1.project/All_Tissue_RPCA_sample.rds"
  outputFileName = "0004_plot.pdf"
  colors=list("Stromal"="#3de1ad",             "Unciliated Epithelial"="#83ccd2",
              "NK"="#ff4777",                  "Macrophages"="#a98175",
              "CD8+ T"="#1685a9",              "Endothelial"="#426666",
              "CD4+ T"="#e4c6d0",              "Proliferative cells"="#ffa631",
              "NKT"="#afdd22",                 "Monocytes"="#c97586",
              "Ciliated Epithelial"="#4c221b", "Bcells"="#884898",
              "ILC3"="#9d2933",                "DC"="#424c50",
              "Mast cells"="#ffec47")
  cellTypeClassified = list("Stromal_cells"=c("Stromal", "Unciliated Epithelial",
                                              "Endothelial", #"Proliferative cells",
                                              "Ciliated Epithelial"),
                            "Myeloid_cells"=c("Macrophages", "Monocytes",
                                              "DC", "Mast cells"),
                            "Lymphocyte"=c("NK", "CD8+ T",
                                           "CD4+ T", "NKT",
                                           "Bcells", "ILC3"))
  cellClustersClassified = list("immune_cells"=c(cellTypeClassified[["Myeloid_cells"]],
                                                 cellTypeClassified[["Lymphocyte"]]),
                                "stromal_cells"=cellTypeClassified[["Stromal_cells"]]
                                )
}

# 初始化
if(TRUE){
  plan(multisession, workers=MAXCORES)
  options(future.globals.maxSize=10*1024*1024*1024,  # 10GB
          scipen=10)  
}



# main
if(TRUE){
  # 加载数据
  if(TRUE){
    scRNAMerged = readRDS(inputFileName)
  }
  # meta.data中添加cellType信息
  scRNAMerged@meta.data[["cellType"]] = scRNAMerged@active.ident
  
  # 统计在不同类型疾病中的不同类型细胞的比例
  if(TRUE){
    # 准备数据
    if(TRUE){
      df = scRNAMerged@meta.data
      df = df[df[["type"]]!="CD45+", ]
      df = table(df[["cellType"]], df[["group"]]) %>%
        as.data.frame() %>%
        spread(., key="Var2", value="Freq")
      rownames(df) = df[["Var1"]]
      df = df[,-1]
      
      dfCellType = df %>% mutate_at(.vars=vars(everything()), list(function(x){x/sum(x)}))
      dfCellType = dfCellType %>% select(order(dfCellType["Stromal", ], decreasing=TRUE))
      #dfCellType = dfCellType %>% select("Ctrl", everything())
      dfCellType = dfCellType %>% select("Ctrl", "Thin", "Asherman’s Syndrome", "endometriosis","Endometrioid", "RPL")
      dfCellType = dfCellType * 100
    }
    # 绘图
    if(TRUE){
      plotCellType <- plot_ly(x=colnames(dfCellType),
                              y=as.numeric(dfCellType[1, ]),
                              type="bar",
                              name=rownames(dfCellType)[1],
                              marker=list(color=colors[[rownames(dfCellType)[1]]]),
                              width=400, height=500
                              )
      for(i in seq(2, dim(dfCellType)[1])){
        plotCellType = plotCellType %>% add_trace(y=as.numeric(dfCellType[i, ]),
                                                  name=rownames(dfCellType)[i],
                                                  marker=list(color=colors[[rownames(dfCellType)[i]]]))
      }
      plotCellType <- plotCellType %>% layout(barmode='stack',
                                              xaxis=list(categoryorder="trace",
                                                         tickangle="270",
                                                         showline=TRUE,
                                                         showticklabels=TRUE,
                                                         ticklen=5),
                                              yaxis=list(title='Percentage(%)',
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
  # 保存图片
  if(TRUE){
    reticulate::py_run_string("import sys")
    save_image(plotCellType, file="temp1.pdf")
    save_image(plotClusters, file="temp2.pdf")
    save_image(plotLine, file="temp3.pdf")
    
    
    qpdf::pdf_combine(c("temp1.pdf", "temp2.pdf", "temp3.pdf"),
                      output=outputFileName)
    file.remove(c("temp1.pdf", "temp2.pdf", "temp3.pdf"))
  }
}



# debug



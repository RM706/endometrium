# 函数
# getCellTypeColor    获取细胞类型所对应的颜色
# makeHeatmap         绘制在不同group中gene表达的热图

# 一级函数
getCellTypeColor = function(cellType, colorList){
  # input:
  #   cellType|\t   factor|\t     细胞类型
  #   colorList|\t  list|\t       key为细胞类型, value为颜色
  # return:
  #   result|\t     character|\t  细胞类型的levels对应的颜色
  cellType = cellType
  colorList = colorList
  
  result = c()
  
  # 数据类型检查
  dataType = class(cellType)
  if(dataType!="factor"){stop(paste0("The data type of parament cellType should be factor, it is ", dataType))}
  
  for(ct in levels(cellType)){
    if(ct %in% names(colorList)){
      result = append(result, colorList[[ct]])
    }else{
      stop(paste0(ct, " not recorded in the colorList"))
    }
  }
  
  return(result)
}

makeHeatmap = function(df, rowSet=NULL, rowColor=NULL, colColor=NULL, rowLevels=NULL, colLevels=NULL){
  # input:
  #   df
  #   rowSet|\t     list|\t             key为row的分组名, value为对应分组中的rowName
  #   rowColor|\t   list|\t             key为row的组名, value为该组的颜色
  #   colColor|\t   list|\t             key为col的组名, value为该组的颜色
  #   rowLevels|\t  characterVector|\t  =names(rowSet), 行组的顺序
  #   colLevels|\t  characterVector|\t  =unique(colnames(df)), 列组的顺序
  #library(ComplexHeatmap)
  #library(circlize)
  #library(grDevices)
  df = df
  rowSet = rowSet
  rowColor = rowColor
  colColor = colColor
  rowLevels = rowLevels
  colLevels = colLevels
  
  # 初始化
  rowSet = ifelse(is.null(rowSet), list(list("row"=rownames(df))), list(rowSet))[[1]]
  rowColor = ifelse(is.null(rowColor), list(stats::setNames(grDevices::rainbow(length(rowSet)), names(rowSet))), list(unlist(rowColor)[names(rowColor) %in% names(rowSet)]))[[1]]
  colColor = ifelse(is.null(colColor), list(stats::setNames(grDevices::rainbow(dim(df)[2]), colnames(df))), list(unlist(colColor)[names(colColor) %in% colnames(df)]))[[1]]
  rowLevels = ifelse(is.null(rowLevels), list(names(rowSet)), list(rowLevels[rowLevels %in% names(rowSet)]))[[1]]
  colLevels = ifelse(is.null(colLevels), list(unique(colnames(df))), list(colLevels[colLevels %in% colnames(df)]))[[1]]

  df = as.data.frame(t(df))
  
  #顶部注释
  dfAnno = df
  dfAnno$column = rownames(dfAnno)
  dfAnno <- data.frame(dfAnno$column)
  dfAnno[, 1] = factor(dfAnno[, 1], levels=colLevels)
  colnames(dfAnno) <- "group"
  top_anno = HeatmapAnnotation(df=dfAnno,
                               border=F,
                               show_annotation_name=F,
                               gp=gpar(col=NA),
                               col=list("group"=colColor),
                               show_legend=T,
                               which='column')
  
  # 左侧注释
  dfAnno = df
  dfAnno = as.data.frame(t(dfAnno))
  dfAnno$rowName = rownames(dfAnno)
  dfAnno <- data.frame(dfAnno$rowName)
  dfAnno$group = NA
  for(i in seq(1, dim(dfAnno)[1])){
    rowName = dfAnno[i, 1]
    for(group in names(rowSet)){
      if(rowName %in% rowSet[[group]]){dfAnno[i, 2] = group}
    }
  }
  dfAnno = dfAnno[, 2, drop=FALSE]
  colnames(dfAnno) <- "gene set"
  dfAnno[, 1] = factor(dfAnno[, 1], levels=rowLevels)
  split = dfAnno[, 1]
  left_anno = HeatmapAnnotation(df=dfAnno,
                                border=F,
                                show_annotation_name=F,
                                gp=gpar(col=NA),
                                col=list("gene set"=unlist(rowColor)),
                                show_legend=T,
                                which='row')
  
  # 数据标准化缩放
  marker_exp <- t(scale(df[, ], scale=T, center=T))
  cat(paste0("min(int): ", floor(min(marker_exp)), ' max(int): ', ceiling(max(marker_exp)), '\n'))
  
  #颜色设置
  col_fun <- colorRamp2(c(floor(min(marker_exp)), 0, ceiling(max(marker_exp))),
                        c("#313695", "white", "#A50026"))
  col_fun(seq(floor(min(marker_exp)), ceiling(max(marker_exp))))
  
  # 主图
  plot <- Heatmap(marker_exp,
                  row_split=split,
                  cluster_row_slices=FALSE,
                  cluster_rows=T,  # 是否进行行聚类
                  show_row_dend=FALSE,
                  cluster_columns=F,  # 是否进行列聚类
                  row_title=NULL,
                  #row_title=gt_render(c("", "", "Immune cells", "", "", "GC B cells"), r=unit(2, "pt"), padding=unit(c(2,2,2,2), "pt")),
                  show_row_names=T,  # 是否显示行名
                  row_names_side='right',  # 行名的位置
                  row_names_gp=gpar(fontsize=10),  # 行名的字体
                  show_column_names=F,  # 是否显示列名
                  column_title=NULL,  # 所有行的title
                  column_names_gp=gpar(fontsize=10),  # 列名的字体
                  heatmap_legend_param=list(title=' '),
                  col=col_fun,  # 色条的设置
                  border=NULL,  # 热图最外侧边界的颜色
                  top_annotation=top_anno,
                  left_annotation=left_anno)
  
  return(plot)
}

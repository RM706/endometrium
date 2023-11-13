.libPaths("/home/weiyihu/package/share/4.0")
setwd("/data1/weiyihu/endometrium/data")
library(Seurat, lib.loc="/home/weiyihu/package/share/4.0")
library(dplyr, lib.loc="/home/weiyihu/package/share/4.0")
library(patchwork, lib.loc="/home/weiyihu/package/share/4.0")
library(Matrix, lib.loc="/home/weiyihu/package/share/4.0")
library(future)
library(future.apply)


MAXCORES = 6
scRNA_path = "/data1/weiyihu/endometrium/data/raw_data"
share_path = "/data1/weiyihu/endometrium/data/share"
studyNameFile3 = c("GSE164449")
studyNameH5 = c("1.CRA002181", "3.HRA000237")
studyNameRData = c("GSE130560")
geneList = list(cellT = c("CD3D", "CD8A", "CD4"),
                #cellTref = c("FOXP3"),
                cellB = c("CD79A","CD20"),
                cellNK = c("NKG7", "KLRC1", "KLRD1", "NCAM1"),
                cellNKT = c("FGFBP2"),
                #cellNeutrophils = c("CST3", "LYZ", "FCGR3B", "CSF3R"),
                cellMacrophage = c("CD14", "FOLR2", "CD163"),
                cellMonocytes = c("FCN1", "VCAN", "CD14", "CD16", "S100A12", "FCGR3A"),
                cellDC = c("LILRA4", "CLEC9A", "CD1C"),
                cellMast = c("TPSAB1", "MS4A2", "KIT")
                #cellILC3 = c("NCR2")
)
markersRef = list(
  MSC=toupper(c("Dcn", "Col3a1", "Sparc", "Gsn", "Col1a2", "Gpx3", "Serping1", "Egr1", "Inmt", "Mfap4", "Igfbp7,Dkk3,Tbx20", "S100a6", "Sult1e1", "Junb", "Ifitm3", "Efemp1", "Fth1", "Igfbp6", "Rnase4", "Id3", "Penk", "Timp2")),
  SMC=toupper(c("Acta2", "Tpm2", "Myh11", "Tagln", "Ppp1r14a", "Prrx2", "Prss12", "Gxylt2", "Ptprz1", "Tpm1", "Btg2,Sept6", "Itga1", "Rbpms", "Dusp1", "Cpe", "Dsp", "Dhrs3", "Nr1d1", "Ldhb", "Rock2", "Vim", "Pls3", "Pttg1,Stk11", "Oat", "Tgm2", "Dpysl3")),
  Macrophage=toupper(c("Csf1r", "Apoe", "C1qa", "C1qb", "C1qc", "Ccl6", "Cd5l", "Cd68", "Cebpb", "Cfp", "Ctsc", "Ctss", "Cybb", "Fcer1g,Gpx1", "Hmox1", "Lgals3", "Lst1", "Ly86", "Lyz2", "Ms4a6c", "Pld4", "Ccl9", "Cybb", "Ms4a4c", "Psap,Plac8", "Ifitm3")),
  NK=toupper(c("Gzmb", "Klra8 ,Klre1", "Klra7 ", "Klrk1", "Klrb1c", "Klrd1", "Klra4", "Serpinb6b", "Il2rb", "Ctsw", "Gzma,Nkg7", "Xcl1", "Ccl5", "Ms4a4b", "Id2", "Ccl4")),
  EC=toupper(c("Igfbp7", "Aqp1", "Ptprb", "Gpihbp1", "Clec14a", "Clec4g", "Kdr", "Ehd3", "Eng", "Egfl7", "Adgrf5", "Timp3,Bgn", "Bmp2", "Il6st", "Tinagl1", "Sdpr", "Fabp4", "Rspo3", "Fcgr2b", "Dnase1l3", "Adgrf5", "Aqp1", "Cldn5")),
  Monocyte=toupper(c("S100a9", "S100a8", "Lcn2", "G0s2", "Il1b", "Ccrl2", "Cxcl2", "Hdc", "Mxd1", "Irg1", "Hcar2", "Gm5483,Neat1", "Slfn4", "Nlrp3", "Slc7a11", "Csf3r", "Clec4e", "Stfa2l1", "Tnfaip2", "Il1f9", "Asprv1")),
  Tcell=toupper(c("Cd3d", "Cd3e", "Cd3g", "Cd7", "Cd8a", "Ccl5", "Cd28", "Ctsw", "Cxcr6", "Gimap3", "Gimap7", "Gzmb", "Hcst,Il2rb", "Il7r")),
  B=toupper(c("Cd79a", "Cd19", "Cd79b", "Bank1", "Bcl11a", "Ms4a1", "Ebf1", "Fcmr", "Mzb1", "Ly6d", "Siglecg", "H2DMb2", "H2-Ob", "Ltb", "Cd2", "Gm43603", "Ralgps2", "Ptprcap", "Cd37", "Napsa", "Cd24a", "Foxp1")),
  DC=toupper(c("Xcr1", "Cd24a", "Ifi205", "Flt3", "Tnni2", "Amica1", "Naaa", "Septin5", "Cst3", "Plbd1", "Tmsb10", "Irf8,Lsp1", "Fscn1", "Ccr7", "Ccl22", "Tbc1d4", "Il12b", "Il4i1", "Traf1", "Eno3")),
  pDC=toupper(c("Cox6a2", "Klk1", "Siglech", "Rnase6", "D13Ertd608e", "Ccr9", "Klk1b27", "Cd7", "Atp1b1", "Smim5,Upb1", "Cd209d", "Pacsin1", "Ppfia4", "Sh3bgr", "P2ry14", "Cadm1", "Gm5547", "Cd209a")),
  Neuron=toupper(c("Plp1", "Kcna1", "fra3", "Egfl8", "Sostdc1", "Mbp", "Gpr37l1", "Cadm4"))
)

# 一级函数
tempRead10X = function(data.dir, barcodesFile, featuresFile, matrixFile,
                       gene.column=2, cell.column=1, unique.features=TRUE, 
                       strip.suffix=FALSE){
  full.data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, barcodesFile)
    gene.loc <- file.path(run, "genes.tsv")
    features.loc <- file.path(run, featuresFile)
    matrix.loc <- file.path(run, matrixFile)
    pre_ver_3 <- file.exists(gene.loc)
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing. Expecting ", 
           basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", 
           basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    if (has_dt) {
      cell.barcodes <- as.data.frame(data.table::fread(barcode.loc, 
                                                       header = FALSE))
    }
    else {
      cell.barcodes <- read.table(file = barcode.loc, header = FALSE, 
                                  sep = "\t", row.names = NULL)
    }
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    }
    else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    if (has_dt) {
      feature.names <- as.data.frame(data.table::fread(ifelse(test = pre_ver_3, 
                                                              yes = gene.loc, no = features.loc), header = FALSE))
    }
    else {
      feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                                yes = gene.loc, no = features.loc), header = FALSE, 
                                  stringsAsFactors = FALSE)
    }
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as.sparse(x = list_of_data[[j]])
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}


# 一级函数
makeData = function(dataPath, studyName, mode, data=NULL){
  dataPath = dataPath
  studyName = studyName
  mode = mode
  
  if(is.null(data)){data=list()}
  
  # 获取每个研究中不同样本的路径
  if(mode == "3"){
    for(study in studyName){
      dir = file.path(dataPath, study)
      sampleName = list.files(dir)
      for(sample in sampleName){
        path = file.path(dir, sample)
        data[[sample]] = list(study=study, sampleName=sample, path=path)
      }
    }
  }else if(mode == "h5"){
    for(study in studyName){
      dir = file.path(dataPath, study)
      sampleName = basename(list.dirs(dir, recursive=FALSE))
      for(sample in sampleName){
        path = file.path(dir, sample, "outs", "filtered_feature_bc_matrix")
        data[[sample]] = list(study=study, sampleName=sample, path=path)
      }
    }
    
  }

  return(data)
}


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


# 二级函数
loadFile = function(data){
  # input:
  #   data, list
  # change:
  #   读取单细胞数据
  # output:
  #   list
  
  # 进度条
  pb = txtProgressBar(min=1, max=length(data), style=3)
  i=1
  # 读取单细胞数据
  for(sample in names(data)){
    i = i+1
    setTxtProgressBar(pb, i)
    files = list.files(data[[sample]][["path"]])
    barcodesFile = files[grepl(".*barcodes.tsv.*", files)]
    featuresFile = files[grepl(".*features.tsv.*", files)]
    matrixFile = files[grepl(".*matrix.mtx.*", files)]
    dataTemp = tempRead10X(data[[sample]][["path"]],
                           barcodesFile=barcodesFile,
                           featuresFile=featuresFile,
                           matrixFile=matrixFile)
    data[[sample]][["data"]] = CreateSeuratObject(counts=dataTemp,
                                                  project=data[[sample]][["sampleName"]])
  }
  close(pb)
  
  return(data)
}


# 前置参数
if(TRUE){
  plan(multisession, workers=MAXCORES)
  options(future.globals.maxSize=10*1024*1024*1024)  # 10GB
}


# 第一部分, 整理样本
if(TRUE){
  data = list()
  # 准备3文件单细胞数据
  data = makeData(dataPath=scRNA_path,
                  studyName=studyNameFile3,
                  mode="3",
                  data=data)
  data = makeData(dataPath=share_path,
                  studyName=studyNameH5,
                  mode="h5",
                  data=data)
  # 手动更改样本名
  data[["GSM5009716"]][["sampleName"]] = "GSE164449_Ctrl1"
  data[["GSM5009717"]][["sampleName"]] = "GSE164449_RPL1"
  data[["Ctrl_11"]][["sampleName"]] = "CRA002181_Ctrl1"
  data[["Ctrl_12"]][["sampleName"]] = "CRA002181_Ctrl2"
  data[["Ctrl_13"]][["sampleName"]] = "CRA002181_Ctrl3"
  data[["RPL_11"]][["sampleName"]] = "CRA002181_RPL1"
  data[["RPL_12"]][["sampleName"]] = "CRA002181_RPL2"
  data[["ND_1"]][["sampleName"]] = "HRA000237_Ctrl1"
  data[["ND_2"]][["sampleName"]] = "HRA000237_Ctrl2"
  data[["ND_3"]][["sampleName"]] = "HRA000237_Ctrl3"
  data[["RPL_1"]][["sampleName"]] = "HRA000237_RPL1"
  data[["RPL_2"]][["sampleName"]] = "HRA000237_RPL2"
  data[["RPL_3"]][["sampleName"]] = "HRA000237_RPL3"
  # 读取3文件单细胞数据
  data = loadFile(data=data)

  
  # 处理特殊格式的数据
  ## 处理GSE130560
  temp = load("/data1/weiyihu/endometrium/data/raw_data/GSE130560/GSE130560_matrix.RData")
  temp2 = read.csv("/data1/weiyihu/endometrium/data/raw_data/GSE130560/GSE130560_phenotype.csv")
  rownames(temp2) = temp2[["X"]]
  temp2 = temp2[,-1]
  temp = CreateSeuratObject(counts=matrix,
                            project="GSE130560",
                            meta.data=temp2)
  temp@meta.data$orig.ident = gsub("ND_", "GSE130560_Ctrl", temp@meta.data$orig.ident)
  temp@meta.data$orig.ident = gsub("PD_", "GSE130560_RSA", temp@meta.data$orig.ident)
  temp@meta.data$sample = gsub("ND_", "GSE130560_Ctrl", temp@meta.data$sample)
  temp@meta.data$sample = gsub("PD_", "GSE130560_RSA", temp@meta.data$sample)
  ## 拆分GSE130560
  temp = SplitObject(temp, split.by="sample")
  for(sampleName in names(temp)){
    seuratObject = temp[[sampleName]]
    seuratObject@project.name = sampleName
    Idents(seuratObject) = sampleName
    temp[[sampleName]] = seuratObject
  }
  ## 将拆分后的数据添加到data中
  for(sampleName in names(temp)){
    seuratObject = temp[[sampleName]]
    data[[sampleName]] = list(study="GSE130560", sampleName=sampleName, data=seuratObject)
  }
  
  
  # 提取多个seuratObject至一个list中
  scRNA = list()
  for(sampleName in names(data)){
    seuratObject = data[[sampleName]][["data"]]
    sampleName = data[[sampleName]][["sampleName"]]
    scRNA[[sampleName]] = seuratObject
  }
  
  # 准备质控
  scRNA = lapply(scRNA, function(seuratObject){
    seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern="^MT-")
    seuratObject[["percent.rb"]] <- PercentageFeatureSet(seuratObject, pattern="^RP[SL]")
    return(seuratObject)
  })
  ## 备份质控前的数据
  if(TRUE){
    scRNAMergedBeforeQC = merge(scRNA[[1]], scRNA[-1])  # 仅合并但不处理批次效应, 用于对比
    pdf("0001_mergedBeforeQC.pdf", width=15, height=5)
    VlnPlot(scRNAMergedBeforeQC, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, pt.size=0)
    VlnPlot(scRNAMergedBeforeQC, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
    dev.off()
    saveRDS(scRNAMergedBeforeQC, file="0001_scRNAMergedBeforeQC.rds")
  }
  
  # 质控
  scRNA = lapply(scRNA, function(seuratObject){
    seuratObject <- subset(seuratObject, subset=nFeature_RNA>500 & nFeature_RNA<3000 & percent.mt<20)
    return(seuratObject)
  })
  ## 备份质控后的数据
  if(TRUE){
    scRNAMergedAfterQC = merge(scRNA[[1]], scRNA[-1])  # 仅合并但不处理批次效应, 用于对比
    pdf("0001_mergedAfterQC.pdf", width=15, height=5)
    VlnPlot(scRNAMergedAfterQC, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, pt.size=0)
    VlnPlot(scRNAMergedAfterQC, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
    dev.off()
    saveRDS(scRNAMergedAfterQC, file="0001_scRNAMergedAfterQC.rds")
  }
  
  # normalize and identify variable features for each dataset independently
  scRNA = future_lapply(scRNA, function(seuratObject){
    seuratObject = NormalizeData(seuratObject)
    seuratObject = FindVariableFeatures(seuratObject, selection.method="vst", nfeatures=2000)
    return(seuratObject)
  })
  
  # 合并多个seuratObject并处理批次效应
  features = SelectIntegrationFeatures(object.list=scRNA)
  scRNA = future_lapply(scRNA, function(seuratObject){
    seuratObject = ScaleData(seuratObject, features=features)
    seuratObject = RunPCA(seuratObject, features=features)
    return(seuratObject)
  })
  anchors = FindIntegrationAnchors(object.list=scRNA, anchor.features=features, reduction="rpca")
  scRNAMerged = IntegrateData(anchorset=anchors)
  DefaultAssay(scRNAMerged) = "integrated"
  saveRDS(scRNAMerged, file="0001_scRNAMerged.rds")
  
  # 去除批次效应
  if(TRUE){
    # PCA降维预处理
    scRNAMerged = ScaleData(scRNAMerged)
    # PCA降维
    scRNAMerged = RunPCA(scRNAMerged)
    # 确定PCA降维的维度
    p1 = ElbowPlot(scRNAMerged, ndims=25)
    p2 = ElbowPlot(scRNAMerged, ndims=50)
    pdf("0001_elbowPlot.pdf")
    print(p1+p2)
    dev.off()
    
    pc.dim = 1:10
    
    # 对细胞进行聚类，把细胞从一个整体中分出亚群
    scRNAMerged <- RunUMAP(scRNAMerged, reduction="pca", dims=pc.dim)
    scRNAMerged <- FindNeighbors(scRNAMerged, reduction="pca", dims=pc.dim)
    scRNAMerged <- FindClusters(scRNAMerged, resolution=seq(0.1, 2, 0.1))  # 聚类的算法
    p1 <- DimPlot(scRNAMerged, reduction="umap", label=TRUE, repel=TRUE)
    pdf("0001_UMAPOfSample.pdf")
    plot = DimPlot(scRNAMerged, reduction="umap", group.by="orig.ident")
    print(plot)
    dev.off()
    
    DefaultAssay(scRNAMerged) <- "integrated"
    #saveRDS(scRNAMerged, file="0001_scRNAMergedAfterClusters.rds")
  }
  ## 处理仅合并但无QC的数据, 用于对比去批次效应的效果
  if(TRUE){
    scRNAMergedAfterQC = ScaleData(scRNAMergedAfterQC)
    scRNAMergedAfterQC <- FindVariableFeatures(scRNAMergedAfterQC, selection.method="vst", nfeatures=2000)
    scRNAMergedAfterQC = RunPCA(scRNAMergedAfterQC)
    pc.dim = 1:10
    scRNAMergedAfterQC <- RunUMAP(scRNAMergedAfterQC, reduction="pca", dims=pc.dim)
    scRNAMergedAfterQC <- FindNeighbors(scRNAMergedAfterQC, reduction="pca", dims=pc.dim)
    scRNAMergedAfterQC <- FindClusters(scRNAMergedAfterQC, resolution=0.4)
    
    pdf("0001_mergedAfterQCUMAP.pdf")
    plot = DimPlot(scRNAMergedAfterQC, reduction="umap", group.by="orig.ident")
    print(plot)
    dev.off()
  }
  
  # 删除错误的列
  scRNAMerged@meta.data = scRNAMerged@meta.data[, -which(names(scRNAMerged@meta.data)=="diease")]
  # 补充不同样本中的disease信息
  scRNAMerged@meta.data[scRNAMerged@meta.data[["orig.ident"]]=="GSE164449_Ctrl1", "disease"] = "Healthy"
  scRNAMerged@meta.data[scRNAMerged@meta.data[["orig.ident"]]=="GSE164449_RPL1", "disease"] = "RPL"  # Recurrent pregnancy loss
  scRNAMerged@meta.data[scRNAMerged@meta.data[["orig.ident"]] %in% c("CRA002181_Ctrl1", "CRA002181_Ctrl2", "CRA002181_Ctrl3"), "disease"] = "Healthy"
  scRNAMerged@meta.data[scRNAMerged@meta.data[["orig.ident"]] %in% c("CRA002181_RPL1", "CRA002181_RPL2"), "disease"] = "RPL"
  scRNAMerged@meta.data[scRNAMerged@meta.data[["orig.ident"]] %in% c("HRA000237_Ctrl1", "HRA000237_Ctrl2", "HRA000237_Ctrl3"), "disease"] = "Healthy"
  scRNAMerged@meta.data[scRNAMerged@meta.data[["orig.ident"]] %in% c("HRA000237_RPL1", "HRA000237_RPL2", "HRA000237_RPL3"), "disease"] = "RPL"
  scRNAMerged@meta.data[scRNAMerged@meta.data[["orig.ident"]] %in% c("GSE130560_Ctrl1", "GSE130560_Ctrl2", "GSE130560_Ctrl3"), "disease"] = "Healthy"
  scRNAMerged@meta.data[scRNAMerged@meta.data[["orig.ident"]] %in% c("GSE130560_RSA1", "GSE130560_RSA2", "GSE130560_RSA3"), "disease"] = "RSA"
  
  # 最终确定resolution
  # 此处需要使用clustree包分析要用的resolution值
  Idents(scRNAMerged) = "integrated_snn_res.0.4"
  
  
}


# 第二部分, 细胞类型注释
if(TRUE){
  DefaultAssay(scRNAMerged) = "integrated"
  pdf("0001_UMAPOfClusters.pdf")
  DimPlot(scRNAMerged, reduction="umap", group.by="integrated_snn_res.0.4", label=TRUE)
  dev.off()
  # 寻找不同cluster之间的差异基因, 注意, 此时的DefaultAssay(scRNAMerged)=="integrated"
  scRNAMergedMarkers = FindAllMarkers(scRNAMerged, only.pos=TRUE)
  #saveRDS(scRNAMergedMarkers, file="0001_scRNAMergedMarkers0.4.rds")

  # 展示不同cluster之间的差异基因
  top10 = dplyr::top_n(dplyr::group_by(scRNAMergedMarkers, cluster), n=10, wt=avg_log2FC)
  scRNAMergedTemp = subset(scRNAMerged, downsample=300)
  plot = DoHeatmap(scRNAMergedTemp, features=top10$gene)+ NoLegend()
  pdf("0001_heatmapOfTopGene.pdf", width=21, height=21)
  print(plot)
  dev.off()

  # 筛选cluster差异基因集中的免疫细胞标记基因
  for(i in 1:dim(top10)[1]){
    gene = top10$gene[i]
    cluster = top10$cluster[i]
    for(cellType in names(geneList)){
      if(gene %in% geneList[[cellType]]){print(paste(gene, "of", cellType, "in", cluster))}
    }
  }
  
  # 堆叠小提琴图
  geneTemp = c()
  for(g in geneList){
    g = g[g %in% rownames(scRNAMerged)]
    geneTemp = append(geneTemp, g)
  }
  temp = scRNAMerged
  Idents(temp) = "integrated_snn_res.0.4"
  plot = StackedVlnPlot(obj=temp, features=geneTemp)
  
  pdf("0001_geneExpressionOfCelltype.pdf", height=35, width=21)
  print(plot)
  for(cellType in names(geneList)){
    geneSet = geneList[[cellType]]
    plot1 = FeaturePlot(scRNAMerged, features=geneSet, min.cutoff=0, pt.size=0.5)
    print(plot1)
  }
  dev.off()
  
  # 细胞类型注释
  new.cluster.ids <- c("T", "NK", "NK", "T", "NK", "Macrophage", "NK", "NKT", "NKT",
                       "Monocyte", "Unknown", "DC", "B", "Mast", "Mast", "NK",
                       "Macrophage", "Macrophage")
  names(new.cluster.ids) <- levels(scRNAMerged)
  scRNAMerged <- RenameIdents(scRNAMerged, new.cluster.ids)
  plot = DimPlot(scRNAMerged, reduction="umap", label=TRUE, pt.size=0.5)
  pdf("0001_UMAPOfCelltype.pdf")
  print(plot)
  dev.off()
  saveRDS(scRNAMerged, "0001_scRNAMergedAfterAnnotation.rds")
}


# test




# debug





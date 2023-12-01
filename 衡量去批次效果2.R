.libPaths(c("/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/",
            "/home/weiyihu/App/conda/envs/RStudio/lib/R/library"
            #"/home/weiyihu/package/share/4.0"
))
setwd("/data1/weiyihu/endometrium/data")
if(FALSE){
  #library("Seurat", lib.loc="/home/weiyihu/package/share/4.0")
  #library("dplyr", lib.loc="/home/weiyihu/package/share/4.0")
  #library("patchwork", lib.loc="/home/weiyihu/package/share/4.0")
  #library("Matrix", lib.loc="/home/weiyihu/package/share/4.0")
  #library("future")
  #library("future.apply")
  #library("nabor")
  #library("tidyverse", lib="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
  #library("lisi", lib="/home/weiyihu/App/conda/envs/RStudio/lib/R/library")
  #library("SeuratDisk", lib="/home/weiyihu/package/share/4.0")
  #library("sceasy", lib="/home/shumeng/R/x86_64-redhat-linux-gnu-library/4.0/")
}

if(TRUE){
  library("Seurat")
  library("dplyr")
  library("patchwork")
  library("Matrix")
  library("future")
  library("future.apply")
  library("nabor")
  library("tidyverse")
  library("lisi")
  library("SeuratDisk")
  library("sceasy")
}


BatchKL=function(df, dimensionData=NULL, batch="BatchID", replicates=200, n_neighbors=100, n_cells=100){
  # input:
  #   df|\t             data.frame|\t  seuratObject@meta.data
  #   dimensionData|\t  matrix|\t      seuratObject@@reductions$umap@cell.embeddings
  #   batch|\t          str|\t         colname of seuratObject@meta.data
  #   replicates|\t     int|\t         =200, the number of boostrap times
  #   n_neighbors|\t    int|\t         =100, the number of nearest neighbours of cell(from all batchs)
  #   n_cells|\t        int|\t         =100, the number of randomly picked cells
  # output:
  #   float|\t  The lower value denotes the better mixing performance.
  df = df
  dimensionData = dimensionData
  batch = batch
  replicates = replicates
  n_neighbors = n_neighbors
  n_cells = n_cells
  
  set.seed(1)
  if (is.null(dimensionData)){
    umapdata=as.matrix(df[,c("UMAP_1","UMAP_2")])
  }else{
    umapdata=as.matrix(dimensionData)
  }
  batchdata=factor(as.vector(df[, batch]))
  table.batchdata=as.matrix(table(batchdata))[,1]
  tmp00=table.batchdata/sum(table.batchdata)  #proportation of population
  n=dim(df)[1]
  KL=sapply(1:replicates,function(x){
    bootsamples=sample(1:n,n_cells)
    nearest=nabor::knn(umapdata,umapdata[bootsamples,],k=min(5*length(tmp00),n_neighbors))
    KL_x=sapply(1:length(bootsamples),function(y){
      id=nearest$nn.idx[y,]
      tmp=as.matrix(table(batchdata[id]))[,1]
      tmp=tmp/sum(tmp)
      return(sum(tmp*log2(tmp/tmp00),na.rm = T))
    })
    return(mean(KL_x,na.rm = T))
  })
  return(mean(KL))
}


CalLISI=function(emb, meta, col=c('seurat_clusters', 'sample')){
  # input:
  #   emb|\t     matrix|\t      seuratObject@reductions@cell.embeddings
  #   meta|\t    data.frame|\t  seuratObject@meta.data
  #   col|\t     vector|\t      =c("seurat_clusters", "sample"),
  #                             [1] should be clusters and [2] should be value associated with batch
  # output:
  #   result|\t  numeric|\t     higher [2] indicates better performance for batch mixing
  emb = emb
  meta = meta
  col = col
  
  lisi_index <- lisi::compute_lisi(emb, meta, col)
  
  result = list()
  for(i in col){
    result[[i]] = median(lisi_index[[i]])
  }
  result = unlist(result)
  #clisi = median(lisi_index$Cell_TypeID)
  #ilisi_Batch = median(lisi_index$BatchID)
  
  return(result)
}

# 初始化
if(TRUE){
  use_condaenv("scDML")
}


# 读取数据
if(TRUE){
  #scRNAMerged = readRDS(file.path("rds", "0002_tempBdata2.rds"))
  scRNAMerged = readRDS("0003_scRNAMergedHarmony.rds")
  scRNAMerged = readRDS(file.path("rds", "0002_tempBdata2.rds"))
  temp = subset(scRNAMerged, downsample=1000)
}

# CalLIST
if(TRUE){
  seuratObject = temp
  result = CalLISI(emb=seuratObject@reductions$umap@cell.embeddings,
                   meta=seuratObject@meta.data,
                   col=c('seurat_clusters', 'sample'))
  print(result)
}

# BatchKL
if(TRUE){
  seuratObject = temp
  result = BatchKL(df=seuratObject@meta.data,
                   dimensionData=seuratObject@reductions$umap@cell.embeddings,
                   batch="sample",
                   replicates=200,
                   n_neighbors=100,
                   n_cells=100)
  print(result)
}

# convert seuratObject to annDataObject
if(TRUE){
  seuratObject =temp
  
  adata <- convertFormat(seuratObject, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
  print(adata)
  write_h5ad(adata,paste0(h5ad_out_path,"Merge",".h5ad"))
}

# debug
if(TRUE){
  metaData = read.table("temp_scviDeBatchMetaData.tsv", sep='\t')
  metaData[1, 1] = 0
  rownames(metaData) = metaData$V1
  metaData = metaData[, -1]
  colnames(metaData) = metaData[1, ]
  metaData = metaData[-1, ]
  
  reduction = read.table("temp_scviDeBatchReduction.tsv", sep='\t')
  reduction[1, 1] = -1
  rownames(reduction) = reduction[, 1]
  reduction = reduction[, -1]
  colnames(reduction) = reduction[1, ]
  reduction = reduction[-1, ]
  
  
  d1 = CalLISI(emb=reduction,
               meta=metaData,
               col=c('leiden', 'sample'))
  
  d2 = BatchKL(df=metaData,
               dimensionData=reduction,
               batch="sample",
               replicates=200,
               n_neighbors=100,
               n_cells=100)
}

# scvi
# BatchKL  4.61851667213262
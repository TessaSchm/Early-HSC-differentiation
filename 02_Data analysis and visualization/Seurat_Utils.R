#Initialize Functions
PrepDR <- function( # From Seurat
  object,
  genes.use = NULL,
  use.imputed = FALSE,
  assay.type="RNA"
) {
  
  if (length(VariableFeatures(object = object)) == 0 && is.null(x = genes.use)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector
         of genes names in genes.use and retry.")
  }
  if (use.imputed) {
    data.use <- t(x = scale(x = t(x = object@imputed)))
  } else {
    data.use <- GetAssayData(object, assay = assay.type,slot = "scale.data")
  }
  genes.use <- if(is.null(genes.use)) VariableFeatures(object = object) else genes.use # Changed
  genes.use <- unique(x = genes.use[genes.use %in% rownames(x = data.use)])
  genes.var <- apply(X = data.use[genes.use, ], MARGIN = 1, FUN = var)
  genes.use <- genes.use[genes.var > 0]
  genes.use <- genes.use[! is.na(x = genes.use)]
  data.use <- data.use[genes.use, ]
  return(data.use)
  }

PCA_estimate_nPC<-function(data, whereto, k=10, from.nPC = 2, to.nPC=150, by.nPC=5, maxit=200, seed=617) {
  library(missMDA)
  PC <-seq(from = from.nPC, to = to.nPC, by = by.nPC)
  # Init the error matrices
  error1<-matrix(0, nrow = length(c(1:k)), ncol = length(PC))
  error2<-matrix(0, nrow = length(c(1:k)), ncol = length(PC))
  print(paste0(k,"-fold paritioning..."))
  # K-Fold Partitioning
  dgem.kfold<-dismo::kfold(t(data), k=k)
  # SVD-CV based on https://stats.stackexchange.com/questions/93845/how-to-perform-leave-one-out-cross-validation-for-pca-to-determine-the-number-of
  for(i in c(1:k)) {
    print(paste0("k:",i))
    X.train<-t(data[, dgem.kfold!=i])
    X.test<-t(data[, dgem.kfold==i])
    # Find a few approximate singular values and corresponding singular vectors of a matrix.
    print("Running SVD...")
    # Seurat uses IRLBA to do PCA : https://github.com/satijalab/seurat/blob/cec7cb95c73fd6d605723e9af9a1f96eda5635de/R/dimensional_reduction.R
    pca.results<-irlba::irlba(A = X.train, nv = to.nPC, maxit = maxit) # Otherwise, default maxit=100 do not converge
    gl<-pca.results$v
    for(j in 1:length(PC)) {
      print(paste0("Ndims:",PC[j]))
      P<-gl[,1:PC[j]]%*%t(gl[,1:PC[j]])
      # Naive method
      err1<-X.test %*% (diag(dim(P)[1]) - P)
      # Approximate method
      err2<-X.test %*% (diag(dim(P)[1]) - P + diag(diag(P)))
      error1[i,j]<-sum(err1^2)
      error2[i,j]<-sum(err2^2)
      rm(err1)
      rm(err2)
    }
  }
  errors1<-colSums(error1)
  errors2<-colSums(error2)
  nPC=PC[which(errors2 == min(errors2))]
  saveRDS(nPC,whereto)
  # plot(PC,errors1)
  # plot(PC,errors2)
  return(nPC)
}

ConvertGenes <- function(
  object
) {
  flybase.r6.16.gene.conv.table<-read.table("/ddn1/vol1/staging/leuven/stg_00002/lcb/jjans/Fly_Brain/Resources/FlyBase_r6.16_FBgn_2_GeneSymbol.tsv", quote = '', header=F,sep = "\t",dec=".")
  colnames(flybase.r6.16.gene.conv.table)<-c("FBgi.Submitted", "FBgi.Current", "FBgi", "Symbol")
  object$FBgi <- rownames(object)
  object.FBgn2symbol.converted<-merge(x = object, y = flybase.r6.16.gene.conv.table, by = "FBgi")
  row.names(object.FBgn2symbol.converted)<-object.FBgn2symbol.converted$Symbol
  object.FBgn2symbol.converted<-object.FBgn2symbol.converted[,-which(names(object.FBgn2symbol.converted) %in% c("FBgi","FBgi.Current","FBgi.Submitted","Symbol"))]
  object<-object.FBgn2symbol.converted
  object$genes<-NULL
  
  rm(object.FBgn2symbol.converted,flybase.r6.16.gene.conv.table)
  count.data <- object
  return(count.data)
}


Seurat_Scaling <- function(
  object,
  path,
  genes
) {
genes <- c(genes)
markers <- c("Gad1","VAChT","VGlut","mub","ey","prt","pnt","Eaat1","Gs2","Tk","Awh","Vmat","Tdc2","Tbh","roX1","dati","ab","mamo","SoxN","chinmo","pros","Imp","jim")

sce_seurat <- CreateSeuratObject(raw.data = object, min.cells = 3, min.genes=200)
mito.genes<-grep(pattern = "^mt:", x = rownames(x = sce_seurat@data), value = TRUE)
percent.mito<-colSums(sce_seurat@data[mito.genes, ]) / colSums(sce_seurat@data)
sce_seurat<-AddMetaData(object = sce_seurat, metadata = percent.mito, col.name = "percent.mito")

sce_seurat <- NormalizeData(object = sce_seurat,normalization.method = "LogNormalize",scale.factor = 1e4)
sce_seurat <- FindVariableGenes(object = sce_seurat, do.plot = FALSE,x.low.cutoff = 0,x.high.cutoff = 3, y.cutoff = 0.8)

sce_seurat <- ScaleData(object = sce_seurat,vars.to.regress = c("nUMI","percent.mito"),genes.use=c(sce_seurat@var.genes,markers,genes),model.use="negbinom")
saveRDS(sce_seurat, file = paste0(path,"_sce_seurat.RDS"))
return(sce_seurat)
}

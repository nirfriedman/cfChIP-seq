library(cba)
library(ctc)
library(Biobase)
library(NMF)


ClusterMatrixRows = function(M, ClusterDistance = FALSE, 
                             ClusterCorrelation = TRUE) {
  if( ClusterDistance ) {
    M.row.dist = dist(M)
  }
  if( ClusterCorrelation ) {
    M.cor = cor(t(M), use = "pairwise.complete.obs")
    M.cor[is.na(M.cor)] = 0
    M.row.dist = as.dist(1 - M.cor)
  }
  M.row.clust = hclust(M.row.dist, method="complete")

  M.row.clust.opt = M.row.clust
  if(nrow(M) > 2) {
    M.row.order = order.optimal(M.row.dist, M.row.clust$merge)
    M.row.clust.opt$merge = M.row.order$merge
    M.row.clust.opt$order = M.row.order$order
  }
  M.row.clust.opt
}

PClusterMatrixCols = function(M, tempFile = NULL, 
                              PCluster = "~/Bin/PCluster",
                              Args = "") {
  PClusterMatrixRows(t(M), tempFile, PCluster, Args)
}

PClusterMatrixRows = function(M, tempFile = NULL, 
                              PCluster = "~/Bin/PCluster",
                              Args = "") {
  eraseTmp = TRUE
  if(is.null(tempFile)) {
    tempFile = paste0("/tmp/runPCluster-", floor(runif(1)*10000))
    eraseTmp = TRUE
  }
  writeTSV(M, paste0(tempFile,".tsv"))
  cmd = paste0(PCluster," ",Args, " -C R -R 5 -I 2 -O F -b -d ", tempFile, " ", tempFile,".tsv ", tempFile)  
  system(cmd)
  XX = read.table(paste0(tempFile, ".dump"), sep = "\t", as.is = TRUE, header = TRUE)

  score.max = max(XX[,"Score"])
  score.min = min(XX[,"Score"])
  
  Score.to.Height = function(i) {
    2*(-XX[i,"Score"] + score.max)/(-score.min + score.max)
  }
  
  if(eraseTmp) {
    file.remove(paste0(tempFile, c(".tsv", ".dump", ".cdt", ".gtr")))
  }
  
 # catn("After clean up")
  if(0) {
  dend.leaf = as.dendrogram(hclust(as.dist(matrix(c(1,0,1,0), nc =2))))[[1]]
  Dends = list()
 
  for( i in nrow(XX):1 ) {
    n = XX[i,"Size"]
    catn(i)
    if( n == 1 ) {
#      catn(i,XX[i,"ID"], n, XX[i,"Gene"] )
      dend = dend.leaf
      attr(dend,"label") = XX[i,"Gene"]
      attr(dend, "height") = Score.to.Height(i)
      } else {
      i.l = i+1
      i.r = i+2*XX[i.l, "Size"]
#      catn(i,XX[i,"ID"], n, i.l, i.r )
      h = max(Score.to.Height(i),
              max(attr(Dends[[i.l]],"height"),
                  attr(Dends[[i.r]],"height") ))
      dend = merge( Dends[[i.l]], Dends[[i.r]], height = h)
    }
    Dends[[i]] = dend
  }
  catn(summary(Dends[1]))
  
  as.hclust(Dends[[1]])
  } 
  
  hc = hclust(as.dist(matrix(c(1,0,1,0), nc =2)))
  hc$merge = matrix(nr = nrow(M)-1, nc = 2)
  hc$height = rep(0, nrow(M)-1)
  hc$labels = rownames(M)
  k = 1
  NodeToLabel = rep(NA, nrow(XX))
  NodeHeights = rep(0,nrow(XX))
  for( i in nrow(XX):1 ) {
    n = XX[i,"Size"]
#    catn(i,k)
    if( n == 1 ) {
      NodeToLabel[i] =  -which(rownames(M) == XX[i,"Gene"])
      NodeHeights[i] = Score.to.Height(i)
    } else {
      i.l = i+1
      i.r = i+2*XX[i.l, "Size"]
      #      catn(i,XX[i,"ID"], n, i.l, i.r )
      NodeHeights[i] = max(c(Score.to.Height(i),
                             NodeHeights[i.l],
                             NodeHeights[i.r]))
      NodeToLabel[i] = k
      hc$height[k] = NodeHeights[i]
      hc$merge[k,1] = NodeToLabel[i.l]
      hc$merge[k,2] = NodeToLabel[i.r]
      k = k+1
    }
  }
  hc$order = -NodeToLabel[NodeToLabel < 0]
#  return(as.hclust(as.dendrogram(hc)))
  hc
}
  
clusterMatrix = function(M, ClusterDistance = FALSE, 
                         ClusterCorrelation = TRUE) {
  if( ClusterDistance ) {
    M.row.dist = dist(M)
    M.col.dist = dist(t(M))
  }
  if( ClusterCorrelation ) {
    M.row.dist = as.dist(1 - cor(t(M), use = "pairwise.complete.obs"))
    M.col.dist = as.dist(1 - cor(M, use = "pairwise.complete.obs"))
  }
  M.row.clust = hclust(M.row.dist, method="complete")
  M.col.clust = hclust(M.col.dist, method="complete")
  
  M.col.clust.opt = M.col.clust
  if(ncol(M) > 2) {
    M.col.order = order.optimal(M.col.dist, M.col.clust$merge)
    M.col.clust.opt$merge = M.col.order$merge
    M.col.clust.opt$order = M.col.order$order
  }
  M.row.clust.opt = M.row.clust
  if(nrow(M) > 2) {
    M.row.order = order.optimal(M.row.dist, M.row.clust$merge)
    M.row.clust.opt$merge = M.row.order$merge
    M.row.clust.opt$order = M.row.order$order
  }
  list( row = M.row.clust.opt, col = M.col.clust.opt)
}

writeCDT = function(M,hr = NULL, hc = NULL, fileName) {
  D.r = 1+1-1*is.null(hc)
  D.c = 2+1-1*is.null(hr)
  
  M1 = formatC(M)
  M1[is.na(M)] = ""
  M2 = matrix("", nr = nrow(M)+D.r, nc = ncol(M)+D.c)
  M2[(D.r+1):nrow(M2), (D.c+1):ncol(M2)] = M1
  M2[(D.r+1):nrow(M2), D.c-1] = rownames(M)
  M2[(D.r+1):nrow(M2), D.c] = rownames(M)
  M2[1,D.c-1] = "UID"
  M2[1,D.c] = "Name"
  M2[1,(D.c+1):ncol(M2)] = colnames(M)
  if( !is.null(hc) ) {
    M2[2,1] = "AID"
    M2[2,(D.c+1):ncol(M2)] = paste0("ARRY", 0:(ncol(M)-1),"X")
    M2[,(D.c+1):ncol(M2)] = M2[,D.c+hc$order]
  }
  if( !is.null(hr) ) {
    M2[1,1] = "GID"
    M2[(D.r+1):nrow(M2),1] = paste0("GENE", 0:(nrow(M)-1),"X")
    M2[(D.r+1):nrow(M2),] = M2[D.r+hr$order,] 
  }
 
  write.table(M2, file=fileName, quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)
}

writeCluster = function(M,row.clust = NULL, col.clust = NULL, ClusterFilename) {
  catn("Writing cluster to", ClusterFilename)
  if(!is.null(col.clust))
    r2atr(col.clust, file=paste0(ClusterFilename, ".atr"), distance = "euclidean")
  if(!is.null(row.clust))
    r2gtr(row.clust, file=paste0(ClusterFilename, ".gtr"), distance = "euclidean")
  writeCDT(M, row.clust, col.clust, fileName = paste0(ClusterFilename, ".cdt"))
}

writeTSV = function(M1,fileName, EWeight = NULL) {
  catn("Writing TSV to", fileName)
  M = formatC(M1)
  if( is.null(EWeight))
    EWeight = rep(1,ncol(M1))
  M[is.na(M1)] = ""
  M2 = matrix("", nr = nrow(M1)+2, nc = ncol(M1)+2)
  M2[1,1] = "UID"
  M2[1,2] = "Name"
  M2[1,3:ncol(M2)] = colnames(M1)
  M2[2,1] = "EWEIGHT"
  M2[2,2] = "EWEIGHT"
  M2[2,3:ncol(M2)] = EWeight
  M2[3:nrow(M2), 3:ncol(M2)] = M
  M2[3:nrow(M2), 1] = rownames(M1)
  M2[3:nrow(M2), 2] = rownames(M1)
  
  write.table(M2, file=fileName, quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)
}


writeNMFOutput = function(B, nm, NMFFileName, B.norm = NULL, 
                          ClusterRows = TRUE, ClusterCols = TRUE,
                          NormalizeRows = TRUE ) {
  if( is.null(NMFFileName))
    return()

  catn("Saving to", NMFFileName)  
  if( grepl("/$", NMFFileName) )
    if( !dir.exists(NMFFileName))
      dir.create(NMFFileName)
  
  saveRDS(nm, paste0(NMFFileName,"binary.rds"))
  nm.coef = coef(nm)
  nm.basis = basis(nm)
  
  Norm = mean(colSums(nm.coef))
  nm.coef = nm.coef/Norm
  nm.basis = nm.basis*Norm
  if( is.null(B.norm))
    B.norm = rep(1,nrow(nm.basis))
    
  nm.basis = nm.basis * B.norm
  
  C = nm.basis %*% nm.coef

  if(nm@method == "offset")
    C = C + (offset(nm)*B.norm) %*% t(rep(1, ncol(C)))
  
  nm.basis = log2(1+nm.basis)
  nm.basis = nm.basis - rowMedians(nm.basis)
  row.names(nm.coef) = 1:nrow(nm.coef)
  colnames(nm.basis) = row.names(nm.coef)
  
  if(nm@method != "pNMF") {
    cs = consensus(nm)
    l = clusterMatrix(cs)
    p = ggplot(melt(cs[l$row$order, l$row$order]), aes(Var1,Var2, fill = value) )+geom_raster()
    p = p + coord_equal() + scale_fill_gradient(low = "white", high="darkred")
    ggsave(p, filename = paste0(NMFFileName,"consensus.png"))
  }
  
  C = log2(1+C)
  if(NormalizeRows)
    C = C - rowMedians(C)
  writeTSV(C,paste0(NMFFileName,"predict.tsv"))
  writeTSV(nm.basis, paste0(NMFFileName,"basis.tsv"))
  writeTSV(t(nm.coef), paste0(NMFFileName,"coef.tsv"))
  
  if(ClusterRows || ClusterCols) {
    l.basis = clusterMatrix(nm.basis)
    row.order = l.basis$row$order
    if( !ClusterRows) {
      l.basis$row = NULL
      row.order = 1:nrow(nm.basis)
    }
    writeCluster(nm.basis, l.basis$row, l.basis$col,  paste0(NMFFileName,"basis"))
    l.coef = clusterMatrix(t(nm.coef))
    col.order = l.coef$row$order
    if( !ClusterCols) {
      l.coef$row = NULL
      col.order = 1:ncol(nm.coef)
    }
    writeCluster(t(nm.coef), l.coef$row, l.basis$col,  paste0(NMFFileName,"coef"))
    writeCluster(C,l.basis$row,l.coef$row, paste0(NMFFileName,"predict"))
    
    D = C[row.order,col.order]
    E = B[rownames(D),colnames(D)]
    E = log2(E+1)
    E = E - rowMeans(E)
    D = cbind(D, rep(NA,nrow(D)), E )
    writeTSV(D,paste0(NMFFileName,"compare.tsv"))
  }
}

runNMF = function(B, fn, rank, NMFMethod = "offset", B.norm = 1, ClusterRows = TRUE, ClusterCols = TRUE) {
  I = !is.na(colSums(B))
  nm = nmf(B[,I],rank = rank, nrun = 30, method = NMFMethod)
  writeNMFOutput(B, nm, fn, B.norm, ClusterRows = ClusterRows, ClusterCols = ClusterCols)
  return(nm)
}

extendMNF = function(B, rank, fn, nm, NMFMethod= "offset", B.norm = 1, ClusterRows = TRUE, ClusterCols = TRUE) {
  
  nm.coef = coef(nm)
  nm.coef[nm.coef < 1E-6 * max(nm.coef)] = 0 
  nm.basis = basis(nm)
  nm.basis[nm.basis < 1E-8 * max(nm.basis)] = 0
  old.rank = nrow(nm.coef)
  if( rank <= 0 )
    rank = old.rank
  
  if(!all(colnames(nm.coef) %in% colnames(B))  )   
    catn("Mismatch in column names")
  if(!all(rownames(nm.basis) %in% rownames(B))  )   
    catn("Mismatch in row names")
  
  H =  matrix(1,nr = rank, nc = ncol(B))
  colnames(H) = colnames(B)
  H[1:old.rank,colnames(nm.coef)] = nm.coef
  
  W = matrix(1,nr = nrow(B), ncol = rank)
  rownames(W) = rownames(B)
  W[rownames(nm.basis),1:old.rank] = nm.basis
  
  if(nm@method == "offset") {
    O = runif(nrow(B), max = median(offset(nm)))
    names(O) = rownames(B)
    O[rownames(nm.basis)] = offset(nm)
    init = nmfModel(rank, B, W = W, H = H, offset = O, model = "NMFOffset")
    
  } else 
    init = nmfModel(rank, B, W = W, H = H)
  
  nm = nmf(B,rank = rank,method = NMFMethod, seed = init)
  
  writeNMFOutput(B, nm, fn, B.norm, ClusterRows = ClusterRows, ClusterCols = ClusterCols)
  return(nm)
}

runMultiNMF = function(B, OutputFile, NMF.range = 3:10, NMFMethod = "offset",
                       B.norm = 1, ClusterRows = TRUE, ClusterCols = TRUE) {
  nm.list = list()
  for( rank in NMF.range ) {
    catn("NMF rank ", rank)
    NMFFileName = paste0(OutputFile,"-",NMFMethod,"-",rank, "/")
    
    nm = runNMF(B,NMFFileName, rank, NMFMethod, B.norm, ClusterRows, ClusterCols)
    nm.list[[rank]] = nm
  }
  NMF.quality = list()
  NMF.quality[["Dispersion"]] = sapply(NMF.range, function(i) dispersion(nm.list[[i]]))
  NMF.quality[["Cophcor"]] = sapply(NMF.range, function(i) cophcor(nm.list[[i]]))
  NMF.quality[["silhouette.samples"]] = sapply(NMF.range, 
                                               function(i) mean(silhouette(nm.list[[i]])[,3]))
  NMF.quality[["silhouette.features"]] = sapply(NMF.range, 
                                                function(i) mean(silhouette(nm.list[[i]],
                                                                            what = "features")[,3]))
  p = plot_grid(plotlist = lapply(names(NMF.quality),
                                  function(a) qplot(NMF.range,NMF.quality[[a]], geom="line") + labs(title = a)), 
                nrow = 2, ncol = 2 )
  NMFFileName = paste0(OutputFile,"-",NMFMethod,"-")
  ggsave(plot = p, filename = paste0(NMFFileName, "-ranks.png"))
  ranks = sapply(NMF.quality,function(x) NMF.range[which.max(x)])
  catn("Best NMF rank =", ranks)
  
  list(nm = nm.list[[ranks[[1]]]],nm.list = nm.list, ranks = ranks, qual = NMF.quality)
}
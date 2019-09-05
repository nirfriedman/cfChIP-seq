
Max.Matrix.Mult = 20000

ComputeGeneCounts = function (A, W2G = Win2Gene.matrix) {
  if(!is.list(A)) 
    A = list(A)
  breaks = c(seq(0,dim(W2G)[1], Max.Matrix.Mult), dim(W2G)[1])
  L = lapply(1:(length(breaks)-1), function(i) {
    Is = (breaks[i]+1):breaks[i+1]
    M = W2G[Is,] 
    do.call("cbind", lapply(A, function(v) t(v[Is] %*% M)))
  })
  G = Reduce("+", L)
  matrix(G, nr = dim(G)[1], nc = dim(G)[2], dimnames = dimnames(G))
}

if( !exists("Win2Gene.matrix") ) {
  Win2Gene.Matrix.filename = paste0(SetupDIR,"Win2Gene.rds")
  print("Loading Window to Gene mapping")
  L = readRDS(Win2Gene.Matrix.filename)
  Win2Gene.matrix = L$Matrix
  MultiPromoterGenes = L$Multi
  rm(L)
#  colnames(Win2Gene.matrix)[is.na(colnames(Win2Gene.matrix))] = ""
  Genes = colnames(Win2Gene.matrix)
  GeneWindows = as.integer(rownames(Win2Gene.matrix))
  W = width(TSS.windows[GeneWindows])/1000
  GeneLength = ComputeGeneCounts(W)
  names(GeneLength) = Genes
  rm(W)
}

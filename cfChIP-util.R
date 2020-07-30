Genes2Win = function(G) {
  n = length(G)
  m = 1
  nr = nrow(Win2Gene.matrix)
  Ws = c()
  while( m <= n ) {
    m1 = min(n, m+500)
    W = which(as.logical(Win2Gene.matrix[,G[m:m1]] != 0))
    WW = data.frame(W %% nr, floor(W/nr)+1)
    Ws = c(Ws,GeneWindows[WW[,1]])
    m = m + 500
  }
  Ws
}

GenesBrowserWin = function(G) {
#  WW = which( as.logical(Win2Gene.matrix[,G] != 0), arr.ind = TRUE )
#  WW[,1] = as.numeric(rownames(WW))
 
  GW = rep("", length(G))
  n = length(G)
  m = 1
  nr = nrow(Win2Gene.matrix)
  while( m <= n ) {
    m1 = min(n, m+500)
    W = which(as.logical(Win2Gene.matrix[,G[m:m1]] != 0))
    WW = data.frame(W %% nr, floor(W/nr)+1)
    S = split(WW[,1], WW[,2])
    names(S) = G[m:m1]
    Is = GeneWindows[sapply(S, min)]
    Js = GeneWindows[sapply(S, max)]
    St = start(TSS.windows[Is])-10000
    St[St < 0] = 0
    En =  end(TSS.windows[Js])+10000
    GW[m:m1] = paste0(chrom(TSS.windows[Is]),":", St,"-", En)
    m = m + 500
  }
  GW
}


WinsBrowserWin = function(W) {
  W = as.numeric(W)
  paste0(chrom(TSS.windows[W]), ":", pmax(0,start(TSS.windows[W])-10000), "-", 
         end(TSS.windows[W])+10000)
}

GRBrowserWin = function(GR) {
  GR = resize(GR, width=width(GR)+20000, fix = "center")
  Starts = formatC(start(GR))
  Ends = formatC(end(GR))
  Chr = as.character(chrom(GR))
  paste0(Chr,":", Starts, "-", Ends)
}

ArraySubColumns = function(A, Cs, f = "mean") {
  Cs = Cs[Cs %in% colnames(A)]
  if(length(Cs)>1) {
    apply(A[,Cs], 1, f)
  } else {
    if(length(Cs) == 1)
      A[,Cs]
    else
      0
  }
}

GenesIncreasedPVal = function(G,s) {
  Y = GeneCounts[G,s]
  Lam = HealthyNormAvg[G] / QQNorm[s] + GeneBackground[G,s] 
  Pv = ppois(Y-1,Lam, lower.tail = FALSE, log.p = T)
  pvalue = -ppois(sum(Y)-1, sum(Lam),  lower.tail = FALSE, log.p = T)
  list( pvalue = pvalue, PVals = Pv, Obs = Y, Exp = Lam)
}

GeneSigsIncreasedPVal = function(Sigs) {
  sapply(Sigs, function(sig) sapply(colnames(GeneCounts), function(s) GenesIncreasedPVal(sig,s)[[1]]))
}

GeneSigsCounts = function(Sig)  sapply(Sig, 
                                       function(sig) 
                                         if( length(sig) > 1 ) colSums(GeneCounts.QQnorm[sig,])/length(sig) 
                                       else GeneCounts.QQnorm[sig,])

GeneSigsIncreasedCounts = function(Sigs) {
  GeneSigsCounts(Sigs) - 
    sapply(Sigs, function(sig) rep(mean(HealthyNormAvg[sig]), length(Samples)))
}


sigPValue = function(s, sig, W = WINS, B = Background) {
  Bg = sigBackground(s,sig, B)
  Fg = W[sig,s]
  computePValue(Fg, Bg)
}

evalSig = function(sig,  W = WINS, B = Background) {
  Bg = sigBackgroundParallel(sig, B)
  Fg = W[sig,]
  sapply(colnames(W), function(s) computePValue(Fg[,s], Bg[s,]))
}

evalSigPoisson  = function(sig) evalSig(sig)

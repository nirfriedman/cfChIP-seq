
GenesBrowserWin = function(G) {
#  WW = which( as.logical(Win2Gene.matrix[,G] != 0), arr.ind = TRUE )
#  WW[,1] = as.numeric(rownames(WW))
 
  nr = nrow(Win2Gene.matrix)
  W = which(as.logical(Win2Gene.matrix[,G] != 0))
  WW = data.frame(W %% nr, floor(W/nr)+1)
  S = split(WW[,1], WW[,2])
  names(S) = G
  Is = sapply(S, min)
  Js = sapply(S, max)
  St = start(TSS.windows[Is])-10000
  St[St < 0] = 0
  En =  end(TSS.windows[Js])+10000
  paste0(chrom(TSS.windows[Is]),":", St,"-", En)
}


WinsBrowserWin = function(W) {
  sapply(W,
         function(w) {
           paste0(chrom(TSS.windows[w]), ":", max(0,start(TSS.windows[w])-10000), "-", end(TSS.windows[w])+10000)
         })
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
  Pv = ppois(Y,Lam, lower.tail = FALSE, log.p = T)
  pvalue = -ppois(sum(Y), sum(Lam),  lower.tail = FALSE, log.p = T)
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

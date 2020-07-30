#  Definitions


fitNoise = function(X, thresh = 0.95, MinNumber = 50, Prior = NA, PriorStrength = 100) {
  X = X[!is.na(X)]
  if(length(X) < MinNumber) 
    return(NA)
  t = quantile(X,thresh)
  X = X[X <= t]
  T = table(X)
  N = as.integer(names(T))
  M = sum(T)
  
  alpha = 1
  beta = 0
  if( !is.na(Prior) ) {
    beta = PriorStrength
    alpha = Prior*PriorStrength
  }
  
  LL = function(l)  {
    sum( T*dpois(N,l,log = TRUE)) - M*ppois(max(N),l,log.p = TRUE) + (alpha-1)*log(l) - beta*l
  }
  m = mean(X)
  up = max(0.1,m*10)
  down = max(0.0001,m/10)
  op = optimize(LL, c(up,down),maximum = TRUE)
  op$maximum
}

fitNoise.nbin = function(X,thresh = 0.95) {
  t = quantile(X,thresh,na.rm=TRUE)
  X = X[X <= t]
  fit = fitdist(X,"nbinom")
  coef(fit)
}

fitNoise.mean = function(X,thresh = 0.95) {
  t = quantile(X,thresh)
  X = X[X <= t]
  mean(X)
}


if( !exists("Background.Regions")) {
  print("Loading Background model")
  BackgroundModel.filename = paste0(SetupDIR,"BackgroundModel.rds")
  L = readRDS(BackgroundModel.filename)
  Background.Regions.width = L$Background.Regions.width 
  Background.Regions = L$Background.Regions
  Background.Regions.num = L$Background.Regions.num
  Background.Uniq.Regions = L$Background.Uniq.Regions
  Background.Inds.width = L$Background.Inds.width
  Background.Inds = L$Background.Inds 
  Background.Inds.num = L$Background.Inds.num 
  Background.Uniq.Inds = L$Background.Uniq.Inds
  Background.RegionChr = L$Background.RegionChr
  Background.IndRegion = L$Background.IndRegion
  Background.WindowInd = L$Background.WindowInd
  Background.WindowRegion = L$Background.WindowRegion
  WinChr = L$WinChr
  Background.RegionWindows = L$Background.RegionWindows
  Background.IndWindows = L$Background.IndWindows
  rm(L)
}


Background.Ind.Strech = 
Background.windows = TSS.windows$type == "background" & width(TSS.windows) > 4000

Background.windows.width = 5

if (TargetMod == "H3K4me3-scer") { # TODO: set a cutoff for human/yeast in params
  Background.windows = TSS.windows$type == "background" & width(TSS.windows) > 400
  Background.windows.width = 2
}

Edge.Penalty = 0



getBackground = function( mu, w ) {
  mu.genome = mu$genome
  W = TSS.windows[w]
  chr = as.character(chrom(W))
  mu.chr = mu$chr[chr]
  
  mu.region = mu$region.uniq[as.numeric(Background.WindowRegion[w])]
    
  mu.ind = mu$ind.uniq[as.numeric(Background.WindowInd[w])]
  
  len = (width(W)+Edge.Penalty)/1000
  c( ind = mu.ind*len, region = mu.region*len, chr = mu.chr*len, genome = mu.genome*len)
} 



getMultiBackground = function( mu, w ) {
  mu.genome = mu$genome
  W = TSS.windows[w]
  chr = as.character(chrom(W))
  mu.chr = mu$chr[chr]
 
  mu.region = mu$region.uniq[as.numeric(Background.WindowRegion[w])]
  
  mu.ind = mu$ind.uniq[as.numeric(Background.WindowInd[w])]
  
  len = (width(W)+Edge.Penalty)/1000
  rbind( ind = mu.ind*len, region = mu.region*len, chr = mu.chr*len, genome = mu.genome*len)
} 

getMultiBackgroundEstimate = function( mu, w ) {
  W = TSS.windows[w]
  mu.ind = mu$ind.uniq[as.numeric(Background.WindowInd[w])]
  
  len = (width(W)+Edge.Penalty)/1000
  mu.ind*len
} 

getBackgroundParallel = function( MUs, w, terse=TRUE ) {
  W = TSS.windows[w]
  len = (width(W)+Edge.Penalty)/1000
  
  i = as.numeric(Background.WindowInd[w])
  if(!terse) {
    r = as.numeric(Background.WindowRegion[w])
    chr = as.character(chrom(W))
    sapply( MUs, function(mu) {
      mu.genome = mu$genome
      mu.chr = mu$chr[chr]
      mu.region = mu$region.uniq[r]
      
      mu.ind = mu$ind.uniq[i]
      
      c( ind = mu.ind*len, region = mu.region*len, chr = mu.chr*len, genome = mu.genome*len)
    }) 
  } else
    sapply( MUs, function(mu) mu$ind.uniq[i])*len
} 

sigBackground = function(s,sig, B = Background) {
  mu = B[[s]]
  sapply(sig, function(w) getBackground(mu,w)[1])
}

sigBackgroundParallel.Jump = 1000
sigBackgroundParallel = function(sig, B = Background) {
  doWindow = function(w) getBackgroundParallel(B,w, terse = TRUE)
  if( length(sig) <= sigBackgroundParallel.Jump ) { 
    sapply(sig, doWindow)
  } else {
    breaks = c(seq(0,length(sig)-1, by = sigBackgroundParallel.Jump),length(sig))
    do.call(cbind,
            lapply(1:(length(breaks)-1), function(i) {
              print(breaks[i])
              sapply(sig[(breaks[i]+1):breaks[i+1]],doWindow)
              }))
  }
}


#
# Test significance of the observation being Above (or below) the expectation. 
#
# Obs - real counts (integer)
# Bg - background rate
# Exp.mean and Exp.sd - mean and variance of expectation (in normalized units)
# Scale - scale from normalized units to observations (1/Qnorm)
# Above = TRUE if testing for value above expected, FALSE if testing for below
computePValue = function(Obs, Bg, Exp.mean = NULL, Exp.var = NULL, Scale = 1, Above = TRUE) {
  obs = sum(Obs)
  bg = sum(Bg)
  if( is.null(Exp.mean) ) {
    exp.mean = 0
  } else
    exp.mean = sum(Exp.mean)
  
  if( is.null(Exp.var)) {
    exp.var = 0
  } else
    exp.var = sum(Exp.var)
  
  names(Scale) = NULL
  nb.mean = Scale*exp.mean+bg
  nb.var = exp.var*(Scale**2)+nb.mean
  if( exp.var > 0 ) {
    nb.size = (nb.mean**2) /((Scale**2)*exp.var)
    pval = pnbinom(obs - 1*Above, size = nb.size, mu = nb.mean, lower.tail = !Above, log.p = TRUE)
  } else {
    pval = ppois(obs - 1*Above, nb.mean, lower.tail = !Above, log.p = TRUE)    
  }
  c(pv = pval, 
    obs = obs, 
    bg = bg, 
    lambda = nb.mean,  
    exp = exp.mean, 
    exp.var = exp.var, 
    scale = Scale,
    zscore = (obs-nb.mean)/sqrt(nb.var)
  )
}

computePValue0 = function(Obs, Bg) {
  obs = sum(Obs)
  bg = sum(Bg)
  pval = ppois(obs-1, bg, lower.tail = FALSE, log.p = TRUE)    

  c(pv = pval, 
    obs = obs, 
    bg = bg, 
    lambda = bg,  
    exp = 0, 
    exp.var = 0, 
    scale = 1,
    zscore = (obs-bg)/sqrt(bg)
  )
}


computeMultiPValue1 = function(obs, bg, exp.mean = NULL, exp.var = NULL, Scale = 1, Above = TRUE) {
  if( is.null(exp.mean) ) 
    exp.mean = 0
  if( is.null(exp.var)) 
    exp.var = 0
  nb.mean = Scale*exp.mean+bg
  nb.var = exp.var*(Scale**2)+nb.mean
  nb.size = (nb.mean**2) /((Scale**2)*exp.var)
  pvals = ppois(obs-1*Above, nb.mean, lower.tail = !Above, log.p = TRUE)   
  
  I = is.finite(nb.size)
  if( any(I) )
    pvals[I] =  pnbinom(obs[I]-1*Above, size = nb.size[I], mu = nb.mean[I], lower.tail = !Above, log.p = TRUE)
  
  rbind(pv = pvals, 
        obs = obs, 
        bg = bg, 
        lambda = nb.mean,  
        nb.var = nb.var,
        exp = exp.mean, 
        exp.var = exp.var, 
        scale = Scale,
        zscore = (obs-nb.mean)/sqrt(nb.var))
}

computeMultiPValue = function(obs, bg, exp.mean = NULL, exp.var = NULL, Scale = 1, Above = TRUE) {
  if( is.null(exp.mean) ) 
    exp.mean = 0
  if( is.null(exp.var)) 
    exp.var = 0
  nb.mean = Scale*exp.mean+bg
  nb.size = (nb.mean**2) /((Scale**2)*exp.var)
  pvals = ppois(obs-1*Above, nb.mean, lower.tail = !Above, log.p = TRUE)   

  I = is.finite(nb.size)
  if( any(I) )
    pvals[I] =  pnbinom(obs[I]-1*Above, size = nb.size[I], mu = nb.mean[I], lower.tail = !Above, log.p = TRUE)
  
  return(pvals)
}


getPotentionalWin = function(s, W = WINS, Bwindows = Background.windows ) {
  Y = W[,s]
  Y[!Bwindows] = NA
  Y
}

buildBackground = function(s = NULL, W = WINS, Y=NULL, TWin = TSS.windows, Bwindows = Background.windows) {
  if( is.null(Y) )
    Y = getPotentionalWin(s, W, Bwindows)  
  else
    Y[!Bwindows] = NA
  
  mu.genome = fitNoise(Y)
  mu.chr = sapply(ChrList,function(chr) 
    fitNoise(Y[which(as.character(chrom(TWin)) == chr)], 
             Prior=mu.genome, PriorStrength = 1000))
  names(mu.chr) = ChrList
  
  mu.region = sapply(1:Background.Regions.num, function(i) {
    Ws =  Background.RegionWindows[[i]]
    fitNoise(Y[Ws],
             Prior = mu.chr[Background.RegionChr[i]], 
             PriorStrength = 1000)
  })
  I = which(is.na(mu.region))
  mu.region[I] = mu.chr[Background.RegionChr[I]]

  mu.region.uniq = sapply(Background.Uniq.Regions$revmap, function(l) max(mu.region[l]))
  
  mu.ind = sapply(1:Background.Inds.num, function(i) {
    Ws = Background.IndWindows[[i]]
    fitNoise(Y[Ws], Prior = mu.region[Background.IndRegion[i]], PriorStrength = 500)
  })
  I = which(is.na(mu.ind))
  mu.ind[I] = mu.region[Background.IndRegion[I]]
  mu.ind.uniq = sapply(Background.Uniq.Inds$revmap, function(l) max(mu.ind[l]))
  
  mu = list( genome = mu.genome/Background.windows.width, 
             chr = mu.chr/Background.windows.width, 
             region = mu.region/Background.windows.width, 
             ind = mu.ind/Background.windows.width,
             region.uniq = mu.region.uniq/Background.windows.width,
             ind.uniq = mu.ind.uniq/Background.windows.width)
}


ChrEnd = cumsum(seqlengths(genome.seqinfo)[ChrList]/1e+6)
ChrStart = c(0, ChrEnd[1:(length(ChrList)-1)])
names(ChrStart) = ChrList

plotBackground = function( mu, 
                           PlotInd = TRUE, 
                           PlotRegion = TRUE, 
                           PlotGenome = TRUE, 
                           PlotChr = TRUE ) {
  Ws = seq(1,length(TSS.windows),by=round(length(TSS.windows)/500))
  Xs = start(TSS.windows[Ws])/1e+6 + ChrStart[WinChr[Ws]]
  Ls = (width(TSS.windows[Ws])+Edge.Penalty)/1000
#  Z = sapply(Ws, function(w) getBackground(mu, w))
  Z = getMultiBackground(mu, Ws)
  Z = t(t(Z)/Ls)
  p = ggplot(data.frame(x=Xs, 
                        Ind = Z[1,],
                        Region = Z[2,],
                        Chr = Z[3,], 
                        Genome = Z[4,]))+theme_bw()
  if( PlotInd )
    p = p+geom_point(aes(x=x,y=Ind),color="light blue",size=0.2)
  if( PlotRegion )
    p = p+geom_line(aes(x=x,y=Region),color="blue",size=0.2)
  if( PlotChr )
    p = p+geom_line(aes(x=x,y=Chr), color="black" ,size=.5)
  if( PlotGenome )
    p = p+geom_line(aes(x=x,y=Genome), color="dark gray" ,size=1)
#  p = p + ggtitle(s)
  YMax = quantile(Z, 0.99, na.rm = TRUE)*1.1
  p = p+ylim(c(0,YMax))
  p = p+geom_segment(data=data.frame(x = c(0, ChrEnd)),
                     aes(x=x,xend = x, y=0, yend=YMax),
                     color="gray30", size=0.5)
  p = p+scale_x_continuous(breaks=(ChrStart+ChrEnd)/2, labels = gsub("chr", "", ChrList))
  p = p+labs(x = "Chromosome", y = "reads/KB")
  p
}

checkBackground = function(mu, Z) {
  Ls = width(TSS.windows)/1000
  TT = Background.WindowInd
  runValue(TT) = mu$ind.uniq[runValue(TT)]
  X = as.vector(TT)
  Lam = X *Ls
  Y = Z/Ls
  N = 1:100000
  delta = 0.1
  for( lam in seq(0.1,2,by=delta)) {
    N1 =  which(Lam[N] > lam & Lam[N] < lam+delta)
    N = N[Y[N] <= quantile(Y,0.9)]
    l = -lgamma(Z+1) + Z*(log(Lam)) - Lam
    
  }
  
  DensityScatter(X[N],Y[N])
  DensityScatter(Lam[N],Z[N]-Lam[N])
}


library(Biobase)
source(paste0(SourceDIR, "QQNormalize.R")) 
source(paste0(SourceDIR, "MetaPlot.R")) 
source(paste0(SourceDIR,"cfChIP-util.R"))
source(paste0(SourceDIR, "PlotPVal.R"))
source(paste0(SourceDIR, "EstimateGammaPoisson.R"))
source(paste0(SourceDIR, "YieldEstimation-Functions.R"))
source(paste0(SourceDIR, "NMF-util.R"))

mem.maxVSize(vsize = Inf)

#not the right place, but for now
WinDescription = data.frame(Type = TSS.windows$type, 
                            Gene = TSS.windows$name, 
                            Tissue = TSS.windows$tissue, 
                            stringsAsFactors = FALSE)
rownames(WinDescription) = 1:nrow(WinDescription)

Win.notexcluded.name = which(Win.notexcluded)

logSum = function(a, b) {
  m = max(a,b)
  m + log(exp(a-m)+exp(b-m))
}

logAvg = function(a, b) {
  m = max(a,b)
  m + log(exp(a-m)+exp(b-m))-log(2)
}

twoTailedPValue = function(a,b) {
  pmin(a,b) + log(2)
}

cfChIP.Params <- function() {
  list( 
    Save = TRUE,
    DataDir = NULL,
    Background = FALSE,
    GeneCounts = FALSE,
    GeneBackground = FALSE,
    Normalize = FALSE,
    OverExpressedGenes = FALSE,
    reuseSavedData = TRUE,
    Verbose = TRUE,
    TSS.windows = TSS.windows,
    Win2Gene = Win2Gene.matrix,
    GeneWindows = GeneWindows,
    CommonGenes = CommonGenes,
    NormRef = Healthy.GeneCount,
    NormRef.var = Healthy.GeneCount.var,
    WinNormRef = Healthy.WinCount,
    WinNormRef.var = Healthy.WinCount.var,
    Signatures = Win.Sig,
    Signatures.Ref = Win.Sig.Ref,
    MinFragLen = 50,
    MaxFragLen = 800,
    Heatmap.Type = "dot",
    Heatmap.MaxCount = 15,
    Heatmap.Prune.Cols = TRUE,
    Heatmap.Prune.Rows = FALSE,
    Heatmap.Cluster.Cols = TRUE,
    Heatmap.Cluster.Rows = FALSE,
    Heatmap.Qvalue.threshold = 0.0001,
    Heatmap.Filter.threshold = 6,
    PlotPrograms.GlobalHeatmap = TRUE,
    Programs.Eval.Reference = TRUE,
    Heatmap.Detailed.Type = "plain",
    Heatmap.Detailed.Prune.Cols = FALSE,
    Heatmap.Detailed.Prune.Rows = FALSE,
    Scatter.NormalizeByLength = FALSE,
    Scatter.RefSeqOnly = FALSE,
    Scatter.MarkGenePrograms = NULL,
    Scatter.UpperQuantile = 0.995,
    Scatter.LimitUnits = 2,
    Programs = Gene.Programs,
    Programs.Ref = Gene.Programs.Ref,
    Programs.Partition = Gene.Programs.Partition,
    PlotEnrichments.MaxCount = 50,
    PlotSignatures.IndividualHeatmap = FALSE,
    PlotSignatures.MaxCount = 10,
    PlotSignatures.MaxZscore = 20,
    PlotSignatures.Height = 8,
    PlotSignatures.Width = 10,
    PlotPrograms.MaxCount = 10,
    PlotPrograms.IndividualHeatmap = FALSE,
    PlotPrograms.IndividualBarChart = FALSE,
    PlotPrograms.IndividualCSV = FALSE,
    PlotPrograms.IndividualCDT = FALSE,
    PlotPrograms.MaxZscore = 20,
    PlotPrograms.Height = 8,
    PlotPrograms.Width = 10,
    PlotPrograms.RowOrder = NA,
    MetaGene = MetaGene,
    MetaGene.Offset = 5000,
    MetaGene.Width = 25000,
    MetaGene.Tick = 5000,
    MetaGene.Label = "TSS",
    MetaGene.Max = -1,
    MetaGene.Color = "red",
    MetaGene.BGColor = "black",
    MetaEnhancer = MetaEnhancer,
    MetaEnhancer.Offset = 25000,
    MetaEnhancer.Width = 50000,
    MetaEnhancer.Tick = 10000,
    MetaEnhancer.Label = "Enhancer",
    MetaEnhancer.Max = -1,
    MetaEnhancer.Color = "red",
    MetaEnhancer.BGColor = "black",
    #    BackgroundModel = NULL
    QC.PositiveType = "TSS",
    QC.PositiveNucs = 192015, # H3K4me3 in humans
    QC.GenomeLength = 3.3e9,
    QC.SampleVolume = 2, # 2ml per sample
    QC.QuantileCutoff = 0.75, # cutoff of high peak coverage (see below)
    QC.GenomeWeight = 6.6e-12, # molecular weight of 1 genome # CHECK. see http://www.bio.net/bionet/mm/methods/1999-December/080037.html
    QC.NucleosomeLength = 200,
    QC.cfDNA.ng = 10,
    QCFields = paste("Total","Total.uniq","Total.uniq.est","TSS",
                     "%Signal Total","%Background Total",
                     "%Signal TSS","%Background TSS",
                     "lambda.mix.poiss","Mito","%Mito","Global.signal.yield",
                     "Local.signal.yield", "BG.yield","Global.SNR", "Local.SNR","Seq.factor", sep=";")
  )
}



fixCoverage = function( Cov ) {
  for( c in ChrList )
  {
    l = length(Cov[[c]])
    m = seqlengths(genome.seqinfo[c])
    if( l < m )
      Cov[[c]] = append(Cov[[c]], rep(0,m-l))
  }
  Cov
}


cfChIP.BuildFN <- function(Name, param, suff = ".rdata" ) {
  if( is.null(param$DataDir) )
    f = ""
  else {
    f = DataDir
    if( !grepl("/$", DataDir))
      f = paste0(param$DataDir,"/")
  }
  paste0(f, Name, suff)
}

cfChIP.FindFile <- function( filename, param = cfChIP.Params() ) {
  FileType  = NA
  BED.suffixes = c(".bed", ".bed.gz", ".tagAlign", ".tagAlign.gz")
  if( any(sapply(BED.suffixes, function(s) grepl(paste0(s,"$"), filename))))
    FileType = "BED"
  if( grepl(".bw$|bigWig$", filename) )
    FileType = "BW"
  
  if( is.na(FileType) ) {
    for( s in BED.suffixes )
      if( file.exists(paste0(filename, s))) {
        filename = paste0(filename, s)
        FileType = "BED"
      } 
    if( is.na(FileType) && file.exists(paste0(filename, ".bw"))) {
      filename = paste0(filename, ".bw") 
      FileType = "BW"
    }   
  }
  
  if( is.na(FileType ) ) {
    catn(filename, ": Error, cannot determine file type of ",filename)
    return(NULL)
  }
  
  return(list(filename = filename, FileType = FileType))
}

cfChIP.GetRawData = function(filename, param = cfChIP.Params) {
  ll = cfChIP.FindFile(filename, param)
  filename = ll$filename
  FileType = ll$FileType
  
  dat = list()
  if( FileType == "BED") {     
    if( param$Verbose ) catn(filename, ": Reading BED file")
    
    dat$RawBED = import(filename, format = "BED")
    # remove long/short fragments and non-unique copies
    
    # check for single end reads
    if( max(width(dat$RawBED)) <  param$MinFragLen) {
      dat$RawBED = resize(dat$RawBED, width = 166)
    } else 
      dat$RawBED = dat$RawBED[width(dat$RawBED) <= param$MaxFragLen & width(dat$RawBED) > param$MinFragLen]
    
    dat$BED = unique(dat$RawBED)
    dat$Cov = coverage(dat$BED)
  } 
  if( FileType == "BW" ) {     
    if( param$Verbose ) catn(filename, ": Reading BigWig file", filename)
    dat$BW = import(filename)
    dat$Cov = coverage(dat$BW, weight="score")
  } 
  return(dat)
}

cfChIP.GetCoverage  = function(filename, param=cfChIP.Params()) {
  dat = cfChIP.GetRawData(filename, param)
  return(fixCoverage(dat$Cov))
}

cfChIP.ComputeOverExpressed = function(X.norm, Counts, BG, NormRef, NormRef.var, QQNorm) {
  X = X.norm
  Y = Counts
  logX = log2(X+1)
  logH = log2(NormRef+1)
  Diff.log = logX - logH
  
  Lam = NormRef / QQNorm + BG 
  Lam[is.na(Lam)] = 100
  Pv.up = -computeMultiPValue(Y, BG, NormRef, NormRef.var, 1/QQNorm)
  Pv.down = -computeMultiPValue(Y, BG, NormRef, NormRef.var, 1/QQNorm, Above = FALSE)
  Pv = pmax(Pv.up,Pv.down) - log(2)
  Qv = -log10(p.adjust(exp(-Pv),method="fdr"))
  
  data.frame(healthy = logH, 
             sample = logX, 
             obs = Y, 
             exp = Lam, 
             qvalue = Qv, 
             pvalue = Pv/log(10), 
             X = X, 
             H = NormRef, 
             Significant.up = (Qv > 3 & ((logX-logH > 1) & (Y > Lam))),
             Significant.down = (Qv > 3 & ((logX-logH < -1) & (Y < Lam)) ),
             Zscore =  (Y - Lam)/sqrt(NormRef.var/(QQNorm**2) + Lam),
             stringsAsFactors = FALSE)
}

cfChIP.ProcessFile <- function( filename = NULL,
                               dat = NULL,
                               param = cfChIP.Params(),
                               Force = FALSE,
                               HardForce = FALSE,
                               Change = FALSE ) 
{
  if( is.null( filename ) && is.null( dat ))
  {
    catn("Need one of filename or dat be assigned!")
    return(NULL)
  }
  if( is.null( dat ) )
  {
    Name = BaseFileName(filename)
  } else {
    Name = dat$Name
    if( param$Verbose ) catn(Name, ": Processing precomputed data")
  }
  fn = cfChIP.BuildFN(Name, param )
  
  if( is.null(dat) ) {
    if( param$reuseSavedData && file.exists(fn) ) {
      if( param$Verbose ) catn(Name, ": Reading precomputed data", fn)
      dat <- readRDS(fn)
    } else  {
      dat = list(Name = Name, 
                 Cov = NULL,
                 Counts = NULL,
                 Heights = NULL,
                 Background = NULL,
                 GeneCounts = NULL,
                 GeneHeights = NULL,
                 GeneBackground = NULL, 
                 Counts.QQnorm = NULL,
                 GeneCounts.QQnorm = NULL, 
                 QQNorm = 1,
                 OverExpressedGenes = NULL)
      Change = TRUE
    }
  }
  
  
  #enforce dependencies
  param$Normalize = param$Normalize || param$OverExpressedGenes
  param$GeneCounts = param$GeneCounts || param$Normalize
  param$GeneBackground = param$GeneBackground || param$Normalize
  param$Background = param$Background || param$GeneBackground
  

  
  if( HardForce) {
    dat$Counts = NULL
    Force = TRUE
  }
  
  if( Force ) {
    dat$Cov = NULL
    dat$Background = NULL
    dat$GeneCounts = NULL
    dat$GeneBackground = NULL 
    dat$Counts.QQnorm = NULL
    dat$GeneCounts.QQnorm = NULL 
    dat$QQNorm = 1
  }
  
  # Get BED reads into a GenomicRanges object
  if( is.null(dat$Counts) ) {
    dd = cfChIP.GetRawData(filename, param)
    dat$BED = dd$BED
    dat$BW = dd$BW
    dat$Cov  = dd$Cov
    
    if(!is.null(dd$RawBED)) {
      RawBED.dups <- countOverlaps(dat$BED, dd$RawBED, type = "equal")
      dat$DupCount = as.data.frame(table(RawBED.dups))
    }
    
    if( !is.null(dat$BED))
      dat$FragCount = as.data.frame(table(width(dat$BED)), stringsAsFactors = FALSE)
      
    Change = TRUE
 
    if(!is.null(dat$BED)) {
      if( param$Verbose ) catn(Name, ": counting fragment overlap")
      dat$Counts = countOverlaps(query = param$TSS.windows, 
                                 subject = resize(dat$BED, width=100, fix="center"))
    } else
      if( !is.null(dat$BW)) {
        if( param$Verbose ) catn(Name, ": counting BigWig overlap")
        dat$Cov = coverage(dat$BW, weight="score")
        dat$Counts = rep(0,length(param$TSS.windows))
       
        if( !( "chrY" %in% names(dat$Cov)) ) 
          dat$Cov[["chrY"]] = Rle(0,seqlengths(param$TSS.windows)["chrY"])
      
        ChrRle = Rle(chrom(param$TSS.windows))
        ChrStarts = start(ChrRle)
        ChrEnds = end(ChrRle)
        ChrName = as.character(runValue(ChrRle))
        for( i in 1:nrun(ChrRle)) {
          catn(ChrName[i])
          ws = ChrStarts[i]:ChrEnds[i]
          dat$Counts[ws] = aggregate(dat$Cov[[ChrName[i]]], 
                                     ranges(param$TSS.windows)[ws], 
                                     sum)
        }
        # assuming a typical read is 200bp
        dat$Counts = dat$Counts/200
        
        dat$Heights = max(dat$Cov[param$TSS.windows])
        
      } else {
        catn(Name, ": error! cannot compute counts")
        return(dat)
      }
    dat$Heights = rep(0, length(param$TSS.windows))
    dat$Cov = fixCoverage(dat$Cov)
    for( chr in unique(chrom(param$TSS.windows))) {
      ww = which(chrom(TSS.windows) == chr)
      dat$Heights[ww] = max(dat$Cov[TSS.windows[ww]])
    }
  }
  # remove BED, BW, and Cov
  if( !is.null(dat$BED))  
    dat$BED = NULL
  if( !is.null(dat$BW) )
    dat$BW = NULL
  if( !is.null(dat$Cov))
    dat$Cov = NULL
  
  # Background
  if( param$Background && is.null(dat$Background) ) {
    if( param$Verbose ) catn(Name, ": Computing background model")
    dat$Background = buildBackground(Y = dat$Counts, TWin = param$TSS.windows)
    dat$GeneBackground = NULL 
    dat$Counts.QQnorm = NULL
    dat$GeneCounts.QQnorm = NULL
    Change = TRUE
  }
  
  # Gene counts
  if( param$GeneCounts && is.null(dat$GeneCounts)) {
    if( param$Verbose ) catn(Name, ": Computing gene counts")
    dat$GeneCounts =  ComputeGeneCounts(dat$Counts[GeneWindows], param$Win2Gene)
    if(is.matrix(dat$GeneCounts))
      dat$GeneCounts = dat$GeneCounts[,1]
    names(dat$GeneCounts) = Genes
    dat$GeneCounts.QQnorm = NULL
    dat$GeneHeights = MaxGeneCounts(dat$Heights[GeneWindows], param$Win2Gene)
    
    Change = TRUE
  }
  
  # GeneBackground
  if( param$GeneBackground && is.null(dat$GeneBackground)) {
    if( param$Verbose ) catn(Name, ": Computing gene background")
    mu = dat$Background
    Z = getMultiBackgroundEstimate(mu,param$GeneWindows)
    dat$GeneBackground =  ComputeGeneCounts(Z,param$Win2Gene)
    if(is.matrix(dat$GeneBackground))
      dat$GeneBackground = dat$GeneBackground[,1]
    names(dat$GeneBackground) = Genes
    dat$GeneCNV = dat$GeneBackground/(GeneLength*dat$Background$genome)
    
    dat$GeneCounts.QQnorm = NULL
    Change = TRUE
  }

  # Normalize
  if( param$Normalize && is.null(dat$GeneCounts.QQnorm ) ) {
    if( param$Verbose ) catn(Name, ": Normalize")
    
    GeneDiff = pmax(dat$GeneCounts - dat$GeneBackground, 0)
#    GeneDiff[GeneDiff < 2*(dat$GeneBackground)**.5] = 0
    
    A = cbind(GeneDiff,param$NormRef )
    colnames(A) = c("Sample","Ref")
    Qs = QQNormalizeGenes(A, CommonG = param$CommonGenes)
    dat$QQNorm = Qs[1]/Qs[2]
    dat$GeneCounts.QQnorm = GeneDiff * dat$QQNorm
    names(dat$GeneCounts.QQnorm) = Genes
    dat$OverExpressedGenes = NULL
    
    #estimate background throughout
    dat$WinBackground = getMultiBackgroundEstimate(dat$Background, 1:length(param$TSS.windows))
    dat$Counts.QQnorm = pmax(dat$Counts - dat$WinBackground, 0)*dat$QQNorm
    Change=TRUE
  }
  
  #estimate over expressed genes
  if( param$OverExpressedGenes && is.null(dat$OverExpressedGenes)) {
    if( param$Verbose ) catn(Name, ": Computing overexpressed genes")
    
    
    dat$OverExpressedGenes =  cfChIP.ComputeOverExpressed(dat$GeneCounts.QQnorm, 
                                                   dat$GeneCounts, 
                                                   dat$GeneBackground, 
                                                   param$NormRef, 
                                                   param$NormRef.var, 
                                                   dat$QQNorm)
    
    Change=TRUE
  }
  
  if( param$OverExpressedGenes && is.null(dat$OverExpressedWins) && !is.null(param$WinNormRef) ) {
    if( param$Verbose ) catn(Name, ": Computing overexpressed windows")

    dat$OverExpressedWins = cfChIP.ComputeOverExpressed(dat$Counts.QQnorm, 
                                                   dat$Counts, 
                                                   dat$WinBackground, 
                                                   param$WinNormRef, 
                                                   param$WinNormRef.var, 
                                                   dat$QQNorm)
    rownames(dat$OverExpressedWins) = 1:nrow(dat$OverExpressedWins)
    Change=TRUE
  }
  
  # we are done!
  if( Change && param$Save ) {
    if( param$Verbose ) catn(Name, ": Saving data")
    saveRDS(dat,fn)
  }
  
  return( dat )
}

cfChIP.BuildOutputName = function( dat, Dir, Prefix = NULL, Suffix = ".pdf") {
  if( is.null(Dir) ) 
    Dir = "./"
  if( !grepl("/$", Dir))
    Dir = paste0(Dir,"/")
  
  if( is.null(Prefix) )
    Prefix = ""
  
  TargetDir = paste0(Dir, Prefix)
  TargetDir = sub("/$","", TargetDir)
  if( !dir.exists(TargetDir) )
    dir.create(TargetDir)
  
  fname = paste0(TargetDir, "/", dat$Name, Suffix)
  return(fname)
}

cfChIP.BackgroundPlot = function(dat, Dir = NULL, Force = FALSE) {
  fname = cfChIP.BuildOutputName(dat, Dir, "Chr", ".pdf")
  if( Force || !file.exists(fname)) {
    p = plotBackground(dat$Background)
    s = dat$Name
    p = p + ggtitle(s)
    ggsave(plot = p, filename = fname, height = 5, width=20, units = "cm")
  }
}

cfChIP.FragmentLenPlot = function(dat, Dir = NULL, Force = FALSE ) {
  fname = cfChIP.BuildOutputName(dat, Dir, "Fraglen", ".pdf")
  if( Force || !file.exists(fname)) {
    df = dat$FragCount
    colnames(df) = c("length", "freq")
    df[,"length"] = as.numeric(as.character(df[,"length"]))
    df[,"freq"] = 1e6 * df[,"freq"] / sum(df[,"freq"])
    p = ggplot(df, aes(x=length, y=freq)) +
      geom_col()
    p = p + xlab("Length") + ylab("Fragment Per Million")
    name = dat$Name
    p = p + ggtitle(name)
    ggsave(filename = fname, plot = p)
  }
}

cfChIP.EstimatectDNA = function(dat,  Dir = NULL, Force = FALSE ) {
  fname = cfChIP.BuildOutputName(dat, Dir, "ctDNA", ".pdf")
  xx = EstimateCtDNAFraction(dat)
  if( Force || !file.exists(fname)) {
    ggsave(plot = xx$plot, filename = fname)
  }
  return(xx$ctDNA)
}



cfChIP.MetaPlot = function(dat, filename, Dir = NULL,  param = cfChIP.Params(), Force = FALSE) {
  dat$Cov = cfChIP.GetCoverage(filename, param)
  
  fname = cfChIP.BuildOutputName(dat, Dir, "Meta", ".pdf")
  
  if( Force || !file.exists(fname)) {
    mp.list = list()
    if( !is.null(param$MetaGene)) 
      mp.list[["Gene"]] = c(BED = param$MetaGene,
                            Offset = as.numeric(param$MetaGene.Offset),
                            Label = param$MetaGene.Label,
                            Max = as.numeric(param$MetaGene.Max),
                            Color = param$MetaGene.Color,
                            BGColor = param$MetaGene.BGColor,
                            Tick = as.numeric(param$MetaGene.Tick),
                            Width = as.numeric(param$MetaGene.Width))
    
    if( !is.null(param$MetaEnhancer)) 
      mp.list[["Enhancer"]] = list(BED = param$MetaEnhancer,
                                   Label = param$MetaEnhancer.Label,
                                   Max = as.numeric(param$MetaEnhancer.Max),
                                   Color = param$MetaEnhancer.Color,
                                   BGColor = param$MetaEnhancer.BGColor,
                                   Offset = as.numeric(param$MetaEnhancer.Offset),
                                   Tick = as.numeric(param$MetaEnhancer.Tick),
                                   Width = as.numeric(param$MetaEnhancer.Width))
    
    if( length(mp.list) > 0 ) {

      if( !is.null(dat$BED) ) {
        # calculate reads/million ratio
        Norm = 1e6 / length(dat$BED)
      } else
        Norm = dat$QQNorm
      
      p = PlotMeta(Cov = dat$Cov, Norm = Norm, PlotList = mp.list)
      p = p + ggtitle(dat$Name)
      ggsave(plot = p, filename = fname, width=8,height=11, units = "in")
    }
  }
}

cfChIP.WriteTrack = function(dat, filename, param, Dir = NULL, Force = FALSE) {
  dat$Cov = cfChIP.GetCoverage(filename, param)
  fname = cfChIP.BuildOutputName(dat, Dir, Suffix = ".bw")
  if( Force || !file.exists(fname)) {
#    XX.gr = GRanges(dat$Cov*dat$QQNorm, seqinfo = genome.seqinfo)
    XX.gr = GRanges(dat$Cov*dat$QQNorm)
    # seqinfo(XX.gr) = genome.seqinfo
    XX.gr = XX.gr[score(XX.gr) > 0]
    export(XX.gr,fname, format = "BigWig")
  }
}

cfChIP.EvalSig = function(sig, dat) {
  mu = dat$Background
  Bg = getMultiBackgroundEstimate(mu, sig)
  Fg = dat$Counts[sig]
  computePValue(Fg, Bg)
}

cfChIP.EvalMultiSig = function(sig, LL) {
  ZZ = 1:length(TSS.windows) %in% sig
  Bg = sapply(LL, function(dat) sum(dat$WinBackground[ZZ]))
  Fg = sapply(LL, function(dat) sum(dat$Counts[ZZ]))
  computeMultiPValue1(Fg, Bg)
}


cfChIP.EvalProg = function(prog, dat) {
  Bg = dat$GeneBackground[prog]
  Fg = dat$GeneCounts[prog]
  computePValue(Fg, Bg)
}


cfChIP.EvalMultiProg = function(prog, LL) {
  ZZ = Genes %in% prog
  Bg = sapply(LL, function(dat) sum(dat$GeneBackground[ZZ]))
  Fg = sapply(LL, function(dat) sum(dat$GeneCounts[ZZ]))
  computeMultiPValue1(Fg, Bg)
}



cfChIP.EvalSigWithReference = function(sig, dat, Exp.mean, Exp.var = 0) {
  mu = dat$Background
  Bg = getMultiBackgroundEstimate(mu, sig)
  Fg = dat$Counts[sig]
  
  Scale = 1/dat$QQNorm
  pvAbove = computePValue(Fg, Bg, exp.mean = Exp.mean, exp.var = Exp.var, Scale = Scale )
  pvBelow = computePValue(Fg, Bg, exp.mean = Exp.mean, exp.var = Exp.var, Scale = Scale, Above = FALSE )
  xx = c(pvAbove, c(pv.above = as.numeric(pvAbove["pv"]), pv.below = as.numeric(pvBelow["pv"])))
  xx["pv"] = twoTailedPValue(pvAbove["pv"], pvBelow["pv"])
  xx
}


cfChIP.EvalMultiSigWithReference = function(sig, LL, Exp.mean, Exp.var = 0) {
  ZZ = 1:length(TSS.windows) %in% sig
  Bg = sapply(LL, function(dat) sum(dat$WinBackground[ZZ]))
  Fg = sapply(LL, function(dat) sum(dat$Counts[ZZ]))
  
  Scale = sapply(LL, function(dat) 1/dat$QQNorm)
  
  pvAbove =  computeMultiPValue1(Fg, Bg, exp.mean = Exp.mean, exp.var = Exp.var, Scale = Scale)
  pvBelow = computeMultiPValue1(Fg, Bg, exp.mean = Exp.mean, exp.var = Exp.var, Scale = Scale, Above = FALSE)
  
  xx = rbind(pvAbove, pv.above = as.numeric(pvAbove["pv",]), pv.below = as.numeric(pvBelow["pv",]))
  xx["pv",] = twoTailedPValue(pvAbove["pv",], pvBelow["pv",])
  xx
}

cfChIP.EvalProgWithReference = function(prog, dat, Exp.mean, Exp.var = 0) {
  Bg = dat$GeneBackground[prog]
  Fg = dat$GeneCounts[prog]

  Scale = 1/dat$QQNorm
  pvAbove = computePValue(Fg, Bg, Exp.mean = Exp.mean, Exp.var = Exp.var, Scale = Scale )
  pvBelow = computePValue(Fg, Bg, Exp.mean = Exp.mean, Exp.var = Exp.var, Scale = Scale, Above = FALSE )
  xx = c(pvAbove, c(pv.above = as.numeric(pvAbove["pv"]), pv.below = as.numeric(pvBelow["pv"])))
  xx["pv"] = twoTailedPValue(pvAbove["pv"], pvBelow["pv"])
  xx
}


cfChIP.EvalMultiProgWithReference = function(prog, LL, Exp.mean, Exp.var) {
  ZZ = names(LL[[1]]$GeneCounts) %in% prog
  Bg = sapply(LL, function(dat) sum(dat$GeneBackground[ZZ]))
  Fg = sapply(LL, function(dat) sum(dat$GeneCounts[ZZ]))
  Scale = sapply(LL, function(dat) 1/dat$QQNorm)
  

  pvAbove = sapply(1:length(LL), function(i) computePValue(Fg[i], Bg[i], 
                                                           Exp.mean = Exp.mean, 
                                                           Exp.var = Exp.var, 
                                                           Scale = Scale[i] ))
  pvBelow = sapply(1:length(LL), function(i) computePValue(Fg[i], Bg[i], 
                                                           Exp.mean = Exp.mean, 
                                                           Exp.var = Exp.var, 
                                                           Scale = Scale[i], Above = FALSE ))
  xx = rbind(pvAbove, pv.above = as.numeric(pvAbove["pv",]), pv.below = as.numeric(pvBelow["pv",]))
  xx["pv",] = twoTailedPValue(pvAbove["pv",], pvBelow["pv",])
  xx
}


CollateSigEvals = function(Evals, WithReference, Sig.Width, QQNorm) {
  PVals = -sapply(names(Evals),function(s) Evals[[s]][1])
  QVals = -log10(p.adjust(exp(-PVals),method="fdr"))
  Fg = sapply(names(Evals),function(s) Evals[[s]]["obs"])
  Bg = sapply(names(Evals),function(s) Evals[[s]]["bg"])
  Zscores = sapply(names(Evals),function(s) Evals[[s]]["zscore"])
  names(Fg) = names(Evals)
  names(Bg) = names(Evals)
  names(Zscores) = names(Evals)
  
  if( WithReference ) {
    Exp = sapply(names(Evals),function(s) Evals[[s]]["exp"])
    names(Exp) = names(Evals)
    PVals.above = -sapply(names(Evals),function(s) Evals[[s]]["pv.above"])
    PVals.below = -sapply(names(Evals),function(s) Evals[[s]]["pv.below"])
  } else
    Exp = Bg*0
  
  Counts = Fg - Bg
  Counts[Counts<0] = 0
  Counts.norms = Counts * QQNorm
 
  Counts.norms = t(t(Counts.norms) / Sig.Width)
  Exp.norm = t(t(Exp)/Sig.Width)
  
  if( WithReference ) {
    df = data.frame(ObservedCounts = Fg, 
                    Background = formatC(Bg), 
                    NormalizedCounts = formatC(Counts.norms),
                    NormalizedExp = formatC(Exp.norm),
                    pValue = formatC(PVals/log(10)), 
                    pValue.above = formatC(PVals.above/log(10)), 
                    pValue.below = formatC(PVals.below/log(10)), 
                    qValue = formatC(QVals),
                    zscore = formatC(Zscores),
                    stringsAsFactors = FALSE )
  } else
    df = data.frame(ObservedCounts = Fg, 
                    Background = formatC(Bg), 
                    NormalizedCounts = formatC(Counts.norms),
                    pValue = formatC(PVals/log(10)), 
                    qValue = formatC(QVals),
                    zscore = formatC(Zscores),
                    stringsAsFactors = FALSE )
  df
}

cfChIP.EvaluateSignatures = function( dat, Dir = NULL, 
                                      param = cfChIP.Params(), 
                                      Write = TRUE, 
                                      Force = FALSE,
                                      WithReference = FALSE ) {
 
  if( Write ) {
    if( WithReference) {
      fname = cfChIP.BuildOutputName(dat, Dir, "SignaturesVsHealthy", ".csv")
    } else
      fname = cfChIP.BuildOutputName(dat, Dir, "Signatures", ".csv")
  }
  if(is.null(param$Signatures)) {
    message("Error: missing signatures")
    return()
  }
  if( Force || !Write || !file.exists(fname) )  {
    if( param$Verbose ) catn(dat$Name)
    Sig = param$Signatures
    if( WithReference ) {
      Sig.Ref = param$Signatures.Ref
      if(is.null(Sig.Ref)) {
        Sig.Ref = list(
          avg = sapply(Sig, function(w) sum(param$WinNormRef[w])),
          var = 0
        )
        if( !is.null( param$WinNormRef.var) )
          Sig.Ref$var = sapply(Sig, function(w) sum(param$WinNormRef.var[w]))
      }
      Evals = lapply(1:length(Sig),function(i) cfChIP.EvalSigWithReference(Sig[[i]],dat,Sig.Ref$avg[i], Sig.Ref$var[i]))
      names(Evals) = names(Sig)
    } else
      Evals = lapply(Sig,function(s) cfChIP.EvalSig(s,dat))

   Sig.Width = sapply(Sig, function(sig) sum(width(param$TSS.windows[sig])))/1000
   df = CollateSigEvals(Evals, WithReference, Sig.Width, dat$QQNorm)
    
    if( Write ) 
      write.csv(df[order(as.numeric(df$pValue), decreasing = TRUE),], file = fname, quote = FALSE)
    return(df)
  }
}


cfChIP.EvaluatePrograms = function( dat, 
                                    Dir = NULL, 
                                    param = cfChIP.Params(), 
                                    Write = TRUE, 
                                    Force = FALSE,
                                    WithReference = FALSE) {

  if( Write ) 
    if( Write ) {
      if( WithReference) {
        fname = cfChIP.BuildOutputName(dat, Dir, "Programs", ".csv")
      } else
        fname = cfChIP.BuildOutputName(dat, Dir, "ProgramsVsNull", ".csv")
    }

  if( Force || !Write || !file.exists(fname) )  {  
    if(param$Verbose ) catn(dat$Name)
#    O = dat$OverExpressedGenes
    Sig = param$Programs
    if( WithReference ) {
      Sig.Ref = param$Programs.Ref
      if( is.null(Sig.Ref) ) {
        Sig.Ref = list(
          avg = sapply(param$Programs, function(p) sum(param$NormRef[p])),
          var = rep(0, length(param$Programs))
        )
        if( !is.null( param$NormRef.var) )
          Sig.Ref$var = sapply(param$Programs, function(p) sum(param$NormRef.var[p]))
      }
     
      Evals = lapply(1:length(Sig),
                     function(i) cfChIP.EvalProgWithReference(Sig[[i]],dat,Sig.Ref$avg[i], Sig.Ref$var[i]))
      names(Evals) = names(Sig)
    } else
      Evals = lapply(Sig,function(s) cfChIP.EvalProg(s,dat))
    
    Sig.Width = sapply(Sig, length)
    df = CollateSigEvals(Evals, WithReference, Sig.Width, dat$QQNorm)

    if( Write ) 
      write.csv(df[order(as.numeric(df$pValue), decreasing = TRUE),], file = fname, quote = FALSE)
    
    return(df)
  }
}

Max.Qv = 20

if(0) {
cfChIP.OverExpressedGenes = function( dat, Dir = NULL, param = cfChIP.Params(), Force = FALSE ) {
  fname = cfChIP.BuildOutputName(dat, Dir, "DiffGenes", ".csv")
  fname2 = cfChIP.BuildOutputName(dat, Dir, "DiffScatter", ".pdf")
  
  if( Force || !file.exists(fname2)) {
    O = dat$OverExpressedGenes
    O = O[(O$Significant.down|O$Significant.up )& rownames(O) %in% Genes.notexcluded,]
    
    if( nrow(O) > 1 ) {
      if( !is.null(GeneDescription)) {
        O = cbind(O, GeneDescription[rownames(O),])
      } else {
        O$Description = "-"
        O$Tissue = "-"
        O$GTEX.HighExpressed = "-"
        O$GTEX.Expressed = "-"
        O$Suspect = "-"
      }
      df = data.frame(Obs = O$X, 
                      Ref = O$H, 
                      FoldChange = O$sample-O$healthy, 
                      Pvalues = O$pvalue, QValues = O$qvalue,
                      Tissue = O$H3K4me3.Tissue, HighExpressed = O$GTEX.HighExpressed, Expressed = O$GTEX.Expressed,
                      Suspect = O$Suspect, Browser=GenesBrowserWin(rownames(O)),
                      Description = O$Description, stringsAsFactors = FALSE)
      rownames(df) = rownames(O)
      df = df[order(df$FoldChange, decreasing = TRUE),]
      colnames(df) =  c("Observed", "Healthy", "FoldChange (log2)", "p-Value (-log10)",  "Q-Value (-log10)", 
                        "Tissue (Atlas)", "High Expressed (GTEX)", "Expressed (GTEX)",
                        "Suspect?", "Browser window", "Description")
      write.csv(df, file = fname, quote = FALSE)
    }
    
    df = dat$OverExpressedGenes
    df$qvalue = pmin(df$qvalue, Max.Qv )
    
    if( !is.null(param$Scatter.MarkGenePrograms)) {
      MarkPrograms = strsplit(param$Scatter.MarkGenePrograms, ";")[[1]]
      df$program = NA
      for( p in MarkPrograms ) 
        if( p %in% names(param$Programs) ) 
          df[param$Programs[[p]],"program"] = p
    }
    if( as.logical(param$Scatter.RefSeqOnly) ) 
      df = df[which(GeneDescription[rownames(df),"isRefSeq"]),]
    
    
    t.up = as.numeric(param$Scatter.UpperQuantile)
    t.round = as.numeric(param$Scatter.LimitUnits)
    if( as.logical(param$Scatter.NormalizeByLength) ) {
      eps = -4
      df$healthy.norm = log2(2**eps+ (2**df$healthy -1)/GeneLength[rownames(df)])
      df$sample.norm = log2(2**eps+ (2**df$sample -1)/GeneLength[rownames(df)])
#      df$healthy.norm = (2**df$healthy -1)/GeneLength[rownames(df)]
#      df$sample.norm = (2**df$sample -1)/GeneLength[rownames(df)]
      mVal = t.round*ceiling(max(quantile(df$healthy.norm,t.up),quantile(df$sample.norm, t.up))/t.round)
      p = ggplot(df, aes(x=healthy.norm, y=sample.norm, color=qvalue)) 
      p = p+labs( x = paste0("log2(2^", eps," + Healthy/length )"), y = paste0("log2(2^", eps," + ",dat$Name, "/length )"))
      p = p + coord_fixed(ratio=1,xlim=c(eps,mVal), ylim=c(eps,mVal))
      Label.x = eps+2
    } else {
      mVal = t.round*ceiling(max(quantile(df$healthy,t.up),quantile(df$sample, t.up))/t.round)
      p = ggplot(df, aes(x=healthy, y=sample, color=qvalue)) 
      p = p+labs( x = "log2(1 + Healthy )", y = paste("log2(1 +",dat$Name, ")"))
      p = p + coord_fixed(ratio=1,xlim=c(0,mVal), ylim=c(0,mVal))
      Label.x = 2
    }
    p = p + theme_minimal()
    p = p+geom_point(size=0.15, color="gray80")
    p = p+geom_abline(slope = 1, intercept = 0, size=1, color="blue")
    
    Significant = df$Significant.up
    
    if( !is.null(param$Scatter.MarkGenePrograms)) {
      df$significant = Significant
      p = p+geom_point(data=df[!is.na(df$program),], aes(size = significant, color = program)) 
      p = p + scale_size_discrete(range = c(.75,1.5), labels=c("", "significant"))
      p = p + scale_color_brewer(palette = "Set1")
    }  else {
      p = p+geom_point(size=.5, data=df[Significant,])
      p = p+scale_color_gradient(high="red", low="gray30", limits=c(3,Max.Qv), na.value = "black", name="Q-value\n(-log10)" )
    }
    
    p = p+geom_density2d(colour="black", show.legend = FALSE, bins=20)
    p = p + theme(axis.line = element_line(color="black"))
    p = p + theme(axis.ticks = element_line(color="black"))
    p = p + theme(panel.grid =  element_blank())
    p = p + geom_text(data=data.frame(x = Label.x,y=mVal, 
                                      label = paste(sum(rownames(df)[Significant] %in% Genes.notexcluded),"significant genes")),
                      aes(x,y,label=label), color="black")
    ggsave(p, filename = fname2, width=6, height=6)
  }
}
}

cfChIP.WriteOverExpressed = function(O, fname, Description = NULL, NotExcluded = NULL, 
                                     BrowserFun = NULL) {
  O = O[(O$Significant.down|O$Significant.up ),]
  if( nrow(O) > 1 && !is.null(NotExcluded) )
    O = O[rownames(O) %in% NotExcluded,]
  
  if( nrow(O) > 1 ) {
    if( !is.null(Description)) {
      if( ncol(Description) == 1 ) {
        O = cbind(O, Description = Description[rownames(O),])
      } else
        O = cbind(O, Description[rownames(O),])
    } 
    # fill in missing columns
    for(n in c("Description", "Tissue", "GTEX.HighExpressed", 
               "GTEX.Expressed", "Suspect", "Type", "Gene") )
      if( !(n %in% colnames(O)))
        O[,n] = ""

    if( !is.null(BrowserFun) ) {
       BrowserWin = do.call(BrowserFun,list(rownames(O)))
    } else
      BrowserWin = ""
    
    df = data.frame(Obs = O$X, 
                    Ref = O$H, 
                    FoldChange = O$sample-O$healthy, 
                    Pvalues = O$pvalue, QValues = O$qvalue,
                    Tissue = O$Tissue, 
                    HighExpressed = O$GTEX.HighExpressed, 
                    Expressed = O$GTEX.Expressed,
                    Suspect = O$Suspect,
                    Type = O$Type,
                    Gene = O$Gene,
                    Browser=BrowserWin,
                    Description = O$Description, 
                    stringsAsFactors = FALSE)
    rownames(df) = rownames(O)
    df = df[order(df$FoldChange, decreasing = TRUE),]
    colnames(df) =  c("Observed", "Healthy", "FoldChange (log2)", "p-Value (-log10)",  "Q-Value (-log10)", 
                      "Tissue (Atlas)", "High Expressed (GTEX)", "Expressed (GTEX)",
                      "Suspect?","Type","Gene(s)", "Browser window", "Description")
    SelectCol = sapply(colnames(df), function(n) any(df[,n] != "", na.rm = TRUE))
    df = df[,SelectCol]
    for(n in colnames(df))
      if( is.numeric(df[,n]))
        df[,n] = formatC(df[,n])
    write.csv(df, file = fname, quote = FALSE)
  }
}

cfChIP.OverExpressedScatter = function(df, fname2, 
                                       Name,
                                       MarkPrograms = NULL,
                                       Subset = NULL,
                                       t.up = 0.995,
                                       t.round = 2,
                                       NormalizeByLength = NULL,
                                       notExcluded = NULL, 
                                       xlabel = "Healthy") {
  
  df = df[df$healthy > 0 | df$sample > 0,]
  df$qvalue = pmin(df$qvalue, Max.Qv )
  
  if( !is.null(MarkPrograms)) {
    df$program = NA
    for( p in MarkPrograms ) 
      df[MarkPrograms[[p]],"program"] = p
  }
  
  if( !is.null(Subset) )
    df = df[rownames(df) %in% Subset,]
  
  
  mVal = t.round*ceiling(max(quantile(df$H,t.up),quantile(df$X, t.up))/t.round)
  p = ggplot(df, aes(x=H, y=X, color=qvalue)) 
  p = p+labs( x = xlabel, y = Name)
  p = p + coord_fixed(ratio=1,xlim=c(0,mVal), ylim=c(0,mVal))
  if( !is.null(NormalizeByLength) ) {
    eps = 2**-2
    eps.trans = scales::trans_new(name = "logep", 
                                  transform = function(x) log(x+eps),
                                  inverse = function(z) exp(z)-eps)
                                    
#    df$healthy.norm = log2(2**eps+ df$H/NormalizeByLength[rownames(df)])
#    df$sample.norm = log2(2**eps+ df$X/NormalizeByLength[rownames(df)])
#    mVal = t.round*ceiling(max(quantile(df$healthy.norm,t.up),quantile(df$sample.norm, t.up))/t.round)
#    p = ggplot(df, aes(x=healthy.norm, y=sample.norm, color=qvalue)) 
#    p = p+labs( x = paste0("log2(2^", eps," + ", xlabel, "/length )"), y = paste0("log2(2^", eps," + ",Name, "/length )"))
    p = p + scale_x_continuous(trans = eps.trans, breaks=c(0,1,3,10,30,100,300,1000,3000,10000,30000,100000))
    p = p + scale_y_continuous(trans = eps.trans, breaks=c(0,1,3,10,30,100,300,1000,3000,10000,30000,100000))
    Label.x = 2
  } else {
    p = p + scale_x_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300,1000,3000,10000,30000,100000))
    p = p + scale_y_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300,1000,3000,10000,30000,100000))
    Label.x = 1.5
  }
  p = p + theme_minimal()
  p = p+geom_point(size=0.15, color="gray80")
  p = p+geom_abline(slope = 1, intercept = 0, size=1, color="blue")
  
  Significant = df$Significant.up
  
  if( !is.null(MarkPrograms)) {
    df$significant = Significant
    p = p+geom_point(data=df[!is.na(df$program),], aes(size = significant, color = program)) 
    p = p + scale_size_discrete(range = c(.75,1.5), labels=c("", "significant"))
    p = p + scale_color_brewer(palette = "Set1")
  }  else {
    p = p+geom_point(size=.5, data=df[Significant,])
    p = p+scale_color_gradient(high="red", low="gray30", limits=c(3,Max.Qv), na.value = "black", name="Q-value\n(-log10)" )
  }
  
  p = p+geom_density2d(colour="black", show.legend = FALSE, bins=20)
  p = p + theme(axis.line = element_line(color="black"))
  p = p + theme(axis.ticks = element_line(color="black"))
  p = p + theme(panel.grid =  element_blank())
  
  if( !is.null(notExcluded)) {
    nSig = sum(rownames(df)[Significant] %in% notExcluded)
  } else
    nSig = sum(df[,Significant])
  
  p = p + geom_text(data=data.frame(x = Label.x,y=mVal*1.05, 
                                    label = paste(nSig,"significant")),
                    aes(x,y,label=label), vjust = 1, color="black")
  ggsave(p, filename = fname2, width=6, height=6)
}




cfChIP.OverExpressedGenes = function( dat, Dir = NULL, param = cfChIP.Params(), Force = FALSE ) {
  fname = cfChIP.BuildOutputName(dat, Dir, "DiffGenes", ".csv")
  fname2 = cfChIP.BuildOutputName(dat, Dir, "DiffScatter", ".pdf")
  
  if( is.null(dat$OverExpressedGenes ))
    return()
  
  if( Force || !file.exists(fname2)) {
    
    cfChIP.WriteOverExpressed(dat$OverExpressedGenes, 
                              fname, 
                              Description = GeneDescription, 
                              NotExcluded = Genes.notexcluded, 
                              BrowserFun = GenesBrowserWin )
      
    MarkPrograms = NULL 
    if( !is.null(param$Scatter.MarkGenePrograms) )
      MarkPrograms = param$Programs[strsplit(param$Scatter.MarkGenePrograms, ";")]


    GeneSubset = NULL
    if( as.logical(param$Scatter.RefSeqOnly) ) 
      GeneSubset = rownames(GeneDescription)[GeneDescription[,"isRefSeq"]]
    
    t.up = as.numeric(param$Scatter.UpperQuantile)
    t.round = as.numeric(param$Scatter.LimitUnits)
    
    NormGeneLength = NULL
    if( as.logical(param$Scatter.NormalizeByLength) )
      NormGeneLength = GeneLength
    
    cfChIP.OverExpressedScatter(dat$OverExpressedGenes, 
                                fname2,
                                dat$Name,
                                MarkPrograms,
                                GeneSubset,
                                t.up,
                                t.round,
                                NormGeneLength,
                                Genes.notexcluded
    )
        
  }
}


cfChIP.OverExpressedWins = function( dat, Dir = NULL, param = cfChIP.Params(), Force = FALSE ) {
  fname = cfChIP.BuildOutputName(dat, Dir, "DiffWins", ".csv")
  fname2 = cfChIP.BuildOutputName(dat, Dir, "DiffWinScatter", ".png")
  
  if( is.null(dat$OverExpressedWins) )
    return()
  
  if( Force || !file.exists(fname2)) {
    
    cfChIP.WriteOverExpressed(dat$OverExpressedWins, 
                              fname,  
                              NotExcluded = Win.notexcluded.name, 
                              BrowserFun = WinsBrowserWin )
      
    t.up = as.numeric(param$Scatter.UpperQuantile)
    t.round = as.numeric(param$Scatter.LimitUnits)
    
    NormWinLength = NULL
    if( as.logical(param$Scatter.NormalizeByLength) )
      NormWinLength = width(TSS.windows)/1000
    
    cfChIP.OverExpressedScatter(dat$OverExpressedWins, 
                                fname2,
                                dat$Name,
                                NULL,
                                NULL,
                                t.up,
                                t.round,
                                NormWinLength,
                                Win.notexcluded.name
    )
  }
}

cfChIP.EvaluateProgramsHypG = function( dat, Dir = NULL, param = cfChIP.Params(), Write = TRUE, Force = FALSE) {
  if( Write ) 
    fname = cfChIP.BuildOutputName(dat, Dir, "Programs-HypG", ".csv")
  
  if( Force || !Write || !file.exists(fname) )  {  
    O = dat$OverExpressedGenes
    UpGenes = rownames(O)[O$Significant.up & rownames(O) %in% Genes.notexcluded]
    Sig = param$Programs
    
    pValues = list()
    Overlaps = list()
    Obs  = list()
    for( t in names(Sig)) {
      sig = Sig[[t]]
      x = sum(UpGenes %in% sig)
      m = length(UpGenes)
      n = length(Genes) - m
      k = length(sig)
      if( x > 0 && m > 0 ) {
        p = -phyper(x,m,n,k,lower.tail = FALSE, log.p = TRUE)
      } else
        p = 0
      pValues[[t]] = p
      Obs[[t]] = x
      Overlaps[[t]] = 100*x/k
    }
    
    pValues = unlist(pValues)
    Obs = unlist(Obs)
    Overlaps = unlist(Overlaps)
    qValues = -log10(p.adjust(exp(-pValues),method="fdr"))
    
    df = data.frame(ObservedCounts = Obs, N = sapply(Sig, length), K = length(UpGenes), 
                    Fraction = Overlaps, pValue =pValues/log(10), qValue = qValues )
    if( Write ) 
      write.csv(df[order(Overlaps, decreasing = TRUE),], file = fname, quote = FALSE)

    return(df)
  }

}

cfChIP.WriteEnrichR = function( dat, Dir = OutputDIR, Force = FALSE, Write = TRUE ) {
  if( Write ) 
    fname = cfChIP.BuildOutputName(dat, Dir, "EnrichR", ".csv")
  if( Force || !file.exists(fname)) {
    catn("Evaluating EnrichR annotations to", fname)
    O = dat$OverExpressedGenes
    O = O[order(O$pvalue, decreasing = TRUE),]
    UpGenes = rownames(O)[O$Significant.up & rownames(O) %in% Genes.notexcluded]
    UpGenes = UpGenes[GeneDescription[UpGenes, "isRefSeq"]]
    if(length(UpGenes > 10)) {
      # this is funny looking, but we do not want initiate enrichR interface if not used
      source(paste0(SourceDIR,  "EnrichR-Interface.R"))
      
      df = CollectEnrichRAnnotation( UpGenes )
      if(Write && nrow(df) > 0) {
        odf = df
        odf$P.value = formatC(df$P.value)
        odf$Adjusted.P.value = formatC(df$Adjusted.P.value)
        odf$Combined.Score = formatC(df$Combined.Score)
        write.csv(odf, file = fname, row.names = FALSE, quote = FALSE)
      }
    }
  }
  df
}

cfChIP.plotSignatures = function(LL, outputPlotSignatures, param= cfChIP.Params(), 
                                 DetailedPlots = FALSE ) {
  if( param$Verbose ) catn("Evaluating signatures to ", outputPlotSignatures)
  
  EE = lapply(LL, function(dat) cfChIP.EvaluateSignatures(dat, param = param, Write = FALSE) )
  
  PVals = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$pValue))))
  Zscores = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$zscore))))
  Counts.norms = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$NormalizedCounts))))
  rownames(PVals) = names(LL)
  colnames(PVals) = rownames(EE[[1]])
  dimnames(Counts.norms) = dimnames(PVals) 
  dimnames(Zscores) = dimnames(PVals)
  
  PVals[Counts.norms < 0.25] = 0
  
  l = plotPVal(PVals, Counts.norms, 
               max.Count = as.numeric(param$PlotSignatures.MaxCount),
               MapType =  param$Heatmap.Type,
               Filter.Rows = as.logical(param$Heatmap.Prune.Rows),
               Filter.Cols = as.logical(param$Heatmap.Prune.Cols),
               colOrder = "cluster" 
  )
  
  ggsave(plot=l$plot.Q +coord_equal(), 
         filename = paste0(outputPlotSignatures, ".pdf"),
         width=as.numeric(param$PlotSignatures.Width),
         height=as.numeric(param$PlotSignatures.Height))
  
#  ggsave(plot=l$plot.Q + coord_equal(), filename = paste0(outputPlotSignatures, ".pdf"),width=7,height=10)
  write.csv(PVals, file = paste0(outputPlotSignatures, "-pValues.csv") )
  write.csv(Counts.norms, file = paste0(outputPlotSignatures, "-Counts.csv") )
  write.csv(Zscores, file = paste0(outputPlotSignatures, "-Zscores.csv") )
  
  l = plotPVal(PVals, Zscores, 
               max.Count = as.numeric(param$PlotSignatures.MaxZscore),
               MapType =  param$Heatmap.Type,
               Filter.Rows = as.logical(param$Heatmap.Prune.Rows),
               Filter.Cols = as.logical(param$Heatmap.Prune.Cols),
               colOrder = "cluster" 
  )
  
  ggsave(plot=l$plot.Q +coord_equal(), 
         filename = paste0(outputPlotSignatures, "-zscore.pdf"),
         width=as.numeric(param$PlotSignatures.Width),
         height=as.numeric(param$PlotSignatures.Height))
  Sig = param$Signatures
  if( as.logical(param$PlotSignatures.IndividualHeatmap)) 
    for( sigName in names(Sig) ) {
      if( param$Verbose ) catn("Ploting ", sigName, "signature")
      
      W = Sig[[sigName]]
      Bg = sapply(LL, function(dat) getMultiBackground(dat$Background, W)[1,])
      Fg = sapply(LL, function(dat) dat$Counts[W])
      Cg = sapply(LL, function(dat) dat$Counts.QQnorm[W])
      rownames(Bg) = W
      rownames(Fg) = W
      rownames(Cg) = W
      PV = matrix(ppois(Fg-1,Bg,lower.tail = FALSE, log.p = TRUE),
                  nc = ncol(Fg), nr = nrow(Fg), dimnames = dimnames(Fg))
      rowOrder = order(rowSums(Cg))
      rownames(Cg) = paste(rownames(Cg), TSS.windows[W]$type, TSS.windows[W]$name)
      rownames(PV) = rownames(Cg)  
      Cg = rbind(Cg[rowOrder,], signature = Counts.norms[,sigName])
      PV = rbind(PV[rowOrder,], signature = PVals[,sigName])
      ll = plotPVal(-PV/log(10),Cg, 
                    MapType = param$Heatmap.Detailed.Type,
                    Filter.Rows = as.logical(param$Heatmap.Detailed.Prune.Rows),
                    Filter.Cols = as.logical(param$Heatmap.Detailed.Prune.Cols)
      )
      ggsave(plot=ll$plot.Q+guides(fill=FALSE)+coord_equal(), 
             filename = paste0(outputPlotSignatures,"-",sigName,".pdf"),width=7,height=10)
    }
}


cfChIP.plotPrograms = function(LL, 
                               outputPlotPrograms, 
                               param= cfChIP.Params() ) {
  if( param$Verbose ) catn("Evaluating programs to ", outputPlotPrograms)
  
  EE = lapply(LL, function(dat) cfChIP.EvaluatePrograms(dat, param = param, 
                                                        Write = FALSE, 
                                                        WithReference = as.logical(param$Programs.Eval.Reference)) )
  
  df = do.call(rbind, lapply(names(LL), function(n) {
    df = data.frame(Name = LL[[n]]$Name, 
                    Sig = rownames(EE[[n]]), 
                    stringsAsFactors = FALSE)
    cbind(df,EE[[n]]) }
    )
  )
  
#  if( any(df$pValue >= 3) ) {
#    df = df[df$pValue >= 3,]
    write.csv(df, file = paste0(outputPlotPrograms, "-details.csv"), 
              quote=FALSE, row.names = FALSE)
#  }
  
  PVals = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$pValue))))
  Zscores = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$zscore))))
  
  Counts.obs = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$NormalizedCounts))))
  if( as.logical(param$Programs.Eval.Reference) ) {
    Counts.expected = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$NormalizedExp))))
  } else
    Counts.expected = 0*Counts.obs
  
  Counts.norms = Counts.obs - Counts.expected
  rownames(PVals) = names(EE)
  colnames(PVals) = rownames(EE[[1]])
  dimnames(Counts.norms) = dimnames(PVals)
  dimnames(Counts.obs) = dimnames(PVals)
  dimnames(Counts.expected) = dimnames(PVals)
  dimnames(Zscores) = dimnames(PVals)
  
  QVals = PVals
  QVals[] = -log10(p.adjust(exp(-log(10)*PVals), method="fdr"))
  
  if( as.logical(param$PlotPrograms.GlobalHeatmap)) {
    if( is.null(param$Programs.Partition) || length(param$Programs.Partition) == 1) {
      l = plotPVal(PVals, Counts.norms, 
                   max.Count = as.numeric(param$PlotPrograms.MaxCount),
                   MapType = param$Heatmap.Type,
                   Filter.Rows = as.logical(param$Heatmap.Prune.Rows),
                   Filter.Cols = as.logical(param$Heatmap.Prune.Cols),
                   rowOrder = param$PlotPrograms.RowOrder,
                   Filter.threshold = as.numeric(param$Heatmap.Filter.threshold),
                   Q.val.threshold = as.numeric(param$Heatmap.Qvalue.threshold),
                   colOrder = "cluster"
      ) 
    } else {
      l.list = lapply(param$Programs.Partition, function(x)
        plotPVal(PVals[,x], Counts.norms[,x],  max.Count = as.numeric(param$PlotPrograms.MaxCount),
                 MapType = param$Heatmap.Type,
                 Filter.Rows = as.logical(param$Heatmap.Prune.Rows),
                 Filter.Cols = as.logical(param$Heatmap.Prune.Cols),
                 Filter.threshold = as.numeric(param$Heatmap.Filter.threshold),
                 Q.val.threshold = as.numeric(param$Heatmap.Qvalue.threshold),
                 colOrder = "cluster" )
        )
      o.c = unlist(sapply(l.list, function(ll) colnames(ll$PVals)))
      l = plotPVal(PVals[,o.c], Counts.norms[,o.c],  
                   max.Count = as.numeric(param$PlotPrograms.MaxCount),
                   MapType = param$Heatmap.Type,
                   Filter.Rows = as.logical(param$Heatmap.Prune.Rows),
                   Filter.Cols = as.logical(param$Heatmap.Prune.Cols),
                   rowOrder = param$PlotPrograms.RowOrder,
                   Filter.threshold = as.numeric(param$Heatmap.Filter.threshold),
                   Q.val.threshold = as.numeric(param$Heatmap.Qvalue.threshold))
    }
    
    ggsave(plot=l$plot.Q + coord_equal(), 
           filename = paste0(outputPlotPrograms, ".pdf"),
           width=as.numeric(param$PlotPrograms.Width),
           height=as.numeric(param$PlotPrograms.Height))
    
    
    write.csv(QVals, file = paste0(outputPlotPrograms, "-qValues.csv") )
    write.csv(Counts.norms, file = paste0(outputPlotPrograms, "-Counts.csv") )
    write.csv(Zscores, file = paste0(outputPlotPrograms, "-Zscores.csv") )
    
    l = plotPVal(PVals, Zscores, 
                 max.Count = as.numeric(param$PlotPrograms.MaxZscore),
                 MapType = param$Heatmap.Type,
                 Filter.Rows = as.logical(param$Heatmap.Prune.Rows),
                 Filter.Cols = as.logical(param$Heatmap.Prune.Cols),
                 colOrder = "cluster"
    )
    
    ggsave(plot=l$plot.Q + coord_equal(), 
           filename = paste0(outputPlotPrograms, "-zscore.pdf"),
           width=as.numeric(param$PlotPrograms.Width),
           height=as.numeric(param$PlotPrograms.Height))
  }
  # catn("foobar")
  if( as.logical(param$PlotPrograms.IndividualHeatmap) ||
      as.logical(param$PlotPrograms.IndividualCSV) ||
      as.logical(param$PlotPrograms.IndividualCDT) ||
      as.logical(param$PlotPrograms.IndividualBarChart) ) {
    HH  = param$Programs
    for( sigName in names(HH)) {
      catn(sigName)
      if(  as.logical(param$PlotPrograms.IndividualHeatmap) ||  as.logical(param$PlotPrograms.IndividualCSV) ) {
        G = HH[[sigName]]
        PV = sapply(LL, function(dat) dat$OverExpressedGenes[G, "pvalue"])
        rownames(PV) = G
        C = sapply(LL, function(dat) dat$GeneCounts.QQnorm[G] - param$NormRef[G])
        rowOrder = rev(order(rowSums(C)))
        tryCatch({        
          rowOrder = hclust(dist(C))$order
          })
        PV = rbind(program = t(PVals[colnames(PV),sigName]),
                   t(PVals[colnames(PV),sigName]),
                   t(PVals[colnames(PV),sigName]),
                   rep(NA,ncol(PV)),
                   rep(NA,ncol(PV)),
                   rep(NA,ncol(PV)),
                   PV[rev(rowOrder),])
        
        C = rbind(program = t(Counts.norms[colnames(PV),sigName]),
                  t(Counts.norms[colnames(PV),sigName]),
                  t(Counts.norms[colnames(PV),sigName]),
                  rep(NA,ncol(PV)),
                  rep(NA,ncol(PV)),
                  rep(NA,ncol(PV)),
                  C[rev(rowOrder),colnames(PV)]
        )
        if( param$PlotPrograms.IndividualCSV ) {
          write.csv(PV, file = paste0(outputPlotPrograms,"-",sigName,"-pvalues.CSV"), quote = FALSE)
          write.csv(C, file = paste0(outputPlotPrograms,"-",sigName,"-counts.CSV"), quote = FALSE)
        }
        
        catn("CDT - ",as.logical(param$PlotPrograms.IndividualCDT),param$PlotPrograms.IndividualCDT)
        if( as.logical(param$PlotPrograms.IndividualCDT) ) {
          hc.r = ClusterMatrixRows(C)
          hc.c = ClusterMatrixRows(t(C))
          writeCluster(C, row.clust = hc.r,col.clust = hc.c, ClusterFilename = paste0(outputPlotPrograms,"-",sigName))
        }
        
        if(  as.logical(param$PlotPrograms.IndividualHeatmap) ) {
          if(1) {
            l3 = plotPVal(-PV/log(10), PV*sign(C), 
                          max.Count = 10,
                          MapType = param$Heatmap.Detailed.Type,
                          Filter.Rows = as.logical(param$Heatmap.Detailed.Prune.Rows),
                          Filter.Cols = as.logical(param$Heatmap.Detailed.Prune.Cols))
            
          } else
            l3 = plotPVal(-PV/log(10), C, max.Count = as.numeric(param$PlotPrograms.MaxCount), 
                          MapType = param$Heatmap.Detailed.Type,
                          Filter.Rows = as.logical(param$Heatmap.Detailed.Prune.Rows),
                          Filter.Cols = as.logical(param$Heatmap.Detailed.Prune.Cols))
          
#          ggsave(plot=l3$plot.Q, #+guides(fill=FALSE), 
#                 filename = paste0(outputPlotPrograms,"-",sigName,"-heat.pdf"),width=7,height=10)
          ggsave(plot=l3$plot.Q, #+guides(fill=FALSE), 
                 filename = paste0(outputPlotPrograms,"-",sigName,"-heat.png"),width=7,height=10)
        }
      }
      if( param$PlotPrograms.IndividualBarChart ) {
        df = data.frame(Sample = names(LL), Signature = Counts.obs[,sigName])
        df = data.frame(Sample = names(LL), Signature = Counts.norms[,sigName])
        df$Sample <- factor(df$Sample, levels = rev(names(LL)))
        p = ggplot(df, aes(Sample, Signature)) + theme_classic()
        p = p + geom_col(color = NA, fill = "gray50")
#        p = p + geom_hline(yintercept = mean(Counts.expected[,sigName]), color = "black")
        p = p + geom_hline(yintercept = 0, color = "black")
        p = p + coord_flip()
              ggsave(plot=p, 
               filename = paste0(outputPlotPrograms,"-",sigName,"-bars.pdf"),width=7,height=10)
      }
    }
  }
}


cfChIP.plotEnrichments = function(LL, 
                                  outputPlotPrograms, 
                                  param= cfChIP.Params() ) {
  if( param$Verbose ) catn("Evaluating enrichments to ", outputPlotPrograms)
  
  EE = lapply(LL, function(dat) cfChIP.EvaluateProgramsHypG(dat, param = param, Write = FALSE) )
  
  PVals = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$pValue))))
  Counts.obs = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$ObservedCounts))))
  Counts.N = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$N))))
  Counts.K= do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$K))))
  Counts.Overlaps = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$Fraction))))
  
  rownames(PVals) = names(LL)
  colnames(PVals) = rownames(EE[[1]])
  dimnames(Counts.obs) = dimnames(PVals)
  dimnames(Counts.N) = dimnames(PVals)
  dimnames(Counts.K) = dimnames(PVals)
  dimnames(Counts.Overlaps) = dimnames(PVals)
  
  l = plotPVal(PVals, Counts.Overlaps, 
               max.Count = as.numeric(param$PlotEnrichments.MaxCount), 
               MapType = param$Heatmap.Type,
               Filter.Rows = as.logical(param$Heatmap.Prune.Rows),
               Filter.Cols = as.logical(param$Heatmap.Prune.Cols),
               colOrder = "sum" )
  
  #  catn(outputPlotPrograms)
  ggsave(plot=l$plot.Q + guides(fill=FALSE)+coord_equal(), 
         filename = paste0(outputPlotPrograms, ".pdf"),width=7,height=10)
  write.csv(l$QVals, file = paste0(outputPlotPrograms, "-qValues.csv") )
  write.csv(Counts.Overlaps, file = paste0(outputPlotPrograms, "-overlaps.csv") )
}

cfChIP.countQC = function(dat, filename, GR = NULL, param = cfChIP.Params()) {
  if( is.null(GR) )
    GR = param$TSS.windows
  
  dat.bed = cfChIP.GetRawData(filename, param)
  
  Widths = sapply(split(width(GR), GR$type), sum)

  if( !is.null(dat.bed$BED) ) {
    Counts = sapply(split(countOverlaps(query = GR, 
                                        subject = resize(dat.bed$BED, width=1, fix="center")), 
                          GR$type), 
                    sum)
    # total fragments sequenced
    Total = sum(dat$DupCount$Freq * as.numeric(dat$DupCount$RawBED.dups)) 
    Uniq_frags = length(dat.bed$BED)
  } else {
    Cov = fixCoverage(dat.bed$Cov)
    Counts = sapply(split(1:length(GR), GR$type), function(X) sum(sum(Cov[GR[X]])))
    Total =  sum(dat$Counts)
    Uniq_frags = sum(dat$Counts)
  }
  # estimation of library complexity
  freq = dat$DupCount$Freq
  names(freq) = dat$DupCount$RawBED.dups
  lambda1 = FitNonZeroPoisson(freq)
  Nuniq = sum(freq)
  Total1 = Nuniq /(1 - exp(-lambda1))
  p.mix = FitMixPoisson(freq)
  lambda.mix = round(p.mix$lambda * sum(unlist(lapply(1:length(p.mix$p), function (i) i*p.mix$p[[i]]))), digits = 3)
  Total.mix = round(Nuniq /(1 - MixPoisson(0, p.mix$lambda, p.mix$p)), digits = 0)
 
  
  Noise = dat$Background$genome
#  Counts = c(Counts, Total = Total)
  Widths = c(Widths, Total = sum(width(param$TSS.windows)))
  Background = Noise * Widths / 1000
  Signal = c(Counts, Total = Uniq_frags) - Background
  names(Signal) = paste("Signal", names(Background))
  names(Background) = paste("Background", names(Background))
  Counts = c(Counts, Background, Signal)
  Percents = formatC(100*Counts/Uniq_frags)
  names(Percents) = paste0("%", names(Counts))
  
  FullRecord = c( Total = Total,  Total.uniq = Uniq_frags, Total.uniq.est = unname(Total.mix), 
                  lambda.mix.poiss = unname(lambda.mix), Noise = formatC(Noise), round(Counts), Percents)
  
  # estimate yield
  
  cfDNA = as.numeric(param$QC.cfDNA.ng)
  n_cells = cfDNA * 1e-9 / as.numeric(param$QC.GenomeWeight)
  n_genomes = n_cells * 2
  
  PositiveType = param$QC.PositiveType 
  signal_windows = which(GR$type == PositiveType)
  
  exp_tot_nucs = n_genomes * 
    as.numeric(param$QC.GenomeLength) / as.numeric(param$QC.NucleosomeLength)
  exp_marked_nucs = n_genomes * as.numeric(param$QC.PositiveNucs)
  
  uniq_frags = Nuniq
  est_uniq_frags = as.numeric(Total.mix)
  est_seq_yield = uniq_frags / est_uniq_frags
  
  signal_rate = as.numeric(FullRecord["%Signal Total"]) * 1e-2 
  background_rate = as.numeric(FullRecord["%Background Total"]) * 1e-2 
  signal_background_rate = as.numeric(FullRecord[paste("%Background", PositiveType)]) * 1e-2 
  
  est_signal_frags = uniq_frags * signal_rate
  est_background_frags = uniq_frags * background_rate
  global_yield = (est_signal_frags / exp_marked_nucs)
  global_bg_yield = (est_background_frags / exp_tot_nucs)
#  total_yield = (global_yield + non_specific_yield)
  
  #right now relevant only to TSS's
#  WinHeights = dat$Heights[TSS.windows$name %in% Genes[GeneDescription$H3K4me3.HouseKeeping]]
  
  I = GeneDescription[names(GeneLength),"H3K4me3.HouseKeeping"]
  I[is.na(I)] = FALSE
  I = I & GeneLength < quantile(GeneLength[I], 0.99)
  WinHeights = dat$GeneHeights[I]/dat$GeneCNV[I]

  if(0) {
    m = mean(WinHeights)    
    v = var(WinHeights)
    s = m**2/(v-m)
    qs = ppoints(length(WinHeights))
    df = data.frame( x = sort(WinHeights), 
                     y = (1:length(WinHeights))/length(WinHeights))
    df$y.norm = pnorm(df$x, mean = m, sd = sqrt(v))
    df$y.nb = pnbinom(df$x, size = s, mu = m)
    df$q.norm = qnorm(qs, mean = m, sd = sqrt(v))
    df$q.nb = qnbinom(qs, size = s, mu = m)
    
    ggplot(df, aes(x,y))+geom_line()+geom_line(aes(y=y.norm), color="blue")+geom_line(aes(y=y.nb), color="red")
    qplot(df$q.norm, df$x)
    qplot(df$q.nb, df$x)
    qplot(GeneLength[I], WinHeights)+geom_smooth(color="red")+geom_density2d()
    
  }
  
  high_cov = as.numeric(quantile(WinHeights, as.numeric(param$QC.QuantileCutoff)))
  local_yield = signal_rate * high_cov / n_genomes

  
  FullRecord = c(FullRecord, 
                 Global.signal.yield = formatC(global_yield,digits = 6),
                 Local.signal.yield = formatC(local_yield,digits = 6),
                 BG.yield = formatC(global_bg_yield,digits = 6),
                 Global.SNR = formatC(global_yield/global_bg_yield, digits = 4),
                 Local.SNR = formatC(local_yield/global_bg_yield, digits = 4),
                 Seq.factor = formatC(1/est_seq_yield,digits = 6)
  )
                 
  if( !is.null(param$QCFields) ) {
    SelectFields = strsplit(param$QCFields, ";")[[1]]
    SelectFields = SelectFields[SelectFields %in% names(FullRecord)]
    FullRecord = FullRecord[SelectFields]
  }
  return(FullRecord)
}


cfChIP.EstimateMeanBasis = function(SigCounts, SigBackground, QQnorm) {
  SigDiff = SigCounts - SigBackground
  SigDiff[SigDiff < 0] = 0
  
  Sig.avg = rowMeans(SigDiff)    
  I = Sig.avg > 0
  J = QQnorm > 0
  Sig.avg[I] = rowMeans(t(t(SigDiff[I,J])*QQnorm[J]))
  
  Sig.avg
}

cfChIP.EstimateMeanVarianceBasis = function(SigCounts, SigBackground, QQnorm) {
  SigDiff = SigCounts - SigBackground
  SigDiff[SigDiff < 0] = 0
  
  Sig.avg = rowMeans(SigDiff)    
  Sig.var = rowMeans((SigDiff - Sig.avg)**2)
  names(Sig.avg) = rownames(SigCounts)
  names(Sig.var) = rownames(SigCounts)
  I = which(Sig.avg > 0)
  J = QQnorm > 0
  if( length(I) > 0) {
    Sig.est = sapply(I, function(w) EstimateMeanVar(SigCounts[w,J], 1/QQnorm[J], SigBackground[w,J]))
    Sig.avg[I] = Sig.est["meanP",]
    Sig.var[I] = Sig.est["varP",]
  }
  list( avg = Sig.avg, var = Sig.var )
}

cfChIP.EstimateMeanVarianceWinSig = function(LL, WinSig) {

  WinCounts = do.call("cbind", lapply(LL, function(l) sapply(WinSig, function(w) sum(l$Counts[w]))))
  WinBackground = do.call("cbind", lapply(LL, function(l) sapply(WinSig, function(w) sum(l$WinBackground[w]))))
  QQnorm = sapply(LL, function(l) l$QQNorm)
  names(QQnorm) = colnames(WinCounts)

  cfChIP.EstimateMeanVarianceBasis(WinCounts, WinBackground, QQnorm)
}

cfChIP.EstimateMeanVarianceGeneSig = function(LL, GeneSig ) {
  GeneCounts = do.call("cbind", lapply(LL, function(l) sapply(GeneSig, function(w) sum(l$GeneCounts[w]))))
  GeneBackground = do.call("cbind", lapply(LL, function(l) sapply(GeneSig, function(w) sum(l$GeneBackground[w]))))
  QQnorm = sapply(LL, function(l) l$QQNorm)
  names(QQnorm) = colnames(GeneCounts)
  
  cfChIP.EstimateMeanVarianceBasis(GeneCounts, GeneBackground, QQnorm)
}

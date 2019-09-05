source(paste0(SourceDIR, "QQNormalize.R")) 
source(paste0(SourceDIR, "MetaPlot.R")) 
source(paste0(SourceDIR,"cfChIP-util.R"))
source(paste0(SourceDIR, "PlotPVal.R"))
source(paste0(SourceDIR, "YieldEstimation-Functions.R"))

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
    Signatures = Win.Sig,
    MinFragLen = 50,
    MaxFragLen = 500,
    Heatmap.Split = TRUE,
    Heatmap.MaxCount = 15,
    Scatter.NormalizeByLength = FALSE,
    Scatter.RefSeqOnly = FALSE,
    Scatter.MarkGenePrograms = NULL,
    Scatter.UpperQuantile = 0.995,
    Scatter.LimitUnits = 2,
    Programs = Gene.Programs,
    PlotPrograms.MaxCount = 40,
    PlotEnrichments.MaxCount = 50,
    PlotPrograms.IndividualHeatmap = TRUE,
    PlotPrograms.IndividualBarChart = TRUE,
    PlotPrograms.IndividualCSV = TRUE,
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
    MetaEnhancer.BGColor = "black"
    #    BackgroundModel = NULL
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

cfChIP.ProcessFile <- function( filename = NULL,
                               dat = NULL,
                               param = cfChIP.Params(),
                               Force = FALSE,
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
                 Background = NULL,
                 GeneCounts = NULL,
                 GeneBackground = NULL, 
                 Counts.QQnorm = NULL,
                 GeneCounts.QQnorm = NULL, 
                 QQNorm = 1,
                 OverExpressedGenes = NULL)
      Change = TRUE
    }
  }
  
  # Get BED reads into a GenomicRanges object
  if( is.null(dat$Counts) ) {
    FileType = NA
    if( grepl(".bed$", filename) | grepl(".tagAlign$", filename))
      FileType = "BED"
    if( grepl(".bw$", filename) )
      FileType = "BW"

    if( is.na(FileType) ) {
      if( file.exists(paste0(filename, ".bed"))) {
        filename = paste0(filename, ".bed")
        FileType = "BED"
      } else
        if( file.exists(paste0(filename, ".tagAlign"))) {
          filename = paste0(filename, ".tagAlign")
          FileType = "BED"
        } else
        if( file.exists(paste0(filename, ".bw"))) {
          filename = paste0(filename, ".bw")
          FileType = "BW"
        }   
    }

    if( is.na(FileType ) ) {
      catn(Name, ": Error, cannot determine file type of ",filename)
      return(NULL)
    }
     
    if( FileType == "BED") {     
      if( param$Verbose ) catn(Name, ": Reading BED file", filename)

      RawBED = import(filename, format = "BED")
      # remove long/short fragments and non-unique copies
      
      # check for single end reads
      if( max(width(RawBED)) <  param$MinFragLen) {
        RawBED = resize(RawBED, width = 166)
      } else 
        RawBED = RawBED[width(RawBED) <= param$MaxFragLen & width(RawBED) > param$MinFragLen]
      
      dat$BED = unique(RawBED)
      
      # count duplicate segments
      RawBED.dups <- countOverlaps(dat$BED, RawBED, type = "equal")
      dat$DupCount = as.data.frame(table(RawBED.dups))
    } 
    if( FileType == "BW" ) {     
      if( param$Verbose ) catn(Name, ": Reading BigWig file", filename)
      dat$BW = import(filename)
 
    } 
    Change = TRUE
  }
  
  #enforce dependencies
  param$Normalize = param$Normalize || param$OverExpressedGenes
  param$GeneCounts = param$GeneCounts || param$Normalize
  param$GeneBackground = param$GeneBackground || param$Normalize
  param$Background = param$Background || param$GeneBackground
  
  if( Force ) {
    dat$Counts = NULL
    dat$Cov = NULL
    dat$Background = NULL
    dat$GeneCounts = NULL
    dat$GeneBackground = NULL 
    dat$Counts.QQnorm = NULL
    dat$GeneCounts.QQnorm = NULL 
    dat$QQNorm = 1
  }
  
  #Counts
  if( is.null(dat$Counts) ) {
    if(!is.null(dat$BED)) {
      if( param$Verbose ) catn(Name, ": counting fragment overlap")
      dat$Counts = countOverlaps(query = param$TSS.windows, 
                                 subject = resize(dat$BED, width=100, fix="center"))
      dat$Cov = coverage(dat$BED)
    } else
      if( !is.null(dat$BW)) {
        if( param$Verbose ) catn(Name, ": counting BigWig overlap")
        dat$Cov = coverage(dat$BW, weight="score")
        dat$Counts = rep(0,length(param$TSS.windows))
        if(0) {
          Score = rep(0, length(param$TSS.windows))
          Subs = dat$BW
          Q = findOverlaps(Subs,param$TSS.windows)
          X = Subs[queryHits(Q)]
          Y = subjectHits(Q)
          Z = width(X) * X$score
          O = order(Y)
          S = sapply(split(Z[O],Y[O]), sum)
          uY = unique(Y[O])
          Score[uY] = S
          dat$Counts = round(Score/200)
        }
        if( !( "chrY" %in% names(dat$Cov)) ) 
          dat$Cov[["chrY"]] = Rle(0,seqlengths(param$TSS.windows)["chrY"])
        if(0) {
          GR = GRanges( seqnames = chrom(param$TSS.windows),
                        ranges = ranges(param$TSS.windows),
                        seqinfo =  seqinfo(param$TSS.windows)[names(dat$Cov)])
          ns = c(seq(1, length(GR)-1, by=500000),length(GR))
          for( i in 1:(length(ns)-1)) {
            catn(i)
            dat$Counts[ns[i]:ns[i+1]] = dat$Cov[GR[ns[i]:ns[i+1]]]
          }
        }
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
      } else {
        catn(Name, ": error! count compute counts")
        return(dat)
      }
    dat$Cov = fixCoverage(dat$Cov)
  }
  
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
    names(dat$GeneCounts) = Genes
    dat$GeneCounts.QQnorm = NULL
    Change = TRUE
  }
  
  # GeneBackground
  if( param$GeneBackground && is.null(dat$GeneBackground)) {
    if( param$Verbose ) catn(Name, ": Computing gene background")
    mu = dat$Background
    Z = getMultiBackgroundEstimate(mu,param$GeneWindows)
    dat$GeneBackground =  ComputeGeneCounts(Z,param$Win2Gene)
    names(dat$GeneBackground) = Genes
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
    bg = rep(dat$Background$genome, length(param$TSS.windows))
    rl = chrom(param$TSS.windows)
    Vs = as.character(runValue(rl))
    for( i in 1:nrun(rl) )
      bg[start(rl)[i]:end(rl)[i]] = dat$Background$chr[Vs[i]]
    bg = bg*width(param$TSS.windows)/1000
    
    dat$Counts.QQnorm = pmax(dat$Counts - bg, 0)*dat$QQNorm
    Change=TRUE
  }
  
  #estimate over expressed genes
  if( param$OverExpressedGenes && is.null(dat$OverExpressedGenes)) {
    if( param$Verbose ) catn(Name, ": Computing overexpressed genes")
    
    X = dat$GeneCounts.QQnorm
    Y = dat$GeneCounts
    logX = log2(X+1)
    logH = log2(param$NormRef+1)
    Diff.log = logX - logH
    
    Lam = param$NormRef / dat$QQNorm + dat$GeneBackground 
    Lam[is.na(Lam)] = 100
    Pv = -ppois(Y,Lam, lower.tail = FALSE, log.p = T)*log(10)
    Qv = -log10(p.adjust(exp(-Pv/log(10)),method="fdr"))

    dat$OverExpressedGenes = data.frame(healthy = logH, 
                                        sample = logX, 
                                        obs = Y, 
#                                        bg = dat$GeneBackground,
                                        exp = Lam, 
                                        qvalue = Qv, pvalue = Pv, X = X, H = param$NormRef, stringsAsFactors = FALSE)
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

cfChIP.MetaPlot = function(dat, Dir = NULL,  param = cfChIP.Params(), Force = FALSE) {
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

cfChIP.WriteTrack = function(dat, Dir = NULL, Force = FALSE) {
  fname = cfChIP.BuildOutputName(dat, Dir, Suffix = ".bw")
  if( Force || !file.exists(fname)) {
    XX.gr = GRanges(dat$Cov*dat$QQNorm, seqinfo = genome.seqinfo)
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

cfChIP.EvaluateSignatures = function( dat, Dir = NULL, param = cfChIP.Params(), Write = TRUE, Force = FALSE ) {
 
  if( Write ) 
    fname = cfChIP.BuildOutputName(dat, Dir, "Signatures", ".csv")
    
  if(is.null(param$Signatures)) {
    message("Error: missing signatures")
    return()
  }
  if( Force || !Write || !file.exists(fname) )  {
    Sig = param$Signatures
    Evals = lapply(Sig,function(s) cfChIP.EvalSig(s,dat))
    PVals = -sapply(names(Evals),function(s) Evals[[s]][1])
    QVals = -log10(p.adjust(exp(-PVals),method="fdr"))
    Fg = sapply(names(Evals),function(s) Evals[[s]][2])
    Bg = sapply(names(Evals),function(s) Evals[[s]][3])
    names(Fg) = names(Sig)
    names(Bg) = names(Sig)
    
    Counts = Fg - Bg
    Counts[Counts<0] = 0
    Counts.norms = Counts * dat$QQNorm
    Sig.Width = sapply(Sig, function(sig) sum(width(param$TSS.windows[sig])))/1000
    Counts.norms = t(t(Counts.norms) / Sig.Width)
    
    df = data.frame(ObservedCounts = Fg, Background = formatC(Bg), NormalizedCounts = formatC(Counts.norms),
                    pValue = formatC(PVals/log(10)), qValue = formatC(QVals) )
    
    if( Write ) 
      write.csv(df[order(PVals, decreasing = TRUE),], file = fname, quote = FALSE)
    return(df)
  }
}

cfChIP.EvaluatePrograms = function( dat, Dir = NULL, param = cfChIP.Params(), Write = TRUE, Force = FALSE ) {

  if( Write ) 
    fname = cfChIP.BuildOutputName(dat, Dir, "Programs", ".csv")

  if( Force || !Write || !file.exists(fname) )  {  
    catn(dat$Name)
    O = dat$OverExpressedGenes
    Sig = param$Programs
    Fg = sapply(Sig, function(sig) sum(O[sig, "obs"]))
    Bg = sapply(Sig, function(sig) sum(dat$GeneBackground[sig]))
    Exp = sapply(Sig, function(sig) sum(O[sig, "exp"]))
    Counts.norms = sapply(Sig, function(sig) sum(2**(O[sig, "sample"])-1))
    Healthy.norms = sapply(Sig, function(sig) sum(2**(O[sig, "healthy"])-1))
#    PVals = -ppois(Fg, Exp, lower.tail = FALSE, log.p = T)
    PVals = -log(sapply(1:length(Fg), function(i) poisson.test(Fg[i],Exp[i])$p.value))

    QVals = -log10(p.adjust(exp(-PVals),method="fdr"))
    
    Sig.Width = sapply(Sig, length)
    Counts.norms = t(t(Counts.norms) / Sig.Width)
    Healthy.norms = t(t(Healthy.norms) / Sig.Width)
    df = data.frame(ObservedCounts = Fg, 
                    Expected = Exp, 
                    Background = Bg,
                    NormalizedCounts = formatC(Counts.norms), 
                    NormalizedExpected = formatC(Healthy.norms), 
                    pValue = formatC(PVals/log(10)), qValue = formatC(QVals), stringsAsFactors = FALSE)
    
    if( Write ) 
    write.csv(df[order(PVals, decreasing = TRUE),], file = fname, quote = FALSE)
    
    return(df)
  }
}

Max.Qv = 20

cfChIP.OverExpressedGenes = function( dat, Dir = NULL, param = cfChIP.Params(), Force = FALSE ) {
  fname = cfChIP.BuildOutputName(dat, Dir, "DiffGenes", ".csv")
  fname2 = cfChIP.BuildOutputName(dat, Dir, "DiffScatter", ".pdf")
  
  if( Force || !file.exists(fname2)) {
    O = dat$OverExpressedGenes
    O = O[which(O$qvalue > 3 & (O$sample - O$healthy) > 2 & rownames(O) %in% Genes.notexcluded),]
    
    if( nrow(O) > 1 ) {
      if( !is.null(GeneDescription)) {
        O = cbind(O, GeneDescription[rownames(O),])
        df = data.frame(Obs = O$X, 
                        Ref = O$H, 
                        FoldChange = O$sample-O$healthy, 
                        Pvalues = O$pvalue, QValues = O$qvalue,
                        Tissue = O$H3K4me3.Tissue, HighExpressed = O$GTEX.HighExpressed, Expressed = O$GTEX.Expressed,
                        Suspect = O$Suspect, Browser=GenesBrowserWin( rownames(O)),
                        Description = O$Description, stringsAsFactors = FALSE)
        rownames(df) = rownames(O)
        df = df[order(df$FoldChange, decreasing = TRUE),]
        colnames(df) =  c("Observed", "Healthy", "FoldChange (log2)", "p-Value (-log10)",  "Q-Value (-log10)", 
                          "Tissue (Atlas)", "High Expressed (GTEX)", "Expressed (GTEX)",
                          "Suspect?", "Browser window", "Description")
      } else {
        df = data.frame(Obs = O$X, Ref = O$H, FoldChange = O$sample-O$healthy, 
                        Pvalues = O$pvalue, QValues = O$qvalue,
                        Browser=GenesBrowserWin( rownames(O), stringsAsFactors = FALSE)
        )
        df = df[order(df$FoldChange, decreasing = TRUE),]
        colnames(df) =  c("Observed", "Healthy", "FoldChange (log2)", "p-Value (-log10)",  "Q-Value (-log10)", 
                          "Browser window")
      }
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
    
    Significant = df$qvalue>3 & df$sample-df$healthy > 2
    
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
    p = p + geom_text(data=data.frame(x = Label.x,y=mVal, label = paste(sum(Significant %in% Genes.notexcluded),"significant genes")),
                      aes(x,y,label=label), color="black")
    ggsave(p, filename = fname2, width=6, height=6)
  }
}


cfChIP.EvaluateProgramsHypG = function( dat, Dir = NULL, param = cfChIP.Params(), Write = TRUE, Force = FALSE) {
  if( Write ) 
    fname = cfChIP.BuildOutputName(dat, Dir, "Programs-HypG", ".csv")
  
  if( Force || !Write || !file.exists(fname) )  {  
    O = dat$OverExpressedGenes
    UpGenes = Genes[which(O$qvalue > 3 & (O$sample - O$healthy) > 2)]
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

cfChIP.WriteEnrichR = function( dat, Dir = OutputDIR, Force = FALSE ) {
  catn("EnrichR not implemented yet")
}

cfChIP.plotSignatures = function(LL, outputPlotSignatures, param= cfChIP.Params() ) {
  if( param$Verbose ) catn("Evaluating signatures to ", outputPlotSignatures)
  
  EE = lapply(LL, function(dat) cfChIP.EvaluateSignatures(dat, param = param, Write = FALSE) )
  
  PVals = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$pValue))))
  Counts.norms = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$NormalizedCounts))))
  rownames(PVals) = names(LL)
  rownames(Counts.norms) = names(LL)
  colnames(PVals) = rownames(EE[[1]])
  colnames(Counts.norms) = colnames(PVals) 
  PVals[Counts.norms < 0.25] = 0
  
  l = plotPVal(-PVals, Counts.norms, 
               max.Count = as.numeric(param$Heatmap.MaxCount), 
               splitMap = as.logical(param$Heatmap.Split))
  colOrder = order(colSums(l$QVals))
  colOrder = order(colSums(Counts.norms))
  
  l = plotPVal(-PVals, Counts.norms, 
               max.Count = as.numeric(param$Heatmap.MaxCount),
               splitMap = as.logical(param$Heatmap.Split),
               colOrder = colOrder, 
               rowOrder = rev(1:length(LL))
               )
  
#  ggsave(plot=l$plot.Q + guides(fill=FALSE)+coord_equal(), filename = paste0(outputPlotSignatures, ".pdf"),width=7,height=10)
  ggsave(plot=l$plot.Q + coord_equal(), filename = paste0(outputPlotSignatures, ".pdf"),width=7,height=10)
  write.csv(PVals, file = paste0(outputPlotSignatures, "-pValues.csv") )
  write.csv(Counts.norms, file = paste0(outputPlotSignatures, "-Counts.csv") )
  
  Sig = param$Signatures
  for( sigName in names(Sig) ) {
    W = Sig[[sigName]]
    Bg = sapply(LL, function(dat) getMultiBackground(dat$Background, W)[1,])
    Fg = sapply(LL, function(dat) dat$Counts[W])
    Cg = sapply(LL, function(dat) dat$Counts.QQnorm[W])
    rownames(Bg) = W
    rownames(Fg) = W
    rownames(Cg) = W
    PV = matrix(ppois(Fg,Bg,lower.tail = FALSE, log.p = TRUE),
                nc = ncol(Fg), nr = nrow(Fg), dimnames = dimnames(Fg))
    rowOrder = order(rowSums(Cg))
    rownames(Cg) = paste(rownames(Cg), TSS.windows[W]$type, TSS.windows[W]$name)
    rownames(PV) = rownames(Cg)  
    Cg = rbind(Cg[rowOrder,], signature = l$Counts[,sigName])
    PV = rbind(PV[rowOrder,], signature = l$PVals[,sigName])
    ll = plotPVal(PV,Cg, splitMap = param$Heatmap.Split)
    ggsave(plot=ll$plot.Q+guides(fill=FALSE)+coord_equal(), filename = paste0(outputPlotSignatures,"-",sigName,".pdf"),width=7,height=10)
  }
}


cfChIP.plotPrograms = function(LL, 
                               outputPlotPrograms, 
                               param= cfChIP.Params() ) {
  if( param$Verbose ) catn("Evaluating programs to ", outputPlotPrograms)
  
  EE = lapply(LL, function(dat) cfChIP.EvaluatePrograms(dat, param = param, Write = FALSE) )
  
  PVals = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$pValue))))
  Counts.obs = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$NormalizedCounts))))
  Counts.expected = do.call(rbind, lapply(EE, function(e) as.numeric(as.character(e$NormalizedExpected))))
  Counts.norms = Counts.obs - Counts.expected
  rownames(PVals) = names(LL)
  rownames(Counts.norms) = names(LL)
  rownames(Counts.obs) = names(LL)
  rownames(Counts.expected) = names(LL)
  colnames(PVals) = rownames(EE[[1]])
  colnames(Counts.norms) = colnames(PVals) 
  colnames(Counts.obs) = colnames(PVals) 
  colnames(Counts.expected) = colnames(PVals) 
  
  l = plotPVal(-PVals, Counts.norms, 
               max.Count = as.numeric(param$PlotPrograms.MaxCount),
               splitMap = as.logical(param$Heatmap.Split))
  colOrder = order(colSums(l$QVals))
  colOrder = hclust(dist(t(Counts.norms)))$order
  l = plotPVal(-PVals, Counts.norms, 
               max.Count = as.numeric(param$PlotPrograms.MaxCount), 
               splitMap = as.logical(param$Heatmap.Split),
               colOrder = colOrder, rowOrder = rev(1:length(LL)) )
  
#  catn(outputPlotPrograms)
  ggsave(plot=l$plot.Q + guides(fill=FALSE)+coord_equal(), 
         filename = paste0(outputPlotPrograms, ".pdf"),width=7,height=10)
  write.csv(l$QVals, file = paste0(outputPlotPrograms, "-qValues.csv") )
  write.csv(Counts.norms, file = paste0(outputPlotPrograms, "-Counts.csv") )
  
  catn("foobar")
  if( param$PlotPrograms.IndividualHeatmap ||
      param$PlotPrograms.IndividualCSV ||
      param$PlotPrograms.IndividualBarChart ) {
    HH  = param$Programs
    for( sigName in names(HH)) {
      if( param$PlotPrograms.IndividualHeatmap || param$PlotPrograms.IndividualCSV ) {
        G = HH[[sigName]]
        PV = sapply(LL, function(dat) dat$OverExpressedGenes[G, "pvalue"])
        rownames(PV) = G
        C = sapply(LL, function(dat) dat$GeneCounts.QQnorm[G] - param$NormRef[G])
        rowOrder = rev(order(rowSums(C)))
        tryCatch({        
          rowOrder = hclust(dist(C))$order
          })
        PV = rbind(PV[rev(rowOrder),],
                   program = t(l$PVals[colnames(PV),sigName]))
        C = rbind(C[rev(rowOrder),colnames(PV)],
                  program = t(l$Counts[colnames(PV),sigName]))
        if( PlotPrograms.IndividualCSV ) {
          write.csv(PV, file = paste0(outputPlotPrograms,"-",sigName,"-pvalues.CSV"), quote = FALSE)
          write.csv(C, file = paste0(outputPlotPrograms,"-",sigName,"-counts.CSV"), quote = FALSE)
        }
        if( PlotPrograms.IndividualHeatmap ) {
          l3 = plotPVal(PV, C, max.Count = as.numeric(param$PlotPrograms.MaxCount), 
                        splitMap = as.logical(param$Heatmap.Split))
          ggsave(plot=l3$plot.Q+guides(fill=FALSE), 
                 filename = paste0(outputPlotPrograms,"-",sigName,"-heat.pdf"),width=7,height=10)
        }
      }
      if( param$PlotPrograms.IndividualBarChart ) {
        df = data.frame(Sample = names(LL), Signature = Counts.obs[,sigName])
        df = data.frame(Sample = names(LL), Signature = Counts.norms[,sigName])
        df$Sample <- factor(df$Sample, levels = names(LL))
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
  
  
  if(1) {
    l = plotPVal(-PVals, Counts.Overlaps, 
                 max.Count = as.numeric(param$PlotEnrichments.MaxCount),
                 splitMap = as.logical(param$Heatmap.Split))
    colOrder = order(colSums(l$QVals))
  } else 
    colOrder = hclust(dist(t(Counts.Overlaps)))$order
  l = plotPVal(-PVals, Counts.Overlaps, 
               max.Count = as.numeric(param$PlotEnrichments.MaxCount), 
               splitMap = as.logical(param$Heatmap.Split),
               colOrder = colOrder, rowOrder = rev(1:length(LL)) )
  
  #  catn(outputPlotPrograms)
  ggsave(plot=l$plot.Q + guides(fill=FALSE)+coord_equal(), 
         filename = paste0(outputPlotPrograms, ".pdf"),width=7,height=10)
  write.csv(l$QVals, file = paste0(outputPlotPrograms, "-qValues.csv") )
  write.csv(Counts.Overlaps, file = paste0(outputPlotPrograms, "-overlaps.csv") )
}

cfChIP.countQC = function(dat, GR = NULL, param = cfChIP.Params()) {
  if( is.null(GR) )
    GR = param$TSS.windows
  
  Widths = sapply(split(width(GR), GR$type), sum)
  
  if( !is.null(dat$BED) ) {
    Counts = sapply(split(countOverlaps(query = GR, 
                                        subject = resize(dat$BED, width=1, fix="center")), 
                          GR$type), 
                    sum)
    # total fragments sequenced
    Total = sum(dat$DupCount$Freq * as.numeric(dat$DupCount$RawBED.dups)) 
    Uniq_frags = length(dat$BED)
  } else {
    Cov = fixCoverage(dat$Cov)
    Counts = sapply(split(1:length(GR), GR$type), function(X) sum(sum(Cov[GR[X]])))
    Total = sum(dat$DupCount$Freq * as.numeric(dat$DupCount$RawBED.dups)) 
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
  return(c( Total = Total,  Total.uniq = Uniq_frags, Total.uniq.est = unname(Total.mix), 
            lambda.mix.poiss = unname(lambda.mix), Noise = formatC(Noise), round(Counts), Percents))
}

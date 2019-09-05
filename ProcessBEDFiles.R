#!/usr/local/bin/Rscript --vanilla
#

catn = function(...) { cat(...,"\n") }


initial.options <- commandArgs(trailingOnly = FALSE)
trailing.options = commandArgs(trailingOnly = TRUE)
#Find where are the scripts at

if( !any(grepl("--interactive", initial.options)) ) {
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  SourceDIR <- paste0(dirname(script.name),"/")
  DataDir = paste0(getwd(), "/")
  ANNOTDIR = DataDir
  DevelopmentMode = FALSE
} else {
  # development mode (inside RStudio) we can set options using these variables
  DevelopmentMode = TRUE
  Mod = "H3K4me3-scer/"
  SourceDIR = "~/GoogleDrive/Src/BloodChIPAnalysis/"
  # SourceDIR = "~/GoogleDrive/BloodChIPAnalysis/"
  ANNOTDIR = paste0(SourceDIR, "SetupFiles/", Mod)
  DataDir = paste0("~/GoogleDrive/BloodChIP/Analysis/Samples/", Mod)
  BedDir = paste0("~/GoogleDrive/BloodChIP/Analysis/BED/", Mod)
  OutputDir = paste0("~/GoogleDrive/BloodChIP/Analysis/Output/", Mod)
  SetupDIR  = ANNOTDIR
  trailing.options = paste(
    #                   "-r /Users/nir/GoogleDrive/BloodChIP/Analysis-Paper",
    # "-r /Users/nir/GoogleDrive/BloodChIP/Analysis",
    # "-r /Users/nir/GoogleDrive/BloodChIP/Analysis-Blueprint",
    "-r ~/GoogleDrive/BloodChIP/Analysis",
    "-A",
    "-m H3K4me3-scer", "PS003_15_24hr_RT_No_PI.K4me3-9-Aug-2019_S16.yeast.tagAlign"
    #  "-o /Users/nir/GoogleDrive/BloodChIP/Analysis-Paper/Figures/H3K36me3",
    #  "--geneprograms=/Users/nir/GoogleDrive/BloodChIP/Analysis-Paper/Figures/H3K36me3/TCGA-short-gene-sig.csv",
    #  "--config=/Users/nir/GoogleDrive/BloodChIP/Analysis-Paper/Figures/H3K36me3/config.csv",
    #    "--plotsignatures=FigS4D E062 E066 E034",
  )
  trailing.options = strsplit(trailing.options, " ")[[1]]
  catn(trailing.options)
}


OrigDir = getwd()

extendDir <- function(x) {
  a = substr(x,1,1)
  if( a != "/" && a != "~")
    x = paste0(OrigDir,"/",x)
  a = substr(x, nchar(x), nchar(x))
  if( a != "/")
    x = paste0(x,"/")
  x
}


outputFile = function(x) {
  a = substr(x,1,1)
  if( a != "/" && a != "~")
    x = paste0(extendDir(OutputDir), x)
  x
}

suppressMessages(library("optparse"))
option_list = list(
# directory 
  make_option(c("-d", "--datadir"), type = "character", default = NULL,
              help="data directory"),
  make_option(c("-b", "--BEDdir"), type="character", default=NULL,
              help="location of BED files" ),
  make_option(c("-a", "--annotationdir"), type="character", default=NULL,
              help="location of annotation files" ),
  make_option(c("-t", "--trackdir"), type="character", default=NULL,
              help="location of genome browser track files" ),
  make_option(c("-o", "--outputdir"), type="character", default=NULL,
              help="location of output files" ),
  make_option(c("-m", "--mod"), type="character", default="H3K4me3",
              help="name of modification, used in combination with -r" ),
  make_option(c("-r", "--root"), type="character", default=NULL,
              help="name of root directory, input/output is assumed to be in Root/XXX/mod" ),
  
  
# basic processing   
  make_option(c("-B", "--background"),  action = "store_true", type="logical", default=FALSE, 
              help="Compute background model"),
  make_option(c("-C", "--count"),  action = "store_true", type="logical", default=FALSE, 
              help = "Count gene coverage"),
  make_option(c("-N", "--normalize"),  action = "store_true", type="logical", default=FALSE, 
              help = "Normalize gene counts"),
  make_option(c("-O", "--overexpressed"),  action = "store_true", type="logical", default=FALSE, 
              help = "Evaluate overexpressed genes"),
  make_option(c("-A", "--all"),  action = "store_true", type="logical", default=FALSE, 
            help="Compute all relevant data for each file (equivalent to -BCN)"),
  make_option(c("-F", "--force"),  action = "store_true", type="logical", default=FALSE, 
            help="Forced recomputing"),

# action per file 
  make_option(c("-T", "--tracks"),  action = "store_true", type="logical", default=FALSE, 
            help = "Create normalized genome browser tracks"),
  make_option(c("-M", "--meta"),  action = "store_true", type="logical", default=FALSE, 
            help = "Generate meta plots"),
  make_option(c("-E", "--enrichR"),  action = "store_true", type="logical", default=FALSE, 
            help = "Evaluate enrichments of overexpressed genes"),
  make_option(c("-S", "--signatures"),  action = "store_true", type="logical", default=FALSE, 
            help = "Evaluate signatures of cell types"),
  make_option(c("-P", "--programs"),  action = "store_true", type="logical", default=FALSE, 
            help = "Evaluate expression programs"),
  make_option("--backgroundplot",type="logical", default=FALSE, action="store_true", 
              help = "Plot background coverage in output directory"),

# global output
  make_option(c("--consensus"),  type="character", default=NULL, 
              help = "Compute consensus (avg after re-normalization) and write to specified file"),
  make_option("--genecounts",  type="character", default=NULL, 
            help = "output table of gene counts to specified file"),
  make_option("--genebackground",  type="character", default=NULL, 
            help = "output table of gene background to specified file"),
  make_option("--wincounts",  type="character", default=NULL, 
            help = "output table of window counts to specified file"),
  make_option("--normcounts",  type="character", default=NULL, 
            help = "output table of normalized gene counts to specified file"),
  make_option("--plotsignatures",  type="character", default=NULL, 
            help = "output plot of signature enrichment to specified file"),
  make_option("--plotprograms",  type="character", default=NULL, 
            help = "output plot of program enrichment to specified file"),
  make_option("--plotenrichments",  type="character", default=NULL, 
            help = "output plot of program enrichment (HyperGeometric) to specified file"),
  make_option("--export",  type="character", default=NULL, 
            help = "export data to RDS file"),
  make_option("--QC", type = "character", default = NULL,
              help = "prepare QC report to file"),
  make_option("--cluster", type = "character", default = NULL,
            help = "prepare file for Cluster program"),


# auxilary files (instead of default)
  make_option("--windowsignatures",  type="character", default=NULL, 
              help = "Provide window siganture file (instead of Win-sig.csv)"),
  make_option("--geneprograms",  type="character", default=NULL, 
            help = "Provide gene gene programs file (instead of Gene-sig.csv)"),
  make_option("--metagenes",  type="character", default=NULL, 
            help = "Provide gene annotation file for meta plots (instead of meta-gene.bed)"),
  make_option("--metaenhancers",  type="character", default=NULL, 
            help = "Provide enhancer annotation file for meta plots (instead of meta-enhancer.bed)"),
  make_option("--config", type="character", default = NULL,
            help = "Specify an alternative configuration file")
); 

opt_parser = OptionParser(option_list=option_list);
if( exists("trailing.options") ) {
opt = parse_args(opt_parser, args = trailing.options, positional_arguments = TRUE)
} else
  opt = parse_args(opt_parser, positional_arguments = TRUE)

Files = opt$args

if(!DevelopmentMode && length(Files) == 0)
  parse_args(opt_parser,args = "-h", print_help_and_exit = TRUE)

catn("Initializing")
suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))

TargetMod = opt$options$mod
if( !is.null(opt$options$root)) {
  RootDir =  extendDir(opt$options$root)
  DataDir = paste0(RootDir, "Samples/", TargetMod, "/")
  BedDir = paste0(RootDir, "BED/", TargetMod, "/")
  TracksDIR = paste0(RootDir, "Tracks/", TargetMod, "/")
  OutputDir = paste0(RootDir, "Output/", TargetMod, "/")
} 
ANNOTDIR = paste0(SourceDIR, "SetupFiles/", TargetMod, "/")


if( !is.null(opt$options$datadir) )
  DataDir = extendDir(opt$options$datadir)

DIR = DataDir

if( !is.null(opt$options$annotationdir) )
  ANNOTDIR = extendDir(opt$options$annotationdir)

if( !exists("BedDir") ) 
  BedDir = DataDir


if( !is.null(opt$options$BEDdir) )
  BedDir = extendDir(opt$options$BEDdir)

if( !exists("TracksDIR") ) 
  TracksDIR = DataDir
if( !is.null(opt$options$trackdir) )
  TracksDIR = extendDir(opt$options$trackdir)

if( !is.null(opt$options$outputdir)) {
  OutputDir = extendDir(opt$options$outputdir)
} else
  if( !exists("OutputDir") ) {
    if( !DevelopmentMode  ) {
      OutputDir = OrigDir
    } else
      OutputDir = paste0("~/Data/BloodChIP/Output/", Mod)
  }

MetaDIR = OutputDir

if( !exists("SetupDIR") )
  SetupDIR = ANNOTDIR
# catn("setupDir is: ", SetupDIR)

#initialize variables
{
  TSS.windows = readRDS(paste0(SetupDIR,"Windows.rds"))
  genome.seqinfo = seqinfo(TSS.windows)
  ChrList = paste0("chr", c(1:22,"X", "Y"))
  if (TargetMod == "H3K4me3-scer") {
#    sacCer3.seqinfo = Seqinfo(genome="sacCer3")
    ChrList = (seqnames(genome.seqinfo))
  }

  Win.Sig.filename = paste0(SetupDIR, "WIN-sig.csv")
  Win.Sig = NULL
  if( !is.null(opt$options$windowsignatures))
    Win.Sig.filename = opt$options$windowsignatures
  if( file.exists( Win.Sig.filename) ) {
#    Win.Sig=readRDS(Win.Sig.filename)
    catn("Reading window signatures from", Win.Sig.filename)
    ZZ = read.csv(Win.Sig.filename, as.is = TRUE)
    Win.Sig = split(ZZ$window, ZZ$signature)
  }
  
  source(paste0(SourceDIR,"Background.R"))
  source(paste0(SourceDIR,"ComputeGeneCounts.R"))
  source(paste0(SourceDIR,"CommonGenes.R"))
  Genes.notexcluded = Genes[!Genes.excluded]
  source(paste0(SourceDIR, "cfChIP-Functions.R"))  
  source(paste0(SourceDIR, "YieldEstimation-Functions.R"))
  
  
  GenePrograms.filename =  paste0(SetupDIR, "Gene-sig.csv")
  Gene.Programs = NULL
  if( !is.null(opt$options$geneprograms))
    GenePrograms.filename = opt$options$geneprograms
  if( file.exists(GenePrograms.filename) ) {
    catn("Reading gene programs from", GenePrograms.filename)
    ZZ = read.csv(GenePrograms.filename, as.is = TRUE)
    ZZ = ZZ[ZZ$gene %in% Genes,]
    Gene.Programs = split(ZZ$gene,ZZ$program)
  }

  MetaGene.filename = paste0(SetupDIR, "Meta-genes.bed")
  MetaGene = NULL
  if( !is.null(opt$options$metagenes))  
    MetaGene.filename = opt$options$metagenes
  if( file.exists( MetaGene.filename)){
    catn("Reading meta gene info from", MetaGene.filename )
    MetaGene = import(MetaGene.filename)
  }
  
  MetaEnhancer.filename = paste0(SetupDIR, "Meta-enhancers.bed")
  MetaEnhancer = NULL
  if( !is.null(opt$options$metaenhancer))  {
    MetaEnhancer.filename = opt$options$metaenhancer
  }
  if( file.exists( MetaEnhancer.filename)) {
    catn("Reading meta enhancer info from", MetaEnhancer.filename )
    MetaEnhancer = import(MetaEnhancer.filename)
  }
  
  GeneDescription.filename = paste0(SetupDIR,"GeneDescription.csv")
  GeneDescription = NULL
  if( file.exists(GeneDescription.filename) ) {
    GeneDescription = read.csv(GeneDescription.filename, as.is = TRUE)
    rownames(GeneDescription) = GeneDescription[,1]
    GeneDescription=GeneDescription[,-1]
    GeneDescription$Description = gsub(",",";",GeneDescription$Description)
  }
  
  QC.filename = paste0(SetupDIR,"QC.bed")
  QC.bed = TSS.windows
  if( file.exists(QC.filename) ) {
    QC.bed = import(QC.filename)
    QC.bed$type = QC.bed$name
  }
  
  
}

params = cfChIP.Params()
params$DataDir = DataDir

RetainInMemory = FALSE

if( !is.null(opt$options$config) ) {
  config.filename = opt$options$config
} else
  config.filename = paste0(SetupDIR, "config.csv")

if( file.exists(config.filename) ) {
  catn("Reading configuration options from ",config.filename)
  config.data = read.csv(config.filename,as.is = TRUE, header = FALSE)
  for( i in 1:nrow(config.data) ) {
    params[[config.data[i,1]]] = config.data[i,2]
  }
}

params$Background = opt$options$background
params$GeneCounts = opt$options$count
params$Normalize = opt$options$normalize
params$OverExpressedGenes = opt$options$overexpressed
doOverExpressedGenes = opt$options$overexpressed
if( opt$options$all ) {
  params$Background = TRUE
  params$GeneBackground = TRUE
  params$GeneCounts = TRUE
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
}

doForce = opt$options$force

doMeta = opt$options$meta
doTracks = opt$options$tracks
doEnrichR = opt$options$enrichR
doSignatures = opt$options$signatures
doPrograms = opt$options$programs
doBackgroundPlot = opt$options$backgroundplot

if( DevelopmentMode ) {
#  params$Background = TRUE
#  params$GeneCounts = TRUE
#  params$Normalize = TRUE
#  params$OverExpressedGenes = TRUE
#  doMeta = TRUE
#  doOverExpressedGenes = TRUE
#  doForce = TRUE
  RetainInMemory = TRUE
  
}

if( doMeta || doTracks || doSignatures )
  params$Normalize = TRUE

if( doEnrichR  || doPrograms )
  params$OverExpressedGenes = TRUE

if( doTracks && !dir.exists(TracksDIR))
  dir.create(TracksDIR)

if( doOverExpressedGenes || doMeta )
  if( !dir.exists(OutputDir)) 
    dir.create(OutputDir)

if( doMeta && !dir.exists(MetaDIR))
  dir.create(MetaDIR)

doConsensus = !is.null(opt$options$consensus)
if(doConsensus ) {
  outputConsensus = outputFile(opt$options$consensus)
  params$GeneBackground = TRUE
  params$GeneCounts = TRUE
  RetainInMemory = TRUE
}

doGeneCountsTable = !is.null(opt$options$genecounts)
if(doGeneCountsTable ) {
  outputGeneCountsTable = outputFile(opt$options$genecounts)
  params$GeneCounts = TRUE
  RetainInMemory = TRUE
}

doGeneBackgroundTable = !is.null(opt$options$genebackground)
if(doGeneBackgroundTable ) {
  outputGeneBackgrouodTable = outputFile(opt$options$genebackground)
  params$GeneBackground = TRUE
  RetainInMemory = TRUE
}

doWinCountsTable = !is.null(opt$options$wincounts)
if(doWinCountsTable ) {
  outputWinCountsTable = outputFile(opt$options$wincounts)
  params$WinCounts = TRUE
  RetainInMemory = TRUE
}

doNormCountsTable = !is.null(opt$options$normcounts)
if( doNormCountsTable ) {
  outputNormCountsTable = outputFile(opt$options$normcounts)
  params$Normalize = TRUE
  RetainInMemory = TRUE
}

doPlotSignatures = !is.null(opt$options$plotsignatures)
if(doPlotSignatures  ) {
  outputPlotSignatures = outputFile(opt$options$plotsignatures)
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
  RetainInMemory = TRUE
}

doPlotPrograms = !is.null(opt$options$plotprograms)
if( doPlotPrograms ) {
  outputPlotPrograms = outputFile(opt$options$plotprograms)
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
  RetainInMemory = TRUE
}

doPlotEnrichments = !is.null(opt$options$plotenrichments)
if( doPlotEnrichments ) {
  outputPlotEnrichments = outputFile(opt$options$plotenrichments)
  params$Normalize = TRUE
  params$OverExpressedGenes = TRUE
  RetainInMemory = TRUE
}

doExport = !is.null(opt$options$export)
if(doExport ) {
  outputExport = outputFile(opt$options$export)
  RetainInMemory = TRUE
}

doQC = !is.null(opt$options$QC)
if( doQC ) {
  outputQC = outputFile(opt$options$QC)
  QC = list()
  params$Background = TRUE
}

doCluster = !is.null(opt$options$cluster)
if( doCluster ) {
  outputCluster = outputFile(opt$options$cluster)
  Cluster = list()
  params$Background = TRUE
  params$Normalize = TRUE
}

if( doMeta && is.null(params$MetaGene) )
  catn("Missing meta gene information")


BaseFileName <- function( fname, extList = c(".bed$", ".rdata$", ".bw$", ".tagAlign$", "-H3K4me3") ) {
#  x = file_path_sans_ext(fname)
  x = fname
  for( ext in extList )
    x = sub(ext, "", x)
  
  y = strsplit(x,"/")[[1]]
  n = length(y)
  z = y[n]
  return(z)
}

ProcessBEDFile = function( BFile ) {
  catn(BFile)
  if( !file.exists(BFile) && !file.exists(paste0(BFile,".bed")) && !file.exists(paste0(BFile,".tagAlign")) && !file.exists(paste0(BFile,".bw")) )
    BFile = paste0(BedDir, BFile)
  
  # catn(BFile)
  
  dat = cfChIP.ProcessFile(filename = BFile, param = params, Force = doForce )
 
  if( doQC )
    QC[[dat$Name]] <<- cfChIP.countQC(dat,GR = QC.bed, param= params)
  
  if( doCluster )
    Cluster[[dat$Name]] <<- dat$GeneCounts.QQnorm
  
  if(doBackgroundPlot) {
    if( params$Verbose ) catn(dat$Name, ": Ploting background estimate")
    cfChIP.BackgroundPlot(dat,Dir = OutputDir, Force = doForce)
  }
  
  if(doOverExpressedGenes) {
    if( params$Verbose ) catn(dat$Name, ": Printing overexpressed genes")
    cfChIP.OverExpressedGenes( dat, Dir = OutputDir, param = params, Force = doForce)
  }
  
  if( doMeta && !is.null(params$MetaGene) ) {
    if( params$Verbose ) catn(dat$Name, ": Ploting meta coverage")
    cfChIP.MetaPlot(dat, Dir = MetaDIR, param = params, Force = doForce)
  }
  
  if( doTracks ) {
    if( params$Verbose ) catn(dat$Name, ": saving normalized tracks")
    cfChIP.WriteTrack(dat, Dir = TracksDIR,  Force = doForce)
  }
  
  if( doSignatures ) {
    if( params$Verbose ) catn(dat$Name, ": Evaluating signatures")
    cfChIP.EvaluateSignatures( dat, Dir = OutputDir, param = params, Write = TRUE,  Force = doForce )
  }
  
  if( doPrograms ) {
    if( params$Verbose ) catn(dat$Name, ": Evaluating gene expression programs")
    cfChIP.EvaluatePrograms( dat, Dir = OutputDir, param = params, Write = TRUE,  Force = doForce )
    cfChIP.EvaluateProgramsHypG( dat, Dir = OutputDir, param = params, Write = TRUE,  Force = doForce )
  }
  
  if( doEnrichR ) {
    if( params$Verbose ) catn(dat$Name, ": Checking EnrichR")
    cfChIP.WriteEnrichR( dat, Dir = OutputDir,  Force = doForce )
  }
  
  
 
  if( RetainInMemory ) {
    # get read of the coverage since it takes too much memory....
    if(!DevelopmentMode) {
      dat$Cov = NULL
      dat$BED = NULL
      dat$BW = NULL
    }
    return(dat) 
  } else
    rm(dat)
}


ProcessBEDFileList = function( BFlist ) {
  L = lapply(BFlist, ProcessBEDFile)
  names(L) = sapply(L, function(l) l$Name)

  L  
}


LL = ProcessBEDFileList(Files)

if( doQC ) {
  catn("Writing QC numbers to ", outputQC)
  
  df = data.frame(do.call(rbind, QC))
  
  #  catn((df)
  write.csv(df, file = outputQC, quote = FALSE)
}

if( doCluster ) {
  catn("Saving file for clustering", outputCluster)
  X = do.call(cbind, Cluster)
  colnames(X) = names(Cluster)
  rownames(X) = Genes
  M = matrix(nr = nrow(X)+1, nc = ncol(X)+2)
  M[1,] = c("UID", "Name", colnames(X))
  M[2:(nrow(X)+1),1] = rownames(X)
  M[2:(nrow(X)+1),2] = rownames(X)
  M[2:(nrow(X)+1),3:(ncol(X)+2)] = formatC(log2(1+X))
  M[is.na(M)] = ""
  write.table(M, file = outputCluster, sep = "\t", 
              row.names = FALSE, col.names = FALSE,quote = FALSE )
}

if( RetainInMemory ) {
  SampleNames = sapply(LL, function(l) l$Name)
  names(LL) = SampleNames
  SampleOrder = names(LL)
  
  if( doExport ) {
    catn("Exporting data in RDS format")
    dat = LL[[1]]
    dat$Cov = NULL
    if(!is.null(dat$Counts))
      dat$Counts = sapply(LL, function(l) l$Counts)
    if(!is.null(dat$Background))
      dat$Background = lapply(LL, function(l) l$Background)
    if(!is.null(dat$GeneCounts))
      dat$GeneCounts = sapply(LL, function(l) l$GeneCounts )
    if(!is.null(dat$QQNorm))
      dat$QQNorm = sapply(LL, function(l) l$QQNorm )
    if(!is.null(dat$OverExpressedGenes))
      dat$OverExpressedGenes = lapply(LL, function(l) l$OverExpressedGenes)
    if(!is.null(dat$GeneBackground))  
      dat$GeneBackground = sapply(LL, function(l) l$GeneBackground)
    if(!is.null(dat$GenesCounts.QQnorm)) {
      # dat$GeneCounts.QQnorm = sapply(LL, function(l) l$GeneCounts.QQnorm )
      dat$GeneCounts.QQnorm = sapply(LL, function(l) l$GeneCounts.QQnorm[,1] )
    }
    if(!is.null(dat$Counts.QQnorm)) 
      dat$Counts.QQnorm = sapply(LL, function(l) l$Counts.QQnorm )
    
    saveRDS(dat, outputExport)
  }
  
  if( doGeneCountsTable ) {
    if( params$Verbose ) catn("Writing GeneCounts table to ", outputGeneCountsTable)
    GeneCounts = do.call("cbind", lapply(LL, function(l) l$GeneCounts))
    colnames(GeneCounts) = sapply(LL, function(l) l$Name)
    write.csv(GeneCounts,file = outputGeneCountsTable, quote = FALSE)
  }
  
  if( doGeneBackgroundTable ) {
    if( params$Verbose ) catn("Writing GeneBackground table to ", outputGeneBackgrouodTable)
    GeneBackground = do.call("cbind", lapply(LL, function(l) l$GeneBackground))
    colnames(GeneBackground) = sapply(LL, function(l) l$Name)
    write.csv(GeneBackground,file = outputGeneBackgroundTable, quote = FALSE)
  }
  
  if( doWinCountsTable ) {
    if( params$Verbose ) catn("Writing WinCounts table to ", outputWinCountsTable)
    WinCounts = do.call("cbind", lapply(LL, function(l) l$Counts))
    colnames(WinCounts) = sapply(LL, function(l) l$Name)
    write.csv(WinCounts,file = outputWinCountsTable, quote = FALSE)
  }
  
  if( doNormCountsTable ) {
    if( params$Verbose ) catn("Writing normalized GeneCounts table to ", outputNormCountsTable)
    GeneCounts = do.call("cbind", lapply(LL, function(l) l$GeneCounts.QQnorm))
    colnames(GeneCounts) = sapply(LL, function(l) l$Name)
    write.csv(GeneCounts,file = outputNormCountsTable, quote = FALSE)
  }
  
  if( doConsensus ) {
    if( params$Verbose ) catn("Writing normalized consensus to ", outputConsensus)
    GeneCounts = do.call("cbind", lapply(LL, function(l) l$GeneCounts))
    GeneBackground = do.call("cbind", lapply(LL, function(l) l$GeneBackground))
    colnames(GeneCounts) = sapply(LL, function(l) l$Name)
    colnames(GeneBackground) = colnames(GeneCounts)
    GeneDiff = GeneCounts - GeneBackground
    GeneDiff[GeneDiff < 0] = 0
    
    QQnorm = QQNormalizeGenes(GeneDiff, CommonG = CommonGenes )
    for( i in 1:length(QQnorm) )
      GeneDiff[,i] = GeneDiff[,i] * QQnorm[i]
    
    Avg = apply(GeneDiff,1,median)
    Avg = Avg * 10**6 / sum(Avg)
    
    saveRDS(Avg, outputConsensus)
  }
  
  if( doPlotSignatures ) 
    cfChIP.plotSignatures(LL, outputPlotSignatures, params)
  
  if( doPlotPrograms ) 
    cfChIP.plotPrograms(LL, outputPlotPrograms, params)

  if( doPlotEnrichments ) 
    cfChIP.plotEnrichments(LL, outputPlotEnrichments, params)
  
}
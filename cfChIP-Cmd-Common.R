if(!exists("LoadMeta") )
  LoadMeta = TRUE

if(!exists("LoadSig") )
  LoadSig = TRUE

if(!exists("LoadQC") )
  LoadQC = TRUE

if(!exists("LoadDescription") )
  LoadDescription = TRUE

catn = function(...) { cat(...,"\n") }

initial.options <- commandArgs(trailingOnly = FALSE)
if( !DevelopmentMode ) {
  DataDir = paste0(getwd(), "/")
  ANNOTDIR = DataDir
} else if(0) {
  Mod = "H3K4me3/"
  ANNOTDIR = paste0(SourceDIR, "SetupFiles/", Mod)
  DataDir = paste0("~/GoogleDrive/BloodChIP/Analysis-Paper/Samples/", Mod)
  BEDDIR = paste0("~/GoogleDrive/BloodChIP/Analysis-Paper/BED/", Mod)
  OutputDir = paste0("~/GoogleDrive/BloodChIP/Analysis-Paper/Output/", Mod)
  Files = c("H001.1.K4me3", "H001.2.K4me3")
  SetupDIR  = ANNOTDIR
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

standard_option_list = list(  
# directory 
  make_option(c("-d", "--datadir"), type = "character", default = NULL,
              help="data directory"),
  make_option(c("-b", "--BEDdir"), type="character", default=NULL,
              help="location of BED files" ),
  make_option(c("-a", "--annotationdir"), type="character", default=NULL,
              help="location of annotation files" ),
  make_option(c("-o", "--outputdir"), type="character", default=NULL,
              help="location of output files" ),
  make_option(c("-m", "--mod"), type="character", default="H3K4me3",
              help="name of modification, used in combination with -r" ),
  make_option(c("-r", "--root"), type="character", default=NULL,
              help="name of root directory, input/output is assumed to be in Root/XXX/mod" )
  
)


opt_parser = OptionParser(option_list=c(standard_option_list,option_list));
if( exists("trailing.options") ) {
  opt = parse_args(opt_parser, args = trailing.options, positional_arguments = TRUE)
} else
  opt = parse_args(opt_parser, positional_arguments = TRUE)


Files = opt$args

if(!DevelopmentMode && length(Files) == 0 && HelpOnEmptyArguments )
  parse_args(opt_parser,args = "-h", print_help_and_exit = TRUE)

catn("Initializing")
suppressMessages(library(rtracklayer))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))

TargetMod = opt$options$mod
if( !is.null(opt$options$root)) {
  RootDir =  extendDir(opt$options$root)
  if( is.null(opt$options$datadir))
    DataDir = paste0(RootDir, "Samples/", TargetMod, "/")
  if( is.null(opt$options$BEDdir))
    BEDDIR= paste0(RootDir, "BED/", TargetMod, "/")
  if( is.null(opt$options$outputdir) )
    OutputDir = paste0(RootDir, "Output/", TargetMod, "/")
} 
ANNOTDIR = paste0(SourceDIR, "SetupFiles/", TargetMod, "/")

if( !is.null(opt$options$datadir) )
  DataDir = extendDir(opt$options$datadir)

DIR = DataDir

if( !is.null(opt$options$annotationdir) )
  ANNOTDIR = extendDir(opt$options$annotationdir)

if( !DevelopmentMode && !exists("BEDDIR")) {
  BEDDIR = DataDir
}

if( !is.null(opt$options$BEDdir) )
  BEDDIR = extendDir(opt$options$BEDdir)

if( !is.null(opt$options$outputdir)) {
  OutputDir = extendDir(opt$options$outputdir)
} else
  if( !exists("OutputDir") ) {
    if( !DevelopmentMode  ) {
      OutputDir = OrigDir
    } else
      OutputDir = paste0("~/Data/BloodChIP/Output/", Mod)
  }


if( !exists("SetupDIR") )
  SetupDIR = ANNOTDIR

catn("OutputDir", OutputDir)
catn("BEDDIR", BEDDIR)
catn("SetupDIR", SetupDIR)
catn("DataDir", DataDir)

#initialize variables
{
  TSS.windows = readRDS(paste0(SetupDIR,"Windows.rds"))
  genome.seqinfo = seqinfo(TSS.windows)
  ChrList = paste0("chr", c(1:22,"X", "Y"))
  
  Win.Sig.filename = paste0(SetupDIR, "WIN-sig.csv")
  Win.Sig = NULL
  if( LoadSig && !is.null(opt$options$windowsignatures))
    Win.Sig.filename = opt$options$windowsignatures
  if( LoadSig && file.exists( Win.Sig.filename) ) {
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
  
  
  GenePrograms.filename =  paste0(SetupDIR, "Gene-sig.csv")
  Gene.Programs = NULL
  if( LoadSig && !is.null(opt$options$geneprograms))
    GenePrograms.filename = opt$options$geneprograms
  if( LoadSig && file.exists(GenePrograms.filename) ) {
    catn("Reading gene programs from", GenePrograms.filename)
    ZZ = read.csv(GenePrograms.filename, as.is = TRUE)
    ZZ = ZZ[ZZ$gene %in% Genes,]
    Gene.Programs = split(ZZ$gene,ZZ$program)
  }
  
  MetaGene.filename = paste0(SetupDIR, "Meta-genes.bed")
  MetaGene = NULL
  if( !is.null(opt$options$metagenes))  
    MetaGene.filename = opt$options$metagenes
  if(  LoadMeta && file.exists( MetaGene.filename)){
    catn("Reading meta gene info", MetaGene.filename )
    MetaGene = import(MetaGene.filename)
  }
  
  MetaEnhancer.filename = paste0(SetupDIR, "Meta-enhancers.bed")
  MetaEnhancer = NULL
  if(  !is.null(opt$options$metaenhancer))  
    MetaEnhancer.filename = opt$options$metaenhancer
  
  if( LoadMeta && file.exists( MetaEnhancer.filename)) {
    catn("Reading meta enhancer info", MetaEnhancer.filename )
    MetaEnhancer = import(MetaEnhancer.filename)
  }
  
  GeneDescription.filename = paste0(SetupDIR,"GeneDescription.csv")
  GeneDescription = NULL
  if( LoadDescription && file.exists(GeneDescription.filename) ) {
    GeneDescription = read.csv(GeneDescription.filename, as.is = TRUE)
    rownames(GeneDescription) = GeneDescription[,1]
    GeneDescription=GeneDescription[,-1]
    GeneDescription$Description = gsub(",",";",GeneDescription$Description)
  }
  
  QC.filename = paste0(SetupDIR,"QC.bed")
  QC.bed = TSS.windows
  if( LoadQC && file.exists(QC.filename) ) {
    QC.bed = import(QC.filename)
    QC.bed$type = QC.bed$name
  }
  
  params = cfChIP.Params()
  params$DataDir = DataDir
}


BaseFileName <- function( fname, extList = c(".bed$", ".rdata$", ".bw$", ".tagAlign$", "-H3K4me3", "-H3K4me2", "-H3K4me1", "-H3K36me3") ) {
  #  x = file_path_sans_ext(fname)
  x = fname
  for( ext in extList )
    x = sub(ext, "", x)
  
  y = strsplit(x,"/")[[1]]
  n = length(y)
  z = y[n]
  return(z)
}

BuildGlobalMatrices = function(LL) {
  dat = LL[[1]]
  Counts = NULL
  Background = NULL
  GeneCounts = NULL
  QQNorm = NULL
  OverExpressedGenes = NULL
  GeneCounts.QQnorm = NULL
  Counts.QQnorm = NULL
  
  if(!is.null(dat$Counts))
    Counts <<- sapply(LL, function(l) l$Counts)
  if(!is.null(dat$Background))
    Background <<- lapply(LL, function(l) l$Background)
  if(!is.null(dat$GeneCounts))
    GeneCounts <<- sapply(LL, function(l) l$GeneCounts )
  if(!is.null(dat$QQNorm))
    QQNorm <<- sapply(LL, function(l) l$QQNorm )
  if(!is.null(dat$OverExpressedGenes))
    OverExpressedGenes <<- lapply(LL, function(l) l$OverExpressedGenes)
  if(!is.null(dat$GeneBackground))
    GeneBackground <<- sapply(LL, function(l) l$GeneBackground)
  if(!is.null(dat$GeneCounts.QQnorm))
    GeneCounts.QQnorm <<- sapply(LL, function(l) l$GeneCounts.QQnorm )
  if(!is.null(dat$Counts.QQnorm)) 
    Counts.QQnorm <<- sapply(LL, function(l) l$Counts.QQnorm )
}

LoadListofSamples = function(Files, f = cfChIP.ProcessFile ) {
  LL = lapply(Files, f)
  LL = LL[!sapply(LL, is.null)]
  names(LL) = sapply(LL, function(l) l$Name)
  BuildGlobalMatrices(LL)
  return(LL)
}

LoadAllSamples = function(f = cfChIP.ProcessFile ) {
  Fs = list.files(DataDir, pattern = ".*.rdata")
  LL = LoadListofSamples(Fs, f)
  BuildGlobalMatrices(LL)
  return(LL)
}


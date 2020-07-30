require(optparse)

if(!exists("LoadMeta") )
  LoadMeta = TRUE

if(!exists("LoadSig") )
  LoadSig = TRUE

if(!exists("LoadQC") )
  LoadQC = TRUE

if(!exists("LoadDescription") )
  LoadDescription = TRUE

if(!exists("ExpandFileLists"))
  ExpandFileLists = TRUE

catn = function(...) { cat(...,"\n") }


initial.options <- commandArgs(trailingOnly = FALSE)
if( exists("DevelpomentMode") && !DevelopmentMode ) {
  DataDir = paste0(getwd(), "/")
  ANNOTDIR = DataDir
} else  {
  DevelopmentMode = TRUE
#  Mod = "H3K4me3/"
#  ANNOTDIR = paste0(SourceDIR, "SetupFiles/", Mod)
#  DataDir = paste0("~/BloodChIP/Data/Analysis-Paper/Samples/", Mod)
#  BEDDIR = paste0("~/BloodChIP/Data/Analysis-Paper/BED/", Mod)
#  OutputDir = paste0("~/BloodChIP/Data/Analysis-Paper/Output/", Mod)
#  SetupDIR  = ANNOTDIR
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
  make_option(c("-r", "--root"), type="character", default="~/BloodChIP/Data/raw-data",
              help="name of root directory, input is assumed to be in Root/XXX/mod" ),
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="name of project root directory, output is assumed to be in Project/XXX/mod" ),
  make_option("--config", type="character", default = NULL,
              help = "Specify an alternative configuration file")
)

if(!exists("option_list"))
  option_list = list()

opt_parser = OptionParser(option_list=c(standard_option_list,option_list));
if( exists("trailing.options") ) {
  catn("CMD line = ", trailing.options)
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
  DataDir = paste0(RootDir, "Samples/", TargetMod, "/")
  BedDir = paste0(RootDir, "BED/", TargetMod, "/")
  BEDDIR = paste0(RootDir, "BED/", TargetMod, "/")
  #  TracksDIR = paste0(RootDir, "Tracks/", TargetMod, "/")
  #  OutputDir = paste0(RootDir, "Output/", TargetMod, "/")
} 
if( !is.null(opt$options$project)) {
  ProjectDir =  extendDir(opt$options$project)
  TracksDIR = paste0(ProjectDir, "Tracks/", TargetMod, "/")
  OutputDir = paste0(ProjectDir, "Output/", TargetMod, "/")
} 
ANNOTDIR = paste0(SourceDIR, "SetupFiles/", TargetMod, "/")

if( !is.null(opt$options$datadir) )
  DataDir = extendDir(opt$options$datadir)

DIR = DataDir

if( !is.null(opt$options$annotationdir) )
  ANNOTDIR = extendDir(opt$options$annotationdir)

if(!exists("BEDDIR")) {
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
  
  Win.Sig.filename = paste0(SetupDIR, "Win-sig.csv")
  Win.Sig = NULL
  Win.Sig.Ref = NULL
  if( !is.null(opt$options$windowsignatures))
    Win.Sig.filename = opt$options$windowsignatures
  if( file.exists( Win.Sig.filename) ) {
    #    Win.Sig=readRDS(Win.Sig.filename)
    catn("Reading window signatures from", Win.Sig.filename)
    ZZ = read.csv(Win.Sig.filename, as.is = TRUE)
    Win.Sig = split(ZZ$window, ZZ$signature)
    
    Win.Sig.Ref.filename = sub("csv$", "rds", Win.Sig.filename)
    if( file.exists(Win.Sig.Ref.filename)) {
      catn("Reading window signatures reference from", Win.Sig.Ref.filename)
      Win.Sig.Ref = readRDS(Win.Sig.Ref.filename)
    }
  } else
    if( !is.null(opt$options$windowsignatures) )
      catn("Cannot find file", opt$options$windowsignatures)
  
  source(paste0(SourceDIR,"Background.R"))
  source(paste0(SourceDIR,"ComputeGeneCounts.R"))
  source(paste0(SourceDIR,"CommonGenes.R"))
  Genes.notexcluded = Genes[!Genes.excluded]
  source(paste0(SourceDIR, "cfChIP-Functions.R"))  
  
  
  GenePrograms.filename =  paste0(SetupDIR, "Gene-sig.csv")
  Gene.Programs = NULL
  Gene.Programs.Ref = NULL
  Gene.Programs.Partition = NULL
  if( !is.null(opt$options$geneprograms))
    GenePrograms.filename = opt$options$geneprograms
  if( file.exists(GenePrograms.filename) ) {
    catn("Reading gene programs from", GenePrograms.filename)
    ZZ = read.csv(GenePrograms.filename, as.is = TRUE)
    ZZ = ZZ[ZZ$gene %in% Genes,]
    Gene.Programs = split(ZZ$gene,ZZ$program)
    
    GenePrograms.Ref.filename = sub("csv$", "rds", GenePrograms.filename)
    if( file.exists(GenePrograms.Ref.filename)) {
      catn("Reading gene programs reference from", GenePrograms.Ref.filename)
      Gene.Programs.Ref = readRDS(GenePrograms.Ref.filename)
    }
  } else
    if(!is.null(opt$options$geneprograms)) {
      catn("Cannot find file", opt$options$geneprograms)
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


BaseFileName <- function( fname, extList = c(".bed$", ".rdata$", ".bw$", ".tagAlign$" #, "-H3K4me3", "-H3K4me2", "-H3K4me1", "-H3K36me3"
                                             ) ) {
  #  x = file_path_sans_ext(fname)
  x = fname
  x = sub(".gz$", "",x)
  for( ext in extList )
    x = sub(ext, "", x)
  
  y = strsplit(x,"/")[[1]]
  n = length(y)
  z = y[n]
  return(z)
}


expandFiles = function(f) {
  if( grepl(".txt$", f) | grepl(".csv$", f) ) {
    if(!file.exists(f)) {
      catn("Missing list file", f)
      return(NULL)
    }
    catn("Reading sample list from",f)
    File.list = read.table(f, as.is = TRUE)
    File.list = File.list[,1]
    names(File.list) = NULL
    as.list(File.list)
  } else
    f
}

if(ExpandFileLists)
  Files = unique(unlist(lapply(Files,expandFiles)))

BuildGlobalMatrices = function(LL) {
  dat = LL[[1]]
  Counts = NULL
  Background = NULL
  GeneCounts = NULL
  QQNorm = NULL
  OverExpressedGenes = NULL
  GeneCounts.QQnorm = NULL
  Counts.QQnorm = NULL
  Heights = NULL
  GeneHeights = NULL
  GeneCNV = NULL
  if(!is.null(dat$Counts))
    Counts <<- sapply(LL, function(l) if(is.null(l)) {rep(NA,length(dat$Counts))} else l$Counts)
  if(!is.null(dat$Background))
    Background <<- lapply(LL, function(l) l$Background)
  if(!is.null(dat$GeneCounts))
    GeneCounts <<- sapply(LL, function(l) if(is.null(l)) {rep(NA,length(dat$GeneCounts))} else l$GeneCounts )
  if(!is.null(dat$QQNorm)) 
    QQNorm <<- sapply(LL, function(l) if(is.null(l)) {NA} else as.numeric(l$QQNorm) )
  if(!is.null(dat$OverExpressedGenes))
    OverExpressedGenes <<- lapply(LL, function(l) l$OverExpressedGenes)
  if(!is.null(dat$GeneBackground))
    GeneBackground <<- sapply(LL, function(l) if(is.null(l)) {rep(NA,length(dat$GeneBackground))} else l$GeneBackground)
  if(!is.null(dat$GeneCounts.QQnorm))
    GeneCounts.QQnorm <<- sapply(LL, function(l) if(is.null(l)) {rep(NA,length(dat$GeneCounts.QQnorm))} else l$GeneCounts.QQnorm )
  if(!is.null(dat$Counts.QQnorm)) 
    Counts.QQnorm <<- sapply(LL, function(l) if(is.null(l)) {rep(NA,length(dat$Counts.QQnorm))} else l$Counts.QQnorm )
  if(!is.null(dat$Heights)) 
    Heights <<- sapply(LL, function(l) if(is.null(l)) {rep(NA,length(dat$Heights))} else l$Heights)
  if(!is.null(dat$GeneHeights)) 
    GeneHeights <<- sapply(LL, function(l) if(is.null(l)) {rep(NA,length(dat$GeneHeights))} else l$GeneHeights)
  if(!is.null(dat$GeneCNV)) 
    GeneCNV <<- sapply(LL, function(l) if(is.null(l)) {rep(NA,length(dat$CNV))} else l$GeneCNV)
  
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


ReadSignatureList = function(Sig.filename, type = "signature") {
  Sig.filename = unlist(strsplit(Sig.filename, ","))
  
  for( f in Sig.filename )
    if( !file.exists(f) ) 
      catn("Cannot find", type,"file", f)
  
  Sig.list = list()
  
  for( f in Sig.filename ) {
    catn("Reading", type, "from", f)
    ZZ = read.csv(f, as.is = TRUE)
    Sig = split(ZZ[,2], ZZ[,1])
    if( type == "programs")
      Sig = sapply(Sig, function(s) intersect(s, Genes.notexcluded))
    f.Ref = sub("csv$", "rds", f)
    if( file.exists(f.Ref) ) {
      catn("Reading", type, "reference from", f.Ref)
      Sig.Ref = readRDS(f.Ref)
      if( length(Sig) != length(Sig.Ref$avg) || any(names(Sig) != names(Sig.Ref$avg)) ) {
        catn("Mismatch in consensus of signature", f)
        A = intersect(names(Sig), names(Sig.Ref$avg))
        Sig = Sig[A]
        Sig.Ref$avg = Sig.Ref$avg[A]
        Sig.Ref$var = Sig.Ref$var[A]
      }
    } else{
      x = rep(NA, length(Sig))
      names(x) = names(Sig)
      Sig.Ref = list( avg = x, var = x)
    }
    Sig.list[[f]] = list(Sig = Sig, avg = Sig.Ref$avg, var = Sig.Ref$var )
  } 
  names(Sig.list) = NULL
  Sig = do.call("c", lapply(Sig.list, function(x) x$Sig))
  Sig.Ref = list( avg = do.call("c", lapply(Sig.list, function(x) x$avg)),
                  var = do.call("c", lapply(Sig.list, function(x) x$var))
  )
  list(Sig = Sig, Ref = Sig.Ref)
}


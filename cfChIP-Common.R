library(ggplot2)
library(rtracklayer)
library(preprocessCore)
library(Matrix)
library(reshape2)
source(paste0(SourceDIR,'DensityScatter.R'))

if( !exists( "TSS.windows") ) {
  TSS.windows = readRDS(paste0(ANNOTDIR,"Windows.rds"))
#  Enh.windows = readRDS(paste0(ANNOTDIR,"Enh-windows.rds"))
  genome.seqinfo = seqinfo(TSS.windows)
  ChrList = paste0("chr", c(1:22,"X", "Y"))
}

source(paste0(SourceDIR, "LoadSamples.R"))
if( !exists("WINS") ) {
  LoadSamples(FindSamples())
  if( exists("WINS"))
    Samples = colnames(WINS)
}

source(paste0(SourceDIR,"Background.R"))
source(paste0(SourceDIR,"ComputeGeneCounts.R"))
source(paste0(SourceDIR,"CommonGenes.R"))

if( exists("GeneCounts")) {
  GeneDiff = GeneCounts - GeneBackground
  GeneDiff[GeneDiff < 0] = 0
  
  logGeneCounts.raw = log2(GeneCounts.QQnorm+1)
}

HealthyNormAvg = Healthy.GeneCount
HealthyNormAvg.log = log2(HealthyNormAvg)

uniqueGenes = function(G) {

  I = sapply(G, function(g) all(as.character(TSS.windows$name[GeneWindows[grep(g,TSS.windows$name[GeneWindows])]]) == g))
  
  G1 = G[I]
  G = G[!I]
  Ws = do.call(c,lapply(G,function(g) { v= as.vector(grep(g,TSS.windows$name[GeneWindows])); names(v) = rep(g,length(v)); v}))
  G2 = unique(names(Ws)[match(unique(Ws), Ws)])
  c(G1,G2)
}


if( !exists("Atlas.Win.Sig")) {
  Atlas.Win.Sig.filename = paste0(SetupDIR, "WIN-sig.rds")
  if( file.exists(Atlas.Win.Sig.filename))
    Atlas.Win.Sig=readRDS(Atlas.Win.Sig.filename)
} 

if( exists("WINS") && !is.null(WINS)) {
  SampleGroupOrder = c("H", "M", "C", "T")
  SampleOrder = do.call("c",lapply(SampleGroupOrder, function(g) grep(paste0("^",g), colnames(WINS))))
  SampleOrder = c(SampleOrder, which( !(1:dim(WINS)[2] %in% SampleOrder )) )
}
source(paste0(SourceDIR,"cfChIP-util.R"))

GeneDescription.filename = paste0(SetupDIR,"GeneDescription.csv")
GeneDescription = NULL
if( file.exists(GeneDescription.filename) ) {
  GeneDescription = read.csv(GeneDescription.filename, as.is = TRUE)
  rownames(GeneDescription) = GeneDescription[,1]
  GeneDescription=GeneDescription[,-1]
  GeneDescription$Description = gsub(",",";",GeneDescription$Description)
}


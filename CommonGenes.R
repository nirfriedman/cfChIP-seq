

WinChr = as.character(chrom(TSS.windows))
v = rep(0, times=length(WinChr))
#v[WinChr == "chrX" | WinChr == "chrY"] = 1
vy = v
vy[WinChr == "chrY"] = 1
vx = v
vx[WinChr == "chrX"] = 1
vm = v
vm[WinChr == "chrM"] = 1
Genes.chrY =  as.vector(t( vy[GeneWindows] %*% Win2Gene.matrix ) > 0)
Genes.chrX =  as.vector(t( vx[GeneWindows] %*% Win2Gene.matrix ) > 0)
Genes.chrM =  as.vector(t( vm[GeneWindows] %*% Win2Gene.matrix ) > 0)
rm(vx)
rm(vy)
rm(vm)
Genes.excluded = Genes.chrY | Genes.chrM
# grepl("-AS$", Genes) | grepl("-AS1$", Genes) | grepl("-AS2$", Genes) |
  
ExclusionPrefix = read.csv(paste0(SetupDIR,"FilterGenes.csv"), as.is = TRUE)
for( p in ExclusionPrefix[,1]) {
  Genes.excluded = Genes.excluded | grepl(paste0("^",p),Genes)
}
Genes.notexcluded = Genes[!Genes.excluded]

if( exists("GeneCounts"))
  Samples.select = colnames(GeneCounts)

if( !exists("CommonGenes") ) {
  catn("Loading Common Gene List", paste0(SetupDIR, "CommonGenes.rds"))
  CommonGenes = readRDS(paste0(SetupDIR, "CommonGenes.rds"))
  CommonGenes = CommonGenes[CommonGenes%in%Genes]
}

if( !exists("Healthy.GeneCount")) {
  HealthyRef.filename = paste0(SetupDIR,"HealthyRef.rds")
  Healthy.GeneCount = NULL
  Healthy.WinCount = NULL
  Healthy.GeneCount.var = NULL
  Healthy.WinCount.var = NULL
  
  if( file.exists(HealthyRef.filename)) {
    catn("Loading Reference Count values from", HealthyRef.filename)
    L = readRDS(HealthyRef.filename)
    if( is.list(L) ) {
      Healthy.GeneCount = L$Gene.avg
      Healthy.GeneCount.var = L$Gene.var
      Healthy.WinCount = L$Win.avg
      Healthy.WinCount.var = L$Win.var
      names(Healthy.GeneCount) = Genes
      names(Healthy.GeneCount.var) = Genes
    } else {
      # backward compatability if the file is just the averages
      Healthy.GeneCount = L
    }
    rm(L)
  }
}


Win.notexcluded = (WinChr != "chrM") & (WinChr != "chrX" ) &
  (WinChr != "chrY") # & !WinUnmapable & !(1:dim(WINS)[1] %in% HighValue.exclude)
Win.notexcluded = Win.notexcluded & !(TSS.windows$name %in% Genes[Genes.excluded])

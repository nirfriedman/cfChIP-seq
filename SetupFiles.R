#!/usr/bin/env Rscript --vanilla
#


initial.options <- commandArgs(trailingOnly = FALSE)
if( !any(grepl("--interactive", initial.options)) ) {
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  SourceDIR <- paste0(dirname(script.name),"/")
  DataDir = paste0(getwd(), "/")
} else {
  SourceDIR = "~/Google Drive/Src/BloodChIPAnalysis/"
#  DataDir = "~/Data/BloodChIP/SetupFiles/H3K4me3/"
  DataDir = "~/GoogleDrive/Src/BloodChIPAnalysis/SetupFiles/H3K4me3-scer/"
  # SourceDIR = "~/GoogleDrive/BloodChIPAnalysis"
  # DataDir = "~/GoogleDrive/BloodChIPAnalysis/SetupFiles/H3K4me3-scer/"
  TargetMod = "H3K4me3-scer"
}

suppressWarnings(library(ggplot2))
suppressWarnings(library(preprocessCore))
suppressWarnings(library(Matrix))
suppressWarnings(library(reshape2))
suppressWarnings(library(rtracklayer))

ANNOTDIR = DataDir

TSS.windows.filename = paste0(ANNOTDIR,"Windows.rds")
TSS.windows = readRDS(TSS.windows.filename)
genome.seqinfo = seqinfo(TSS.windows)
ChrList = paste0("chr", c(1:22,"X", "Y"))
if (TargetMod == "H3K4me3-scer") {
  sacCer3.seqinfo = Seqinfo(genome="sacCer3")
  ChrList = (seqnames(sacCer3.seqinfo))
}

BackgroundModel.filename = paste0(DataDir,"BackgroundModel.rds")
print("Building Background model")
sortGR = function(GR) {
  GR.order = lapply(unique(chrom(GR)), function(chr) {
    W = which(as.vector(chrom(GR)) == chr)
    W[order(start(GR[W]))]
  })
  GR.order = do.call("c", GR.order)
  GR[GR.order]  
}

buildOverlapingTiles = function(width, jump = width/4) {
  T = tileGenome(hg19.seqinfo, tilewidth = width, cut.last.tile.in.chrom = TRUE )
  T = T[chrom(T) %in% ChrList ]
  T0 = T
  mWidth = min(width(T0))
  for( o in seq(jump, width-1, by=jump)) {
    T = suppressWarnings(c(T, trim(shift(T0,o))))
  }
  T = T[width(T) >= mWidth] 
  T = sortGR(T)
  T
}

Background.Regions.width =  10**7
Background.Regions = buildOverlapingTiles(Background.Regions.width)
Background.Regions.num = length(Background.Regions)
Background.Uniq.Regions = disjoin(Background.Regions, with.revmap=TRUE)

Background.Inds.width = 5*10**6
Background.Inds = buildOverlapingTiles(Background.Inds.width)
Background.Inds.num = length(Background.Inds)

Background.Uniq.Inds = disjoin(Background.Inds, with.revmap=TRUE)

Background.RegionChr = as.character(chrom(Background.Regions))
Background.IndRegion = nearest(Background.Inds,subject = resize(Background.Regions, width = 10, fix = "center"))

WinChr = as.character(chrom(TSS.windows))

Background.RegionWindows =  lapply(1:Background.Regions.num, function(i) subjectHits(findOverlaps(Background.Regions[i], TSS.windows)))
Background.IndWindows = lapply(1:Background.Inds.num, function(i) subjectHits(findOverlaps(Background.Inds[i], TSS.windows)))

Background.WindowInd = Rle(nearest(TSS.windows,Background.Uniq.Inds))
Background.WindowRegion = Rle(nearest(TSS.windows,Background.Uniq.Regions))

saveRDS(list(Background.Regions.width = Background.Regions.width,
             Background.Regions = Background.Regions,
             Background.Regions.num = Background.Regions.num,
             Background.Uniq.Regions = Background.Uniq.Regions,
             Background.Inds.width = Background.Inds.width,
             Background.Inds = Background.Inds,
             Background.Inds.num = Background.Inds.num,
             Background.Uniq.Inds = Background.Uniq.Inds,
             Background.RegionChr = Background.RegionChr,
             Background.IndRegion = Background.IndRegion,
             WinChr = WinChr,
             Background.RegionWindows = Background.RegionWindows,
             Background.IndWindows = Background.IndWindows,
             Background.WindowInd = Background.WindowInd,
             Background.WindowRegion = Background.WindowRegion),
        BackgroundModel.filename)


Win2Gene.Matrix.filename = paste0(DataDir,"Win2Gene.rds")
print("Computing Window to Gene mapping")

GeneWindows = which(TSS.windows$name != "." & !is.na(TSS.windows$name))
GeneLists = as.character(unique(TSS.windows[GeneWindows]$name))
GeneLists = GeneLists[GeneLists != "."]
Genes = unique(do.call("c",strsplit(GeneLists,";") ))
Genes = Genes[!is.na(Genes)]
GL = strsplit(as.character(TSS.windows[GeneWindows]$name),";")
names(GL) = GeneWindows

Ls = lapply(GL, length)
GenesId = 1:length(Genes)
names(GenesId) = Genes
max.L =  max(do.call("c", Ls))
rs = c()
cs = c()
for( i in 1:max.L ) {
  rs = c(rs, rep(which(Ls == i), each = i))
  cs = c(cs,  GenesId[do.call("c",GL[which(Ls == i)])])
}

Win2Gene.matrix  = sparseMatrix(i=as.integer(rs),j=cs,x = 1, dims = c(length(GeneWindows), length(Genes)), dimnames = list(as.character(GeneWindows),Genes))


MultiPromoterGenes = sapply(Genes,function(g) {
  I = which( Win2Gene.matrix[,g] > 0)
  if( length(I) > 0 )
    max(I) - min(I) > length(I)-1
  else
    NA
})

saveRDS(list(Matrix = Win2Gene.matrix, 
             Multi = MultiPromoterGenes), Win2Gene.Matrix.filename)


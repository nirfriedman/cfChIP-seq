FitNonZeroPoisson <- function( ObsStat ) {
  N = sum(ObsStat)
  M = max(as.integer(names(ObsStat)))
  E = sum(ObsStat * as.integer(names(ObsStat))) / N
  ll <- function( p ) {
    -(E * log(p) - p - log(1-exp(-p)))
  }
  tryCatch( {
    opt = optimize(f = ll, interval = c(0,M))
    return(opt$minimum)
  }, error = function(e) print(e) )
  return(1)
}

FitMixPoisson <- function( ObsStat ) {
  i = list( lambda = FitNonZeroPoisson(ObsStat)/2, p1=10, p2=0, p3=0, p4=0,p5=1,p6=0,p7=0,p8=0,p9=0,p10=0,p11=0,p12=0)
  
  ll <- function( p ) do.call(MixPoissonLL, c(list(ObsStat), p))
  
  estparam = i
  tryCatch( {
    est <- optim(i,ll,control=list(maxit=10000))
    #      print(est$value)
    estparam = est$par
    p = c(estparam["p1"], estparam["p2"], estparam["p3"], estparam["p4"], 
          estparam["p5"], estparam["p6"], estparam["p7"],estparam["p8"],
          estparam["p9"], estparam["p10"], estparam["p11"], estparam["p12"])
    p = exp(p)
    p = p/sum(p)
    lambda = estparam["lambda"]
    return(list( lambda = lambda, p = p))
    
  }, error = function(e) print(e) )
  #  print(lambda)
  #  plot(table(dups),type='h'); plot(1:31, MixPoisson(1:31,lambda, p), type = 'h')
  return(list(lambda = 1, p = c(1,0,0,0)))
}

MixPoisson <- function(x,lambda,p) {
  y = x*0
  for( i in 1:length(p) ){
    y = y+ p[i]*dpois(x,lambda*i)
  }
  return(y)
}

MixPoissonLL <- function(Obs, lambda, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12) {
  x = 1:length(Obs)
  p = c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
  z = sum(exp(p))
  q = MixPoisson( x, lambda, exp(p)/z )
  t = 1-MixPoisson( 0, lambda, exp(p)/z )
  #  print(Obs)
  qll = -sum( Obs * log(q/t) )
}

# estimate library complexity.  based on Nir's EstimateYield.R
CalculateSequencingYield <- function (sample_name = NULL, modification = NULL) {
  catn(sample_name, ": estimating library complexity")
  BFile = paste0(tagAlignDIR, modification, "/", sample_name, ".not_uniq.tagAlign")
  params.MaxFragLen = 1000
  params.MinFragLen = 50
  # params$reuseSavedData = FALSE # TODO change back to true or delete 
  # BFile = "Google Drive/BloodChIP/Analysis/tagAligns/H027_1_i_S19.not_uniq.tagAlign"
  # BFile = sub(pattern = ".K4me3.tagAlign", replacement = "", x = BFile)
  sample_name = basename(BFile)
  # dat = cfChIP.ProcessFile(filename = BFile, param = params)
  if (!file.exists(BFile)) {
    catn(sample_name, ": ERROR, file does not exist")
    return(NULL)
  }
  dat = import(BFile, format = "BED")
  dat = dat[width(dat) <= params.MaxFragLen & width(dat) > params.MinFragLen]
  aligned_reads = length(dat)
  dat.uniq = unique(dat)
  # count duplicate segments
  dups <- countOverlaps(dat.uniq, dat, type = "equal")
  dups = as.data.frame(table(dups))
  dups
  Hists = as.numeric(dups$dups)
  aligned_frags = length(dat)
  N.uniq = sum(dups$Freq)
  
  dupStat = dups$Freq
  names(dupStat) = dups$dups
  Nuniq = sum(dupStat)
  lambda1 = FitNonZeroPoisson(dupStat)
  Total1 = Nuniq /(1 - exp(-lambda1))
  Summary = data.frame(N.frags = c(), N.uniq = c(), Poisson.Lambda = c(), Poisson.total = c(), Poisson.yield = c(), Mix.Lambda = c(), 
                       Mix.Total = c(), Mix.yield = c(), Uniq.reads = c() )
  Summary[sample_name, "N.frags"] = aligned_frags
  Summary[sample_name, "N.uniq"] = N.uniq
  Summary[sample_name, "Poisson.Lambda"] = lambda1
  Summary[sample_name, "Poisson.total"] = Total1
  Summary[sample_name,"Poisson.yield"] = Nuniq/Total1
  
  p.mix = FitMixPoisson(dupStat)
  lambda = p.mix$lambda * sum(unlist(lapply(1:length(p.mix$p), function (i) i*p.mix$p[[i]])))
  Total = Nuniq /(1 - MixPoisson(0, p.mix$lambda, p.mix$p))
  Summary[sample_name,"Mix.Lambda"] = lambda
  Summary[sample_name,"Mix.total"] = round(x = Total, digits = 0)
  Summary[sample_name,"Mix.yield"] = Nuniq/Total
  Summary[sample_name,"Uniq.reads"] = Nuniq
  print(Summary)
  return (Summary)
}

# estimate ChIP yield
newCalculateChIPYield <- function(sample_name = NULL, 
                                  cfDNA_concentration = NULL,
                                  modification = NULL,
                                  y_params = yield.Params()) {
  RFile = paste0(SampleDIR, sample_name, ".rdata")
  dat = readRDS(RFile)
  catn(sample_name, ": computing QC")
  QC.params = cfChIP.countQC(dat,GR = QC.bed)
  
  cfDNA = cfDNA_concentration
  n_cells = round(cfDNA * 1e-9 / y_params$genome_w, digits = 0)
  n_genomes = n_cells * 2
  
  wins = readRDS(paste0(SetupDIR, modification, "/Windows.rds"))
  
  if (modification == "H3K4me1") { 
    signal_windows = wins[which(wins$type == "Enhancer")] 
    signal_background_rate = as.numeric(QC.params["%Background Enhancer"]) * 1e-2 # percent
    n_signal_nucs = y_params$blood_K4me1_nucs
  }
  else if (modification == "H3K4me2") {
    signal_windows = wins[which(wins$type == "TSS" | wins$type == "TSS ensembl" | wins$type == "TSS UCSC" | wins$type == "Enhancer")] 
    signal_background_rate = as.numeric(QC.params["%Background TSS"]) * 1e-2 # percent # TODO change
    n_signal_nucs = y_params$blood_K4me2_nucs  # TODO - find more accurate  value for K4me2. 
  }
  else if (modification == "H3K4me3") {
    signal_windows = wins[which(wins$type == "TSS" | wins$type == "TSS ensembl" | wins$type == "TSS UCSC")]
    signal_background_rate = as.numeric(QC.params["%Background TSS"]) * 1e-2 # percent
    n_signal_nucs = y_params$blood_K4me3_nucs
  }
  else if (modification == "H3K36me3") {
    signal_windows = wins[which(wins$type == "Gene")] 
    signal_background_rate = as.numeric(QC.params["%Background Gene"]) * 1e-2 # percent
    n_signal_nucs = y_params$blood_K36me3_nucs
  }
  
  
  catn(sample_name, ": calculating yield - global method")
  exp_tot_nucs = round(n_genomes * (y_params$genome_len / y_params$nucleosome_len), digits = 0)
  exp_marked_nucs = n_genomes * n_signal_nucs 
  est_params = CalculateSequencingYield(sample_name = sample_name, modification = modification)
  if (is.null(est_params)) {
    catn(sample_name, ": ERROR, can't estimate sequensint yield.")
    return (NULL)
  }
  aligned_frags = est_params$N.frags
  uniq_frags = round(est_params$N.uniq / y_params$sample_valume, digits = 0)
  est_uniq_frags = round(est_params$Mix.total / y_params$sample_valume, digits = 0) 
  est_seq_yield = uniq_frags / est_uniq_frags
  
  signal_rate = as.numeric(QC.params["%Signal Total"]) * 1e-2 # percent  TODO make sure this applies to all modifications!!!
  background_rate = as.numeric(QC.params["%Background Total"]) * 1e-2 # percent
  est_signal_frags = round(est_uniq_frags * signal_rate, digits = 0)
  est_background_frags = round(est_uniq_frags * background_rate, digits = 0)
  global_yield = round((est_signal_frags / exp_marked_nucs), digits = 8) 
  non_specific_yield = round((est_background_frags / exp_tot_nucs), digits = 8) 
  total_yield = round((global_yield + non_specific_yield), digits = 8) 
  
  # second method to calculate yield
  catn(sample_name, ": calculating yield - local method")
    
  high_cov = as.numeric(quantile(max(dat$Cov[signal_windows]), y_params$quantile_cutoff))
  local_yield = ((1 - (signal_background_rate)) * high_cov / n_genomes) 
  
  # background_windows = wins[which(wins$type == "background")]
  # bg_cov = mean(dat$Cov[background_windows])
  # local_yield_non_specific = bg_cov / n_genomes
  local_yield_non_specific = NaN
  
  yield_ratio = round(global_yield / local_yield, digits = 2)
  
  yield_params = data.frame(n.cells = c(), n.genomes = c(), aligned.fragments = c(), uniq.fragments = c(), lib.size = c(), 
                    seq.yield = c(), exp.nucs = c(), exp.marked_nucs = c(), signal.rate = c(), signal.frags = c(),
                    global.yield = c(), bg.rate = c(), bg.frags = c(), non_spc.yield = c(), tot.yield = c(),
                    high.cov = c(), sig.bg.rate = c(), local.yield = c(), local.yield.non_spc = c(), yields.ratio = c())
  
  yield_params[sample_name, "n.cells"] = n_cells
  yield_params[sample_name, "n.genomes"] = n_genomes
  yield_params[sample_name, "aligned.fragments"] = aligned_frags
  yield_params[sample_name, "uniq.fragments"] = uniq_frags
  yield_params[sample_name, "lib.size"] = est_uniq_frags
  yield_params[sample_name, "seq.yield"] = est_seq_yield
  yield_params[sample_name, "exp.nucs"] = exp_tot_nucs
  yield_params[sample_name, "exp.marked_nucs"] = exp_marked_nucs
  yield_params[sample_name, "signal.rate"] =  signal_rate 
  yield_params[sample_name, "siganl.frags"] = est_signal_frags
  yield_params[sample_name, "global.yield"] = global_yield
  yield_params[sample_name,"bg.rate"] = background_rate
  yield_params[sample_name,"bg.frags"] = est_background_frags
  yield_params[sample_name,"non_spc.yield"] = non_specific_yield
  yield_params[sample_name,"tot.yield"] = total_yield
  yield_params[sample_name,"high.cov"] = high_cov
  yield_params[sample_name,"sig.bg.rate"] = signal_background_rate
  yield_params[sample_name,"local.yield"] = local_yield
  yield_params[sample_name,"local.yield.non_spc"] = local_yield_non_specific
  yield_params[sample_name,"yields.ratio"] = yield_ratio
  
  return (yield_params)
}

yield.Params <- function() {
  list( 
    genome_len = 3.3e9,
    sample_valume = 2, # 2ml per sample
    quantile_cutoff = 0.95, # cutoff of high peak coverage (see below)
    genome_w = 6.6e-12, # molecular weight of 1 genome # CHECK. see http://www.bio.net/bionet/mm/methods/1999-December/080037.html
    nucleosome_len = 200,
    low_cutoff = 2, # cutoff of signal 
    blood_K4me3_nucs = 192015,
    blood_K4me1_nucs = 648193,
    blood_K4me2_nucs = 78384330,
    blood_K36me3_nucs = 782396
  )
}

if(0) {
# length of all peaks (cat <sample>-H3K4me3.narrowPeak | awk '{sum += $3-$2}END{print sum}')
genome_len = 3.3e9
nucleosome_len = 200
# K4me3 
monocyte_K4me3 = 29312220 # E029
neutrophil_K4me3 = 51333230 # E030
pbmc_K4me3 = 34563641 # E062
liver_K4me3 = 53883151 # E066
vascular_K4me3 = 42315250 # E122
#  number of H3K4me3 nucleosomes in common blood cells
blood_K4me3_nucs = mean(c(monocyte_K4me3, neutrophil_K4me3, pbmc_K4me3)) / nucleosome_len
genomic_K4me3 = blood_K4me3_nucs * nucleosome_len / genome_len # percent of genome with K4me3 signal

# K4me1
monocyte_K4me1 = 132298213 # E029
neutrophil_K4me1 = 126979280 # E030
pbmc_K4me1 = 38846607 # E062
#  number of H3K4me1 nucleosomes in common blood cells
blood_K4me1_nucs = mean(c(monocyte_K4me1, neutrophil_K4me1)) / nucleosome_len # pbmc is very different
genomic_K4me1 = blood_K4me1_nucs * nucleosome_len / genome_len # percent of genome with K36me3 signal

# K4me2
#  median of H3K4me2 nucleosomes in all roadmap cells (no samples of common blood cells)
blood_K4me2_nucs = 78384330
genomic_K4me2 = blood_K4me2_nucs * nucleosome_len / genome_len # percent of genome with K4me2 signal

# K36me3
monocyte_K36me3 = 172679180 # E029
neutrophil_K36me3 = 147868308 # E030
pbmc_K36me3 = 148890417 # E062
#  number of H3K36me3 nucleosomes in common blood cells
blood_K36me3_nucs = mean(c(monocyte_K36me3, neutrophil_K36me3, pbmc_K36me3)) / nucleosome_len
genomic_K36me3 = blood_K36me3_nucs * nucleosome_len / genome_len # percent of genome with K36me3 signal

}
source(paste0(SourceDIR,"SplitHeatmap.R"))

plotPVal = function( PVals, Counts, 
                     max.log.P.val = 20, 
                     max.Count = 20,
                     Q.val.threshold = 0.001,
                     PValueAdjust = TRUE,
                     rowOrder = NA, colOrder = NA, 
                     splitMap = FALSE) {
  Q.Values = PVals
  if( PValueAdjust ) 
    Q.Values[] =  -log10(p.adjust(exp(Q.Values),method="fdr"))
  
  Q.Values[is.infinite(Q.Values)] = max.log.P.val
  Q.Values[Q.Values>max.log.P.val ] = max.log.P.val
  Q.Values[Q.Values < -log10(Q.val.threshold)] = 0
  
  if( splitMap ) {
    p = plotSplitHeatmap(dim(PVals)[1],dim(PVals)[2],Q.Values, Counts, 
                         Upper.name = "PValue", 
                         Lower.name = "Counts",
                         Upper.max = max.log.P.val,
                         Lower.max = max.Count,
                         colorBins = min(max(max.log.P.val, max.Count)+1,20),
                         Order.Cols = colOrder,
                         Order.Rows = rowOrder,
                         Sep.Width = .5)
  } else {
    ## TBD
    catn("Missing option")
  }
    
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5 ))
  if( 0 ) {
    df = cbind(data.frame(Tissue = row.names(Q.Values)),data.frame(Q.Values))
    m = melt(df)
    m$count = pmin(max.Count,as.vector(Counts))
    #m$variable = factor(as.character(m$variable),
    #                   levels= colnames(Q.Values)[order(-colSums(Q.Values))])
    df = data.frame(Sample = m$Tissue, Tissue = m$variable, Qvalue = m$value, Count = m$count)
    pQ = ggplot(df, aes(Tissue,Sample))+geom_tile(aes(fill=Qvalue))
    pQ = pQ + scale_fill_gradient(low = "white",high = "steelblue")
    pC = ggplot(df, aes(Tissue,Sample))+geom_tile(aes(fill=Count))
    pC = pC + scale_fill_gradient(low = "white",high = "orange")
    list( QVals = Q.Values, plot.Q = pQ, plot.C = pC)
  }
  list( QVals = Q.Values, plot.Q = p, PVals = PVals, Counts = Counts)
}

GeneSigsEval = function(sig, method, Quantile) {
  if( length(sig) > 1 ) {
    sapply( colnames(GeneCounts), function(s)
      computePValue(GeneCounts[sig,s], GeneBackground[sig,s], method, Quantile)[1])
  } else {
    bg = GeneBackground[sig,]
    fg = GeneCounts[sig,]
    ppois(fg,bg,lower.tail = FALSE, log.p = TRUE)
  }
}

GeneSigsPVal = function(Sigs, method = "Poisson", Quantile = 0.1) {
  A = sapply(Sigs, function(sig) GeneSigsEval(sig, method, Quantile))
  rownames(A) = colnames(GeneCounts)
  A
}

GeneSigsCounts = function(Sig)  sapply(Sig, 
                                       function(sig) 
                                         if( length(sig) > 1 ) colSums(GeneCounts.QQnorm[sig,])/length(sig) 
                                       else GeneCounts.QQnorm[sig,])

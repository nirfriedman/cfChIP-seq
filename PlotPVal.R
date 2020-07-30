source(paste0(SourceDIR,"SplitHeatmap.R"))
source(paste0(SourceDIR,"DotHeatMap.R"))
source(paste0(SourceDIR,"PlainHeatmap.R"))

plotPVal = function( PVals, Counts, 
                     max.log.P.val = 20, 
                     max.Count = 20,
                     Q.val.threshold = 0.001,
                     PValueAdjust = TRUE,
                     rowOrder = NA, colOrder = NA, # can be an order, "cluster", "sum" or "NA"
                     Filter.Rows = FALSE,
                     Filter.Cols = FALSE,
                     Filter.threshold = 6,
                     Lower.name = "Counts",
                     Upper.name = "q Value", 
                     Pal = "default",
                     MapType = "dot",
                     MaxNameLength = 50,
                     DarkBackground = FALSE) {
  
  rownames(PVals) = strtrim(rownames(PVals), MaxNameLength)
  rownames(Counts) = strtrim(rownames(Counts), MaxNameLength)
  colnames(PVals) = strtrim(colnames(PVals), MaxNameLength)
  colnames(Counts) = strtrim(colnames(Counts), MaxNameLength)
  Q.Values = PVals
  if( PValueAdjust ) 
    Q.Values[] =  -log10(p.adjust(exp(-log(10)*Q.Values),method="fdr"))
  
  Q.Values[is.infinite(Q.Values)] = max.log.P.val
  Q.Values[Q.Values>max.log.P.val ] = max.log.P.val
  Q.Values[Q.Values < -log10(Q.val.threshold)] = 0
  Lower.Pal = brewer.pal(9,"YlOrRd")[2:9]
  Lower.Neg.Pal = brewer.pal(9,"Greens")[1:9]
  if( Pal == "BlueOrange") {
    Lower.Pal = brewer.pal(9,"Oranges")[2:8]
    Lower.Neg.Pal = brewer.pal(9,"Blues")[2:8]
  }
  if( Pal == "GreenRed") {
    Lower.Pal = brewer.pal(9,"Reds")[2:8]
    Lower.Neg.Pal = brewer.pal(9,"Greens")[2:8]
  }
  if( Pal == "DarkGreenRed") {
    Lower.Neg.Pal = c("#043303", "#0D5808", "#10650A", "#2BC91F", "#38FC2A")
    Lower.Pal = c("#320001","#640005", "#95000A", "#C70011", "#F70018")
    DarkBackground = TRUE
  }
  if( Pal == "BlueYellow") {
    Lower.Pal = brewer.pal(9,"Oranges")[2:6]
    Lower.Neg.Pal = brewer.pal(9,"Blues")[2:8]
  }
  if( Pal == "DarkBlueOrange") {
    Lower.Neg.Pal = c("#022232", "#064465", "#0E6596", "#1787C9", "#1FA9FB")
    Lower.Pal = c("#321102","#642205", "#96320C", "#C84314", "#FA541B")
    DarkBackground = TRUE
  }
  if( Filter.Rows ) {
    I.row = rownames(Q.Values) == "ruler" | apply(Q.Values,1, max, na.rm = TRUE) >= Filter.threshold
  } else
    I.row = rep(TRUE, nrow(Q.Values))
  
  if( Filter.Cols ) {
    I.col = colnames(Q.Values) == "ruler" | apply(Q.Values,2, max, na.rm = TRUE) >= Filter.threshold
  } else
    I.col = rep(TRUE, ncol(Q.Values))

  if( any(I.row) && any(I.col) ) {
    Q.Values = Q.Values[I.row,I.col, drop=FALSE]
    Counts = Counts[I.row, I.col, drop=FALSE]
    
    if(is.na(rowOrder) ) {
      rowOrder = rev(1:nrow(Counts))
    } else {
      if( length(rowOrder) == 1 ) {
        if( rowOrder == "cluster"){
          if (nrow(Counts)==1){
            rowOrder=1
          }else {
            Counts1 = Counts
            Counts1[Q.Values == 0] = 0
            rowOrder = hclust(dist(Counts1))$order
          }
        } else 
          if( rowOrder == "sum") {
            rowOrder = order(rowSums(Q.Values))
          }
      }
    }
    if(is.na(colOrder) ) {
      colOrder = 1:ncol(Counts)
    } else {
      if( length(colOrder) == 1 ) {
        if( colOrder == "cluster"){
          if(ncol(Counts)==1){
            colOrder=1
          }else {
            Counts1 = Counts
            Counts1[Q.Values == 0] = 0
            colOrder = hclust(dist(t(Counts1)))$order
          }
        } else 
          if( colOrder == "sum") {
            colOrder = order(colSums(Q.Values))
          }
      }
    }
    
    if( MapType == "split" ) {
      p = plotSplitHeatmap(nrow(Counts),ncol(Counts),
                           Q.Values, Counts, 
                           Upper.name = Upper.name, 
                           Lower.name = Lower.name,
                           Upper.max = max.log.P.val,
                           Lower.max = max.Count,
                           colorBins = min(max(max.log.P.val, max.Count)+1,20),
                           Order.Cols = colOrder,
                           Order.Rows = rowOrder,
                           Sep.Width = .5)
    } else 
      if( MapType == "dot" ) {
        ZeroColor = "white"
        if(DarkBackground)
          ZeroColor = "black"
        p = plotDotHeatmap(nrow(Counts),ncol(Counts),
                           Q.Values, Counts, 
                           Upper.name = Upper.name, 
                           Lower.name = Lower.name,
                           Upper.max = max.log.P.val,
                           Lower.max = max.Count,
                           colorBins = min(max(max.log.P.val, max.Count)+1,20),
                           Order.Cols = colOrder,
                           Order.Rows = rowOrder,
                           ZeroColor = ZeroColor,
                           Lower.Pal = Lower.Pal,
                           Lower.Neg.Pal = Lower.Neg.Pal)
      } else
        if( MapType == "plain") {
          p = plotPlainHeatmap(nrow(Counts),ncol(Counts),
                               Counts, 
                               Lower.name = Lower.name,
                               Lower.max = max.Count,
                               colorBins = min(max(max.log.P.val, max.Count)+1,20),
                               Order.Cols = colOrder,
                               Order.Rows = rowOrder,
                               Sep.Width = .5)
        } else
          catn(paste0("Unknown map type `", MapType,"`"))
    
    p = p + theme(axis.text.x.top = element_text(angle = 90, 
                                                 hjust = 0, vjust=0.5 ))
    if( DarkBackground )
      p = p + theme(panel.background = element_rect(fill = "black"),
                    panel.grid = element_line(color = "gray20"))
  } else {
    warning("PlotPval - empty matrix after filter")
    p = NULL
  }
  
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
    ppois(fg-1,bg,lower.tail = FALSE, log.p = TRUE)
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

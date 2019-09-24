require(cowplot)
require(reshape2)


CollectMeta = function( Cov, # RleList
                        Regions, # GRanges
                        Width,Offset,
                        WindowSize = 10,
                        Norm = 1 )
{
  Regions = resize(Regions, width=Offset + 0.5*width(Regions), fix="end")
  Regions = resize(Regions, width=Width, fix="start")
  
  print(paste("CollectMeta",Regions$name[1], length(Regions), Norm))
  
  ws = IRanges(start=seq(1,Width-1, by=WindowSize), width=WindowSize)
  n = length(Regions)
  m = ceiling(Width/WindowSize)
  X = matrix(0, nr = n, nc = m, dimnames = list(Regions$name, start(ws)))
  if( 0 ) {
    st = start(Regions)
    en = end(Regions)
    I = as.logical(strand(Regions) == "-")
    temp = en[I]
    en[I] = st[I]
    st[I] = temp
    chr = as.character(chrom(Regions))
    #  Y =  sapply(1:n, 
    #              function(i) as.integer(Cov[[chr[i]]][st[i]:en[i]]))
    X[,] = t(sapply(1:n, 
                    function(i) aggregate(Cov[[chr[i]]][st[i]:en[i]], ws, sum)))
  }
  I = as.logical(strand(Regions) == "-")
  X[,] = t(sapply(Cov[Regions], function(x) aggregate(x, ws, sum)))
  if( any(I) )
    X[I,1:m] = X[I,m:1]
  
  print(quantile(X,c(0.99, 0.9, 0.5)))
  X = X*Norm/WindowSize
  print(quantile(X,c(0.99, 0.9, 0.5)))
  X = X[order(rowSums(X)),]
}


ggPlotCovergeGroups = function(A, pos=5000, label = "TSS",
                               Kbytes = 100, 
                               xtickSpace = 5,
                               extraSpace = 5-1,
                               main = "", ylim = NULL)
{
  if( is.null(ylim) ) {
    M = max(0,max(sapply(A,max),na.rm=TRUE)*1.1)
    m = min(0,min(sapply(A,min), na.rm = TRUE) *1.1)
    ylim=c(m,M)
  }
  w = length(A[[1]])
  if( is.nan(w) || is.infinite(w))
    w = 20000
  print(w)
  
  p = ggplot()
  p = p + theme(legend.position = c(0.95,0.95), legend.justification = c("right", "top") )
    #theme(legend.position = c(1,0), legend.justification=c(1,0)) 
  p = p + scale_color_brewer( palette = "Set1",
                              guide = guide_legend(title=NULL, reverse = FALSE))
  
  # plot lines
  p = p+geom_vline(xintercept = 0, colour="gray", size=1)
  if( ylim[[1]] < 0 && ylim[[2]] > 0 )
    p = p+geom_hline(yintercept = 0, colour="gray", size=1)
  
 

  for( i in 1:length(A)) {
    p = p + geom_line( data=data.frame(x = (1:w)-Kbytes*pos/1000, y = A[[i]], color = rep(names(A)[[i]], 1)), 
                       aes(x = x, y = y, colour = color),size = .5, 
                       show.legend = TRUE)
  }

  ipos = as.integer(pos/(xtickSpace*1000))*Kbytes*xtickSpace
  iwidth = as.integer(w/Kbytes)*Kbytes
  labelat = seq(-ipos,iwidth, xtickSpace*Kbytes)
  labelval = paste0(labelat/Kbytes, "Kb")
  labelval[labelval == "0Kb"] = label
  xlim = c(- Kbytes*pos/1000,iwidth - Kbytes*(pos/1000+extraSpace)) 
  
#  print("ggPlotCovergeGroups")
#  print(labelat)
#  print(labelval)
  p = p+scale_x_continuous(name="", breaks = labelat, labels = labelval, limits = xlim, expand=c(0,0) )
  p = p + scale_y_continuous(name = "", expand=c(0,0), limits = ylim)
  p = p + labs(title=main) 
  p = p + theme(text = element_text(colour = "black",size = 8))
  p = p + theme(axis.text = element_text(colour = "black",size = 8))
  
  return(p)
}



ggPlotHeatMap <- function( A, xlab = "", ylab = "Y", zlab = "Coverage", title = "", 
                         zlim = NULL,
                         offset=5000, 
                         label = "TSS",
                         Kbytes = 100, 
                         color = "red",
                         bg.color = "white",
                         xtickSpace = 5) {
  m = do.call(rbind, A)
  n = ncol(m)
  rownames(m) <- 1:nrow(m)
  colnames(m) <- 1:n
  ioff =  Kbytes*offset/1000
  breaks = seq(0,n,xtickSpace*Kbytes) - ioff
  labels = (breaks)/Kbytes
  labels = paste0(labels, "Kb")
  labels[labels == "0Kb"] = label
  
  ybreaks = cumsum(sapply(A,nrow))
  ylabels = names(A)
  
  df  = melt(m, varnames = c( "y", "x"))
  df$x = df$x - ioff
  t1 = quantile(m, probs = .99, na.rm = T)
  t2 = quantile(m, probs = .05, na.rm = T)
  if( t2 > 0 )
    t2 = 0
  # print(c(t2,t1))
  if( is.null(zlim) ) 
    zlim = c(t2,t1)
  df$value[df$value < zlim[[1]]] = zlim[[1]]
  df$value[df$value > zlim[[2]]] = zlim[[2]]
  #  p = ggplot(df, aes(xlab, ylab, fill=value))
  p = ggplot(df, aes(x=x,y=y, fill=value))
  #  p = p + geom_raster(color = "white")
  p = p + geom_raster()
  p = p+geom_vline(xintercept = 0, colour="gray", size=1)
  for( h in ybreaks[1:(length(ybreaks)-1)] )
    p = p+geom_hline(yintercept = h, colour="gray", size=1)
  ybreaks = (ybreaks+c(0,ybreaks[1:(length(ybreaks)-1)]))/2
  p = p + scale_fill_gradient(low=bg.color, high=color, na.value="black", limits = zlim )

#  print("ggPlotHeatMap")
#  print(breaks)
#  print(labels)
  p = p + scale_x_continuous( breaks = breaks, labels = labels, limits=c( - ioff,n - ioff),expand=c(0,0) )
#  print(ybreaks)
#  print(ylabels)
#  p = p + scale_y_continuous( breaks = ybreaks, labels = ylabels, expand=c(0,0) )
  p = p + labs(x = xlab, y = ylab, fill=zlab, title = title)
  p = p + theme_minimal()
  p = p + theme(text = element_text(colour = "black",size = 8))
  p = p + theme(axis.text = element_text(colour = "black",size = 8))
  #  p = p + theme( legend.position = "none",
  #                axis.ticks = element_blank()
  #  )
  return(p)
}

PlotMeta = function(Cov,  
                    PlotList,
                    WindowSize = 25, 
                    Norm = 1) {

  Cov = fixCoverage(Cov)
  print(paste("Norm", Norm))
  p.list = lapply(PlotList, function(pl) {
    RegionGroups = split(pl$BED, pl$BED$name)
    RegionGroups = RegionGroups[order(names(RegionGroups),decreasing = TRUE)]
    print(names(RegionGroups))
    names(RegionGroups) = sub("^[[:lower:]].","", names(RegionGroups))
    Ms = lapply(RegionGroups, function(g) CollectMeta(Cov, g, pl$Width, pl$Offset, WindowSize, Norm))
  
    Ls = lapply(Ms, colMeans)
    if( !is.null(pl$Max) && pl$Max > 0 ) {
      clim = pl$Max
    } else {
      Ls.max = max(max(unlist(Ls)))
      if( Ls.max > 3) {
        clim = c(0, 5*ceiling(Ls.max/5))
      } else
        clim = c(0, ceiling(Ls.max))
    }
    print(clim)
    print(pl$Color)
    p.meta = ggPlotCovergeGroups(Ls, 
                                 pos = pl$Offset, label = pl$Label, xtickSpace = pl$Tick/1000,
                                 Kbytes = 1000/WindowSize, ylim = clim, extraSpace = 0)
    p.HS = ggPlotHeatMap(Ms, offset = pl$Offset,label = pl$Label,
                       Kbytes = 1000/WindowSize, ylab = "", color = pl$Color, bg.color = pl$BGColor,
                       zlim = clim)
    p.HS = p.HS + guides(fill=FALSE) + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  
    list(p.meta, p.HS)
  })
  p.list = c(lapply(p.list, function(l) l[[1]]), lapply(p.list, function(l) l[[2]]))
  p = plot_grid(plotlist = p.list, ncol = length(PlotList), align="vh",  rel_heights = rep(c(1,2), length(PlotList)),
                axis="lrtb")
  return(p)
}

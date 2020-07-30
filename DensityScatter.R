library(ggplot2)
library(MASS)

#TODO - https://slowkow.com/notes/ggplot2-color-by-density/
#dat$density <- get_density(dat$x, dat$y, n = 100)
#ggplot(dat) + geom_point(aes(x, y, color = density)) + scale_color_viridis()


# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


DensityScatter <- function( c1, c2, alpha = 0.5, title = " ", xlabel = "x", ylabel = "y", 
                            threshold = 0.995, xlim=NULL, ylim=NULL, 
                            Groups = NULL, galpha=1,
                            Value = NULL,
                            logscale=FALSE, 
                            log1p = FALSE,
                            contour = FALSE,
                            diagonal=FALSE, horizontal=NULL, coordeq = FALSE, 
                            smooth = FALSE, linearfit = FALSE,
                            jitter = FALSE, stat = FALSE,
                            bin = FALSE, binwidth = 0.1,
                            pointsize = 1,
                            Cor = FALSE,
                            GroupContour = FALSE ) {
  
  if( logscale && log1p ) {
    cat("Cannot choose both logscale and log1p", "\n")
    log1p = FALSE
  }
  
  # remove outliers
  t1 = quantile(c1,threshold, na.rm=TRUE)
  t2 = quantile(c2,threshold, na.rm=TRUE)
  
  I = (c1 < t1) & (c2 < t2) & !is.na(c1) & !is.na(c2)
  print(length(which(I)))

  groups = vector(length= length(c1))
  groups[] = NA
  alphas = vector(length= length(c1))
  alphas[] = alpha
  
  DoColor = 0
  if( !is.null(Groups) ) {
    DoColor = 1
    Names = labels(Groups)
    for( i in 1:length(Groups) ) {
      groups[Groups[[i]]] = Names[[i]];
#      print(c(i, length(which(Groups[[i]]))))
      alphas[Groups[[i]]] = galpha
    }
  } else {
    if( !is.null(Value) )
    {
#      groups = rank(Value,ties.method = "random")
      groups = Value
      alphas[] = galpha
      DoColor = 1;      
    } 
  }
  
  # set up plot data
  df = data.frame( x=c1[I],y=c2[I], g = groups[I], a = alphas[I])
  
  if( DoColor )  {
    p = ggplot(df,aes(x=x,y=y,colour=g,alpha=a)) 
    if(!is.null(Groups) )
      p = p + scale_colour_brewer(palette = "Set1", na.value = "gray")
    if( !is.null(Value))
      p = p + scale_colour_gradient2(mid = "gray75", high="darkred",low="darkblue")
  } else {
    if( logscale ) {
      df$density = get_density(log10(df$x), log10(df$y), n = 100)
    } else 
      if(log1p) {
        df$density = get_density(log10(1+df$x), log10(1+df$y), n = 100)
      } else
        df$density = get_density(df$x, df$y, n = 100)
      
    p = ggplot(df,aes(x=x,y=y, color = density))
    p = p + scale_color_viridis_c(guide = FALSE)
  }
  # labels
  p = p + labs(title=title, x = xlabel, y=ylabel) 
  
  if( diagonal ) 
    p = p + geom_abline(slope=1,intercept=0, color="darkgray", size=2)

  if( !is.null(horizontal) )
    p = p+geom_hline(yintercept = as.numeric(horizontal), color="orange", size=2)

  pointcolor = "blue"

  if( stat ) {
    # still experimental
    p = p + stat_density2d( aes(fill = ..density.. ), 
                            geom="tile", 
                            contour = F, 
                            colour = "black" )
    pointcolor = "red"
  }
  
  
  if( bin ){
    p = p+ geom_bin2d(binwidth=binwidth) + scale_fill_distiller(palette = "YlOrBr", values = c(0.0,0.2,0.5,1))
  } else {
    if( !DoColor ) { 
      if( jitter ) 
        p = p+geom_jitter(alpha = alpha, #colour=pointcolor, 
                          height = .5, width = .5,
                          shape=16, size=pointsize )
      else
        p = p+geom_point(alpha = alpha, #colour=pointcolor,
                         shape=16, size=pointsize )
    } else {
      if( jitter )
        p = p+geom_jitter(width = .5, height = .5,shape=16,size=pointsize)
      else
        p = p+geom_point(shape=16,size=pointsize,aes(x=x,y=y,colour=g, alpha=a));
    }
  }
  
  if( smooth || linearfit ) {
    if( logscale ) {
      logx = log(c1[I])
      logy = log(c2[I])
    } else {
      logx = (c1[I])
      logy = (c2[I])
    }
    J = !is.infinite(logx) & !is.infinite(logy)
    logdf = data.frame( x=logx[J],y=logy[J] )

    if( linearfit ) {
      lo <- rlm(y~x, logdf, method="MM")
      if( logscale ) {
        expx = exp(logx[J])
        expy = exp(lo$fitted)
      } else {
        expx = (logx[J])
        expy = (lo$fitted)
      }
      p = p+geom_line(data = data.frame( x=expx,y=expy), aes(x=x,y=y), color="orange", size=2)
    }
    if( smooth ) {
      #    sy = stat_smooth(data = logdf,aes(x=x,y=y), method="loess", alpha=12)
      lo <- loess(y ~ x, logdf)
    
      if( logscale ) {
        expx = exp(logx[J])
        expy = exp(lo$fitted)
      } else {
        expx = (logx[J])
        expy = (lo$fitted)
      }
      p = p+geom_line(data = data.frame( x=expx,y=expy), aes(x=x,y=y), color="red", size=2)
    }
  }
  
  if( logscale ) {
      p = p+scale_x_log10()
      p = p+scale_y_log10()
  } else 
    if(log1p) {
      p = p + scale_x_continuous(trans="log1p", breaks = c(0, 1,3, 10,30,100,300, 1000))  
      p = p + scale_y_continuous(trans="log1p",breaks = c(0, 1,3,10,30, 100,300, 1000))  
    } 
  if( !is.null(xlim))
    p = p+xlim(xlim)
  if( !is.null(ylim))
    p = p+ylim(ylim)
  
  if( coordeq )
    p = p + coord_fixed(ratio=1) 

  if( !is.null(Groups) && GroupContour ) 
    contour = FALSE
 
  if( contour ) 
    p = p+geom_density2d(colour="black", show.legend = FALSE)

  
  if( !is.null(Groups) && GroupContour ) {
    p = p+geom_density2d(data = df[!is.na(df$g),],aes(x=x,y=y), 
                         color="black", show.legend = FALSE)

    if(0){
      for( i in seq_along(Groups) ) {
        # p+geom_density_2d(aes(x=x,y=y,colour=g),na.rm = T,)
        dd = df[df$g == i & !is.na(df$g),]
        dens = kde2d(dd$x, dd$y)
        densf = data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
        p=p+geom_contour(densf, aes(x=x,y=y,z=z))
      }
    }
  }
  if( Cor ) {
    r = cor(df$x,df$y, use="pairwise.complete.obs")
#    if( is.null(xlim) ) 
    labeltext = sprintf("R^2 = %.3f", r**2)
     p = p + labs(subtitle=labeltext)
#    p = p + annotate(geom = "text", label = labeltext, x = t1*.1, y=t2*.85, hjust = "top", vjust="right", parse=TRUE)
  }
  p = p + theme_minimal()
  
  return(p)
}


Plot3D <- function( x,y,z, xl, yl, zl, alpha = 1, threshold = .995, logscale = FALSE ) {
  T1 = quantile(x,threshold, na.rm=TRUE)
  T2 = quantile(y,threshold, na.rm=TRUE)
  T3 = quantile(z,threshold, na.rm=TRUE)
  t1 = quantile(x,1-threshold, na.rm=TRUE)
  t2 = quantile(y,1-threshold, na.rm=TRUE)
  t3 = quantile(z,1-threshold, na.rm=TRUE)
  
  I = !is.na(x) & !is.na(y) & !is.na(z) & x < T1 & x > t1 & y < T2 & y > t2 & z < T3 & z > t3
  if( length(which(I)) < 10)
    return( ggplot(data.frame(x=1:10,y=1:10),aes(x = x, y=y) ) + geom_line())
  x = x[I]
  y = y[I]
  z = z[I]
  o = order(z)
  df = data.frame( x = x[o], y = y[o], z = z[o])
#  p = ggplot(df, aes(x = x, y=y, color=z, alpha=z)) + geom_point() #(alpha=alpha)
  p = ggplot(df, aes(x = x, y=y, color=z, alpha=z)) + geom_point() #(alpha=alpha)
  if( logscale )
    p = p+ scale_x_log10()+scale_y_log10()
  p = p+ labs(title=zl, x = xl, y=yl) 
  #  p = p+scale_colour_gradientn(colours = topo.colors(10))
  #  p = p+scale_colour_gradientn(colours = c("gray", "pink", "red", "red1", "red2", "red3", "darkred"),
  #                               values = c(0,0.2, 0.6, 0.7, 0.8, 0.9, 1))
  p = p+scale_color_gradientn(colours = c("gray",  "red3"))
  p = p+geom_jitter(width = .5, height = .5)
  return(p)
}


Plot3DSmooth <- function( x,y,z, xl, yl, zl, alpha = 1, threshold = .995, logscale = FALSE, N = 5000 ) {
  T1 = quantile(x,threshold, na.rm=TRUE)
  T2 = quantile(y,threshold, na.rm=TRUE)
  T3 = quantile(z,threshold, na.rm=TRUE)
  t1 = quantile(x,1-threshold, na.rm=TRUE)
  t2 = quantile(y,1-threshold, na.rm=TRUE)
  t3 = quantile(z,1-threshold, na.rm=TRUE)
  
  I = !is.na(x) & !is.na(y) & !is.na(z) & x < T1 & x > t1 & y < T2 & y > t2 & z < T3 & z > t3
  if( length(which(I)) < 10)
    return( ggplot(data.frame(x=1:10,y=1:10),aes(x = x, y=y) ) + geom_line())
  x = x[I]
  y = y[I]
  z = z[I]

  if( length(x) > N)
    I = sample(1:length(x), N)
  else
    I = rep(TRUE, length(x))
  l = loess(z~x*y,data.frame(x=x[I],y=y[I],z=z[I]))
  df = data.frame(x=x[I],y=y[I],z=l$fitted)
  p = ggplot(df, aes(x = x, y=y, color=z, alpha=alpha)) + geom_point() #(alpha=alpha)
  if( logscale )
    p = p+ scale_x_log10()+scale_y_log10()
  p = p+ labs(title=zl, x = xl, y=yl) 
  #  p = p+scale_colour_gradientn(colours = topo.colors(10))
  #  p = p+scale_colour_gradientn(colours = c("gray", "pink", "red", "red1", "red2", "red3", "darkred"),
  #                               values = c(0,0.2, 0.6, 0.7, 0.8, 0.9, 1))
  p = p+scale_color_gradient2(low = "green", high = "red3", mid = "gray50", midpoint = 0)
  return(p)
}

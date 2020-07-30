library(ggplot2)
library(ggforce)
library(RColorBrewer)

# test values
if( 0 ) {
  nr = 2 
  nc = 3 
  Upper = matrix( 1:6, nr = nr, nc = nc, dimnames = list(c("Alpha","Beta"),c("A","B","C"))) 
  Lower =  matrix( c(0,-2,3,4,-10,12), nr = nr, nc = nc, dimnames = list(c("Alpha","Beta"),c("A","B","C"))) 
  Upper.name = "up" 
  Lower.name = "down"
  Upper.max = NA
  Lower.max = NA
  Order.Rows = NA
  Order.Cols = NA
  Upper.Pal = brewer.pal(9,"Blues")
  Lower.Pal = brewer.pal(9,"Reds")
  Lower.Neg.Pal = brewer.pal(9,"Greens")
  colorBins = 10
  Labels.Rows = NA
  Labels.Cols = NA 
  Sep.Width = 2
  Sep.Color = "light gray"
}


plotPlainHeatmap = function( nr, nc, Lower, 
                             Lower.max = NA,
                             Order.Rows = NA,
                             Order.Cols = NA,
                             Lower.Pal = brewer.pal(9,"YlOrRd"),
                             Lower.Neg.Pal = brewer.pal(9,"Greens"),
                             colorBins = 20,
                             Labels.Rows = NA,
                             Labels.Cols = NA,
                             Lower.name = "value",
                             Sep.Width = 1,
                             Sep.Color = "lightgray"
                           ) {
  
  processValues = function(X, X.max) {
    X[X > X.max] = X.max
    X[X < -X.max] = -X.max
    
    a = (colorBins-1)/X.max
    Y = floor(a * t(abs(X)))
    Y = Y * sign(t(X))
    Y[1:length(X)]
  }
  
  if( any(is.na(Lower.max)) )
    Lower.max = max(abs(Lower))
  if( any(is.na(Labels.Rows)) )
    Labels.Rows = rownames(Lower)
  if( any(is.na(Labels.Cols)) )
    Labels.Cols = colnames(Lower)
  if( any(is.na(Order.Rows) ) )
    Order.Rows = 1:nr
  if( any(is.na(Order.Cols)) )
    Order.Cols = 1:nc
  
  # build dataset
  df = data.frame( value = processValues(Lower[Order.Rows,Order.Cols],Lower.max),
                   x = rep(1:nc, nr ),
                   y = rep(1:nr, each=nc),
                   stringsAsFactors = FALSE
        )

  # build color map
  m2 = floor((2*colorBins-1)/11)+1
 
  Colors = c( rev(colorRampPalette(Lower.Pal, space = "Lab")(colorBins-1)), "white",
              colorRampPalette(Lower.Neg.Pal, space = "Lab")(colorBins-1) 
              )
  
  names(Colors) = (colorBins-1):(-colorBins+1)
  breaks = c( (0:4)*m2+1,
              colorBins,
              2*colorBins-(4:0)*m2-1)
  
  Colors.Labels = c( formatC(Lower.max), rep("", 4), "0", rep("",4),  formatC(-Lower.max))
  p = ggplot(data = df, aes(x = x, y = y, fill=value)) 
  p = p+theme_minimal()
  p = p+geom_raster()
  p = p+geom_hline(color = Sep.Color, size = Sep.Width, yintercept = 0.5 +0:nr)
  p = p+geom_vline(color = Sep.Color, size = Sep.Width, xintercept = 0.5 +0:nc)
  p = p+scale_fill_gradientn(
    colors = rev(Colors),
    #          values = rev(as.numeric(names(Colors))),
    breaks = as.numeric(names(Colors)[breaks]), 
    labels = Colors.Labels,
    name = Lower.name,
    limits = c(-colorBins-1, colorBins-1),
    guide = guide_colorbar())
  p = p+scale_x_continuous(breaks = 1:nc, 
                           labels = as.character(Labels.Cols[Order.Cols]),
                           position = "top",
                           expand = c(0,0))
  p = p+scale_y_continuous(breaks = 1:nr, 
                           labels = as.character(Labels.Rows[Order.Rows]),
                           position = "right",
                           expand = c(0,0))
  p = p + labs(x = "", y = "")
 # p = p + theme(panel.grid =  element_blank())
  # p = p + coord_equal()
  p
}

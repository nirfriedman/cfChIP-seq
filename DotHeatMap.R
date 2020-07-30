library(ggplot2)
library(ggforce)
library(RColorBrewer)

# test values
if( 0 ) {
  nr = 2 
  nc = 3 
  Upper = matrix( 0:5, nr = nr, nc = nc, dimnames = list(c("Alpha","Beta"),c("A","B","C"))) 
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
}


plotDotHeatmap = function( nr, nc, Upper, Lower, 
                             Upper.max = NA,
                             Lower.max = NA,
                             Order.Rows = NA,
                             Order.Cols = NA,
                             Lower.Pal = brewer.pal(9,"YlOrRd"),
                             Lower.Neg.Pal = brewer.pal(9,"Greens"),
                           ZeroColor = "white",
                             colorBins = 20,
                             Labels.Rows = NA,
                             Labels.Cols = NA,
                           Upper.name = "",
                           Lower.name = "value",
                           Upper.offset = 2
                           ) {
  
  processValues = function(X, X.max) {
    X[X > X.max] = X.max
    X[X < -X.max] = -X.max
    
    a = (colorBins-1)/X.max
    Y = floor(a * t(abs(X)))
    Y = Y * sign(t(X))
    Y[1:length(X)]
  }
  
  if( any(is.na(Upper.max)) )
    Upper.max = max(Upper)
  if( any(is.na(Lower.max)) )
    Lower.max = max(abs(Lower))
  if( any(is.na(Labels.Rows)) )
    Labels.Rows = rownames(Upper)
  if( any(is.na(Labels.Cols)) )
    Labels.Cols = colnames(Upper)
  if( any(is.na(Order.Rows) ) )
    Order.Rows = 1:nr
  if( any(is.na(Order.Cols)) )
    Order.Cols = 1:nc
  
  # build datasets
  df = data.frame( upper.value = processValues(Upper[Order.Rows,Order.Cols],Upper.max),
                   value = processValues(Lower[Order.Rows,Order.Cols],Lower.max),
                   x = rep(1:nc, nr ),
                   y = rep(1:nr, each=nc),
                   stringsAsFactors = FALSE
        )

  I = which(df$upper.value > 0)
  df$upper.value[I] = df$upper.value[I]+Upper.offset
  Upper.max = Upper.max + Upper.offset
  
  # build color map
  m2 = floor((2*colorBins-1)/11)+1
 
  Colors = c( rev(colorRampPalette(Lower.Pal, space = "Lab")(colorBins-1)), ZeroColor,
              colorRampPalette(Lower.Neg.Pal, space = "Lab")(colorBins-1) 
              )
  
  names(Colors) = (colorBins-1):(-colorBins+1)
  breaks = c( (0:4)*m2+1,
              colorBins,
              2*colorBins-(4:0)*m2-1)
  
  Colors.Labels = c( formatC(Lower.max), rep("", 4), "0", rep("",4),  formatC(-Lower.max))
  upper.breaks = sqrt(( (0:4)*m2+1)/(2*3.14*(Upper.max)))
  upper.labels = c( "0", rep("", 3), formatC(Upper.max))
  
  if(0) {
    p = ggplot(data = df, aes(x = x, y = y, size = upper.value, color=value)) 
    p = p+theme_classic()+theme(panel.grid = element_line(color = "gray75"))
    p = p+geom_point(aes(stroke = 0))
    p = p + scale_size_area(name = Upper.name, 
                            limits=c(0,colorBins-1), 
                            breaks=upper.breaks,
                            labels=upper.labels 
                            #,                          max_size = 2
    )
    p = p+scale_color_gradientn(
      colors = rev(Colors),
      #          values = rev(as.numeric(names(Colors))),
      breaks = as.numeric(names(Colors)[breaks]), 
      labels = Colors.Labels,
      name = Lower.name,
      limits = c(-colorBins-1, colorBins-1),
      guide = guide_colorbar())
  } else {
    p = ggplot(data = df, aes(x0 = x, y0 = y, 
                              r = sqrt((upper.value)/(2*3.14*(Upper.max))), fill=value)) 
#    p = p+theme_gray()+theme(axis.ticks = element_blank())
    p = p+theme_minimal()
    p = p + theme(panel.grid = element_line(color = "gray90"),
                  panel.grid.minor = element_blank(), axis.line = element_blank())
    p = p+geom_circle( #size = 0
      linetype = "blank", n = 16,show.legend = TRUE)
    p = p+scale_fill_gradientn(
      colors = rev(Colors),
      #          values = rev(as.numeric(names(Colors))),
      breaks = as.numeric(names(Colors)[breaks]), 
      labels = Colors.Labels,
      name = Lower.name,
      limits = c(-colorBins-1, colorBins-1),
      guide = guide_colorbar())
    p = p + scale_radius(name = Upper.name, 
                        #    limits=c(0,colorBins-1), 
                        #    breaks=upper.breaks,
                       #     labels=upper.labels,
                            guide = guide_legend()
                            #,                          max_size = 2
    )
    if(0)
      p = p + scale_radius(name = Upper.name, 
                           limits=c(0,colorBins-1), 
                           breaks=upper.breaks,
                           labels=upper.labels,
                           guide = guide_legend()
                           #,                          max_size = 2
      )
  }
  
  if( any(Labels.Cols == "ruler") ) {
    x.ruler = which(Labels.Cols[Order.Cols] == "ruler")
    p = p + geom_vline(xintercept = x.ruler, color = "gray50", size = 1)
  }
  if( any(Labels.Rows == "ruler") ) {
    y.ruler = which(Labels.Rows[Order.Rows] == "ruler")
    p = p + geom_hline(yintercept = y.ruler, color = "gray50", size = 1)
  }
  
  x.labels = as.character(Labels.Cols[Order.Cols])
  x.labels = gsub("ruler","", x.labels)
  y.labels = as.character(Labels.Rows[Order.Rows])
  y.labels =  gsub("ruler","", y.labels)
  p = p+scale_x_continuous(breaks = 1:nc, 
                           labels = x.labels,
                           position = "top",
                           expand = c(0,0))
  p = p+scale_y_continuous(breaks = 1:nr, 
                           labels = y.labels,
                           position = "right",
                           expand = c(0,0))
  p = p + labs(x = "", y = "")
 # p = p + theme(panel.grid =  element_blank())
  # p = p + coord_equal()
  p
}

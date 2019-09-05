library(ggplot2)
library(RColorBrewer)

TriangleOffset = 0.3
lower.triangle.x = c(0,1,1, 1-TriangleOffset, 0)
lower.triangle.y = c(0,0,1, 1, TriangleOffset)
upper.triangle.x = c(0,0,1-TriangleOffset)
upper.triangle.y = c(TriangleOffset,1,1)

tile.x = c(0, 0, 1, 1)
tile.y = c(0, 1, 1, 0)

buildPolyGrid = function (nr = 2, nc=3, 
                          Xs = upper.triangle.x, 
                          Ys = upper.triangle.y,
                          prefix = "") {

  Xcoords =  rep(1:nr, nc)
  Ycoords = rep(1:nc, each=nr )
  Ids = factor(c(paste0( prefix, formatC(Xcoords,digits=2,format="d", flag="0"), ".", formatC(Ycoords,,digits=2,format="d", flag="0"))))
 
  l = length(Xs)
  positions = data.frame(
    id = rep(Ids, each = l),
    x = sapply(Xcoords,function(x) Xs+x)[1:(l*nr*nc)]-.5,
    y = sapply(Ycoords,function(y) Ys+y)[1:(l*nr*nc)]-.5
  )
  list( ids = Ids, pos = positions )
}

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
}


plotSplitHeatmap = function( nr, nc, Upper, Lower, 
                             Upper.name = "up", 
                             Lower.name = "down",
                             Upper.max = NA,
                             Lower.max = NA,
                             Order.Rows = NA,
                             Order.Cols = NA,
                             Upper.Pal = brewer.pal(9,"Blues"),
                             Lower.Pal = brewer.pal(9,"YlOrRd"),
                             Lower.Neg.Pal = brewer.pal(9,"Greens"),
                             colorBins = 20,
                             Labels.Rows = NA,
                             Labels.Cols = NA,
                             Sep.Width = 2,
                             Sep.Color = "light gray") {
  
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
  
  Upper.grid = buildPolyGrid(nc,nr,upper.triangle.x,upper.triangle.y, "u")
  Lower.grid = buildPolyGrid(nc,nr,lower.triangle.x,lower.triangle.y, "l")
  tile = buildPolyGrid(nc,nr,tile.x,tile.y,"t")  
  
  # build datasets
  Upper.values = data.frame( id = Upper.grid$ids, 
                             value = processValues(Upper[Order.Rows,Order.Cols],Upper.max), 
                             Group = Upper.name)
  Upper.data = merge(Upper.values,Upper.grid$pos, c("id"))
  
  Lower.values = data.frame( id = Lower.grid$ids, 
                             value = processValues(Lower[Order.Rows,Order.Cols],Lower.max), 
                             Group = Lower.name)
  Lower.data = merge(Lower.values,Lower.grid$pos, c("id"))
  
  df = rbind(Upper.data,Lower.data)
  
  # build color map
  m1 = floor(colorBins/6)+1
  k1 = 6*m1 - colorBins
  m2 = floor((2*colorBins-1)/11)+1
  k2 = 11*m2 - (2*colorBins-1)
  
  Colors = c( rev(colorRampPalette(Upper.Pal, space = "Lab")(colorBins-1)), "white",
              rep("red",k1),
              rev(colorRampPalette(Lower.Pal, space = "Lab")(colorBins-1)), "white",
              colorRampPalette(Lower.Neg.Pal, space = "Lab")(colorBins-1),
              rep("red",k2))
  
  names(Colors) = c( paste0(Upper.name, ".", (colorBins-1):0), 
                     paste0("white",1:k1), 
                     paste0(Lower.name, ".", (colorBins-1):(-colorBins+1)),
                     paste0("black",1:k2))
  breaks = c( (0:4)*m1+1, colorBins, 
              6*m1+(0:4)*m2+1,
              6*m1+colorBins,
              6*m1+2*colorBins-(4:0)*m2-1)
  
  Colors.Labels = c( as.character(Upper.max), rep("", 4), "0",
                     as.character(Lower.max), rep("", 4), "0", rep("",4),  as.character(-Lower.max))

  
  p = ggplot(data = df, aes(x = x, y = y)) 
  p = p+geom_polygon(aes(fill = interaction(Group,value), group = id))
  p = p+geom_polygon(data=tile$pos, aes(group = id), fill = NA, 
                     color=Sep.Color, 
                     size=Sep.Width)
  p = p+theme_minimal()
  
  p = p+scale_fill_manual(values = Colors, 
                          breaks = names(Colors)[breaks], 
                          labels = Colors.Labels,
                          guide = guide_legend(nrow=12, title="-log10(p-val)"))
  
 
  
  p = p+scale_x_continuous(breaks = 1:nc, labels = as.character(Labels.Cols[Order.Cols]),position = "top")
  p = p+scale_y_continuous(breaks = 1:nr, labels = as.character(Labels.Rows[Order.Rows]),position = "right")
  p = p + labs(x = "", y = "")
  p = p + theme(panel.grid =  element_blank())

  p
}

GenerateColorBars = function(  Upper.Pal = brewer.pal(9,"Blues"),
                               Lower.Pal = brewer.pal(9,"YlOrBr"),
                               colorBins = 20 ) {
  tile = buildPolyGrid(nc = colorBins,3,tile.x,tile.y)  
  Colors =c("white", colorRampPalette(Upper.Pal, space = "Lab")(colorBins-1), rep("white", colorBins), "white", colorRampPalette(Lower.Pal, space = "Lab")(colorBins-1))
  
  p = ggplot(tile$pos, aes(x=x,y=y, group = id, color=id, fill = id)) 
  p = p + theme_void()+guides(fill=FALSE, color=FALSE)
  p = p + geom_polygon() 
  p = p + scale_fill_manual(values = Colors ) + scale_color_manual(values=Colors)
}

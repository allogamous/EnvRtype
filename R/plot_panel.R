
plot_panel = function(ECs,title,palette = 'magma',n=100){
  # palette = magma, inferno or virids
  if (!requireNamespace('superheat')) utils::install.packages("superheat")
  if (!requireNamespace('viridis')) utils::install.packages("viridis")
  require(superheat)
  require(viridis)
  if(palette == 'viridis')
  {
    superheat(ECs,
              pretty.order.rows = TRUE,
              pretty.order.cols = TRUE,title = title,
              grid.vline.col = "transparent",
              grid.hline.col = "transparent",
              #row.dendrogram = T,
              legend.width = 4,
              left.label.size = 0.1,
              bottom.label.text.size = 2.5,
              bottom.label.size = 0.3,
              bottom.label.text.angle = 90,
              #heat.pal = viridis::inferno(100),
              heat.pal = viridis::viridis(n),
              legend.text.size = 10,
              #   X.text = round(as.matrix(a),1),X.text.col="white",
              legend.height=0.08)
  }
  if(palette == 'magma')
  {
    superheat(ECs,
              pretty.order.rows = TRUE,
              pretty.order.cols = TRUE,title = title,
              grid.vline.col = "transparent",
              grid.hline.col = "transparent",
              #row.dendrogram = T,
              legend.width = 4,
              left.label.size = 0.1,
              bottom.label.text.size = 2.5,
              bottom.label.size = 0.3,
              bottom.label.text.angle = 90,
              #heat.pal = viridis::inferno(100),
              heat.pal = viridis::magma(n),
              legend.text.size = 10,
              #   X.text = round(as.matrix(a),1),X.text.col="white",
              legend.height=0.08)
  }
  
  if(palette == 'inferno')
  {
    superheat(ECs,
              pretty.order.rows = TRUE,
              pretty.order.cols = TRUE,title = title,
              grid.vline.col = "transparent",
              grid.hline.col = "transparent",
              #row.dendrogram = T,
              legend.width = 4,
              left.label.size = 0.1,
              bottom.label.text.size = 2.5,
              bottom.label.size = 0.3,
              bottom.label.text.angle = 90,
              heat.pal = viridis::inferno(n),
            #  heat.pal = viridis::magma(100)
              legend.text.size = 10,
              #   X.text = round(as.matrix(a),1),X.text.col="white",
              legend.height=0.08)
  }
  
  
  
}

# For getting different fills & plots to align from: http://stackoverflow.com/questions/3805029/different-legends-and-fill-colours-for-facetted-ggplot

align.plots <- function(..., vertical=TRUE){
  #http://ggextra.googlecode.com/svn/trunk/R/align.r
  dots <- list(...)
  dots <- lapply(dots, ggplotGrob)
  ytitles <- lapply(dots, function(.g) editGrob(getGrob(.g,"axis.title.y.text",grep=TRUE), vp=NULL))
  ylabels <- lapply(dots, function(.g) editGrob(getGrob(.g,"axis.text.y.text",grep=TRUE), vp=NULL))
  legends <- lapply(dots, function(.g) if(!is.null(.g$children$legends))
    editGrob(.g$children$legends, vp=NULL) else ggplot2:::.zeroGrob)
  
  gl <- grid.layout(nrow=length(dots))
  vp <- viewport(layout=gl)
  pushViewport(vp)
  widths.left <- mapply(`+`, e1=lapply(ytitles, grobWidth),
                        e2= lapply(ylabels, grobWidth), SIMPLIFY=F)
  widths.right <- lapply(legends, function(g) grobWidth(g) + if(is.zero(g)) unit(0, "lines") else unit(0.5, "lines")) # safe margin recently added to ggplot2
  widths.left.max <- max(do.call(unit.c, widths.left))
  widths.right.max <- max(do.call(unit.c, widths.right))
  
  for(i in seq_along(dots)){
    pushViewport(viewport(layout.pos.row=ii))
    pushViewport(viewport(x=unit(0, "npc") + widths.left.max - widths.left[[ii]],
                          width=unit(1, "npc") - widths.left.max + widths.left[[ii]] -
                            widths.right.max + widths.right[[ii]],
                          just="left"))
    grid.draw(dots[[ii]])
    upViewport(2)
  }
}



p <- ggplot(datapoly[datapoly$variable=="val1",], aes(x=x, y=y)) + geom_polygon(aes(fill=value, group=id),colour="black")
p1 <- ggplot(datapoly[datapoly$variable=="val2",], aes(x=x, y=y)) + geom_polygon(aes(fill=value, group=id),colour="black")
align.plots( p,p1)
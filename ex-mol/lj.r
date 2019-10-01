lj <- function (dist, scale, width, maxpe)
{ 
  pmin ( 4*scale*((width/dist)^12 - (width/dist)^6),  
         pmax (0, maxpe - 4 * scale * (dist/width)^2 ) )
}

plot.lj <- function (d, scale, width, maxpe)
{
  plot (d, lj(d,scale,width,maxpe), xlab="", ylab="", type="l", lwd=3)
  abline(h=0,v=(-10):10)
}

# DATA-PLOT.R - Functions for plotting measurement data.


# DEFAULT SET OF LEVELS TO USE AS BREAKPOINTS.

def.levels = round (exp(seq(-2,3,by=0.1)), 2)


# MAKE CONTOUR PLOT OF DATA ON A GRID.  The first argument is the
# matrix of measurement locations, of which only the first two columns
# are used.  These must contain a grid of (x,y) locations, in increasing
# order, with x varying more slowly.  The second argument is the vector 
# of measurements at these locations.  Any additional arguments are passed
# on to the "contour" function.

data.contour = function (loc, data, xlab="x", ylab="y", levels=def.levels, ...)
{
  x.grid = unique(sort(loc[,1]))
  y.grid = unique(sort(loc[,2]))

  contour (x.grid, y.grid, matrix(data,length(x.grid),length(y.grid),byrow=T),
           levels=levels, xlab=xlab, ylab=ylab, ...)
}


# MAKE VARYING-SIZE DOT PLOT OF DATA ON A GRID.  The first argument is the
# matrix of measurement locations, of which only the first two columns
# (for x and y) are used.  The second argument is the vector of measurements 
# at these locations.  The third, "levels", argument gives the breakpoints for 
# determining dot shade.  The fourth, "add", argument is F (default) if a 
# new plot should be started.  Any additional arguments are passed on to the 
# "plot" function.

data.dots = function (loc, data, add=F, xlab="x", ylab="y", levels=def.levels,
                      ...)
{
  if (!add)
  { plot(range(loc[,1]),range(loc[,2]),xlab=xlab,ylab=ylab,type="n")
  }

  g = round (100 * apply (outer(data,levels,"<"), 1, sum) / length(levels))

  points (loc[,1], loc[,2], pch=20, cex=4, col=paste("gray",g,sep=""))
}


# MAKE A SET OF PLOTS, ONE FOR EACH HEIGHT.  The first argument is the matrix 
# of measurement locations, of which only the first three columns are used.  
# These must contain  (x,y,z) locations.  When only values for a particular 
# value of z (height) are looked at, the (x,y) values must form a grid, with 
# values in increasing, x varying more slowly than y.  The second argument is 
# the vector of measurements at these locations.  The third argument is a
# list of plotting functions (or a single function, just "data.contour" by
# default) which are used create each plot.  Any additional arguments are 
# passed on to the plotting functions.
#
# Three plots are put on each page.

data.multi = function (loc, data, plotfns="data.contour", ...)
{
  par(mfrow=c(3,1))

  if (!is.list(plotfns))
  { plotfns = list(plotfns)
  }

  z.grid = unique(sort(loc[,3],decreasing=T))

  for (z in z.grid)
  { w = loc[,3]==z
    for (i in 1:length(plotfns))
    { if (i==1)
      { plotfns[[i]] (loc[w,], data[w], ...)
      }
      else
      { plotfns[[i]] (loc[w,], data[w], add=T, ...)
      }
    }
    title(paste("z =",z))
  }
}

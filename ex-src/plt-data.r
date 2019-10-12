
# Plot data generated on the four grids.  The path to the ex-src directory
# may need to be modified (see the line below).

dir = "c:/cygwin/home/radford/dnd/fbm/ex-src"

postscript(paste(dir,"data-plots.ps",sep="/"),horiz=F)

source(paste(dir,"data-plot.r",sep="/"))

g1 = as.matrix(read.table(paste(dir,"grid1",sep="/")))
g2 = as.matrix(read.table(paste(dir,"grid2",sep="/")))
g3 = as.matrix(read.table(paste(dir,"grid3",sep="/")))
g4 = as.matrix(read.table(paste(dir,"grid4",sep="/")))

d1t = scan(paste(dir,"data-grid1-0-1",sep="/"))
d1 = scan(paste(dir,"data-grid1-0.1-1",sep="/"))
d2 = scan(paste(dir,"data-grid2-0.1-1",sep="/"))
d3 = scan(paste(dir,"data-grid3-0.1-1",sep="/"))
d4 = scan(paste(dir,"data-grid4-0.1-1",sep="/"))

data.multi (g1, d1t, list(data.contour))
data.multi (g1, d1t, list(data.dots))
data.multi (g1, d1t, list(data.dots, data.contour))

data.multi (g1, d1, list(data.dots))
data.multi (g2, d2, list(data.dots))
data.multi (g3, d3, list(data.dots))
data.multi (g4, d4, list(data.dots))

dev.off()

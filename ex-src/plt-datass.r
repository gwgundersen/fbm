
# Plot data generated from the start-stop model.  The path to the ex-src 
# directory may need to be modified (see the line below).

dir = "c:/cygwin/home/radford/dnd/fbm/ex-src"

postscript(paste(dir,"data-plotsss.ps",sep="/"),horiz=F)

source(paste(dir,"data-plot.r",sep="/"))

g1ss = as.matrix(read.table(paste(dir,"gridss1",sep="/")))
g2ss = as.matrix(read.table(paste(dir,"gridss2",sep="/")))

d1sst = scan(paste(dir,"data-gridss1-0-1",sep="/"))
d1ss = scan(paste(dir,"data-gridss1-0.1-1",sep="/"))
d2ss = scan(paste(dir,"data-gridss2-0.1-1",sep="/"))

for (t in unique(g1ss[,4]))
{ s = g1ss[,4]==t
  data.multi (g1ss[s,], d1sst[s], list(data.dots))
  mtext(paste("t =",t),side=4,line=0.3)
}

for (t in unique(g2ss[,4]))
{ s = g2ss[,4]==t
  data.multi (g2ss[s,], d2ss[s], list(data.dots))
  mtext(paste("t =",t),side=4,line=0.3)
}

dev.off()

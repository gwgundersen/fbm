
# Plot predictions found from fitting data on the four grids, with the
# number of sources fixed at 3 or variable.  The path to the ex-src
# directory may need to be modified (see the line below).

dir = "c:/cygwin/home/radford/dnd/fbm/ex-src"

postscript(paste(dir,"pred-plots.ps",sep="/"),horiz=F)

source(paste(dir,"data-plot.r",sep="/"))


pc4 = as.matrix(read.table(paste(dir,"predc4",sep="/")))

data.multi (pc4[,1:3], pc4[,5], list(data.contour))
data.multi (pc4[,1:3], pc4[,5], list(data.dots))
data.multi (pc4[,1:3], pc4[,5], list(data.dots, data.contour))

pc3 = as.matrix(read.table(paste(dir,"predc3",sep="/")))

data.multi (pc3[,1:3], pc3[,5], list(data.contour))
data.multi (pc3[,1:3], pc3[,5], list(data.dots))
data.multi (pc3[,1:3], pc3[,5], list(data.dots, data.contour))

pc2 = as.matrix(read.table(paste(dir,"predc2",sep="/")))

data.multi (pc2[,1:3], pc2[,5], list(data.contour))
data.multi (pc2[,1:3], pc2[,5], list(data.dots))
data.multi (pc2[,1:3], pc2[,5], list(data.dots, data.contour))


pv4 = as.matrix(read.table(paste(dir,"predv4",sep="/")))

data.multi (pv4[,1:3], pv4[,5], list(data.contour))
data.multi (pv4[,1:3], pv4[,5], list(data.dots))
data.multi (pv4[,1:3], pv4[,5], list(data.dots, data.contour))

pv3 = as.matrix(read.table(paste(dir,"predv3",sep="/")))

data.multi (pv3[,1:3], pv3[,5], list(data.contour))
data.multi (pv3[,1:3], pv3[,5], list(data.dots))
data.multi (pv3[,1:3], pv3[,5], list(data.dots, data.contour))

pv2 = as.matrix(read.table(paste(dir,"predv2",sep="/")))

data.multi (pv2[,1:3], pv2[,5], list(data.contour))
data.multi (pv2[,1:3], pv2[,5], list(data.dots))
data.multi (pv2[,1:3], pv2[,5], list(data.dots, data.contour))


dev.off()

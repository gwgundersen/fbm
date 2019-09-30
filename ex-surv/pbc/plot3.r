ps.options(width=6,height=7.8,horizontal=F,pointsize=10)
postscript("plot3.ps")

par(mfrow=c(3,2))
par(mar=c(4,3,1,1)+0.1)
par(mgp=c(1.5,0.5,0))
par(cex.lab=1.2)

run <- "nnn-12a"
iters <- seq(150,400,50)

cov.effect("ascites",run,iters[1],pch=".",lwd=2,
 xlab="median survival without ascites",
 ylab="median survival with ascites")

cov.effect("ascites",run,iters[2],pch=".",lwd=2,
 xlab="median survival without ascites",
 ylab="median survival with ascites")

cov.effect("ascites",run,iters[3],pch=".",lwd=2,
 xlab="median survival without ascites",
 ylab="median survival with ascites")

cov.effect("ascites",run,iters[4],pch=".",lwd=2,
 xlab="median survival without ascites",
 ylab="median survival with ascites")

cov.effect("ascites",run,iters[5],pch=".",lwd=2,
 xlab="median survival without ascites",
 ylab="median survival with ascites")

cov.effect("ascites",run,iters[6],pch=".",lwd=2,
 xlab="median survival without ascites",
 ylab="median survival with ascites")

dev.off()

ps.options(width=6,height=7.8,horizontal=F,pointsize=10)
postscript("plot2.ps")

run <- "nnn-12a"
iters <- seq(160,400,10)

par(mfrow=c(3,2))
par(mar=c(4,3,1,1)+0.1)
par(mgp=c(1.5,0.5,0))
par(cex.lab=1.2)

cov.effect("edema",   run,iters,pch=".",lwd=2,
 xlab="median survival without edema",
 ylab="median survival with edema")
cov.effect("bili",    run,iters,pch=".",lwd=2,
 xlab="median survival with low bilirubin",
 ylab="median survival with high bilirubin")
cov.effect("albumin",    run,iters,pch=".",lwd=2,
 xlab="median survival with low albumin",
 ylab="median survival with high albumin")
cov.effect("alkphos",run,iters,pch=".",lwd=2,
 xlab="median survival with low alkphos",
 ylab="median survival with high alkphos")
cov.effect("platelet",run,iters,pch=".",lwd=2,
 xlab="median survival with low platelet count",
 ylab="median survival with high platelet count")
cov.effect("protime",run,iters,pch=".",lwd=2,
 xlab="median survival with low prothrombin time",
 ylab="median survival with high prothrombin time")

dev.off()

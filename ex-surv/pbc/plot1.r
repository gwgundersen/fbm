ps.options(width=6,height=7.8,horizontal=F,pointsize=10)
postscript("plot1.ps")

run <- "nnn-12a"
iters <- seq(160,400,10)

par(mfrow=c(3,2))
par(mar=c(4,3,1,1)+0.1)
par(mgp=c(1.5,0.5,0))
par(cex.lab=1.2)

cov.effect("drug",   run,iters,pch=".",lwd=2,
 xlab="median survival taking placebo",
 ylab="median survival taking drug")
cov.effect("age",    run,iters,pch=".",lwd=2)
cov.effect("sex",    run,iters,pch=".",lwd=2,
 xlab="median survival for males",
 ylab="median survival for females")
cov.effect("ascites",run,iters,pch=".",lwd=2,
 xlab="median survival without ascites",
 ylab="median survival with ascites")
cov.effect("hepatom",run,iters,pch=".",lwd=2,
 xlab="median survival without hepatomegaly",
 ylab="median survival with hepatomegaly")
cov.effect("spiders",run,iters,pch=".",lwd=2,
 xlab="median survival without spiders",
 ylab="median survival with spiders")

dev.off()

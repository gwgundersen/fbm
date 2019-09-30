ps.options(width=6,height=7.8,horizontal=F,pointsize=10)
postscript("plot4.ps")

par(mfrow=c(3,2))
par(mar=c(4,3,1,1)+0.1)
par(mgp=c(1.5,0.5,0))
par(cex.lab=1.2)

cov.effect("bili","nlx-12b",200,pch=".",lwd=2,
 xlab="median survival with low bilirubin",
 ylab="median survival with high bilirubin")

cov.effect("bili","nlx-12b",seq(80,200,5),pch=".",lwd=2,
 xlab="median survival with low bilirubin",
 ylab="median survival with high bilirubin")

cov.effect("bili","nnx-12a",100,pch=".",lwd=2,
 xlab="median survival with low bilirubin",
 ylab="median survival with high bilirubin")

cov.effect("bili","nnx-12a",seq(28,100,3),pch=".",lwd=2,
 xlab="median survival with low bilirubin",
 ylab="median survival with high bilirubin")

cov.effect("bili","xxn-12b",300,pch=".",lwd=2,
 xlab="median survival with low bilirubin",
 ylab="median survival with high bilirubin")

cov.effect("bili","xxn-12b",seq(60,300,10),pch=".",lwd=2,
 xlab="median survival with low bilirubin",
 ylab="median survival with high bilirubin")

dev.off()

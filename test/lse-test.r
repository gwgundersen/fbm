f <- function (x)
{ log (3*dnorm(x,2,1) + 2*dnorm(x,-1,1.5) +
       dnorm(x,-2.2,0.05) + dnorm(x,-1.8,0.04) + dnorm(x,-2.0,0.07) +
       dnorm(x,1.7,0.08) + dnorm(x,1.9,0.06) + dnorm(x,2.2,0.07) )
}

x.trn <- round(seq(-3,3,0.03),2)
x.tst <- round(seq(-3.015,2.991,0.03),2)

y.trn <- round (f(x.trn) + rnorm(length(x.trn),0,0.1), 4)
y.tst <- round (f(x.tst) + rnorm(length(x.tst),0,0.1), 4)

y.trnt <- round (f(x.trn), 4)
y.tstt <- round (f(x.tst), 4)

write.table(cbind(x.trn,y.trn),"lse-trn",quote=F,row.names=F,col.names=F)
write.table(cbind(x.tst,y.tst),"lse-tst",quote=F,row.names=F,col.names=F)
write.table(cbind(x.trn,y.trnt),"lse-trnt",quote=F,row.names=F,col.names=F)
write.table(cbind(x.trn,y.tstt),"lse-tstt",quote=F,row.names=F,col.names=F)

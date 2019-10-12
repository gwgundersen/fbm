set.seed(1)
f <- 10.005
p <- read.table("pr",head=F)[,1]
n <- ceiling (f * 4*p*(1-p))
r <- rep (1:length(p), times=n)
i <- sample(r)
write.table(i,"ix",quote=F,row.names=F,col.names=F)

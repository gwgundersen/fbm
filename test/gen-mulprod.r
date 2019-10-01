set.seed(2)
x <- round (runif(300,-1,1), 5)
y <- round (cos(x+0.3) + 0.1 * rnorm(300) + sqrt(0.3^2*x^2) * rnorm(300), 5)
write.table(cbind(x,y),"mulprod.dat",quote=F,row.names=F,col.names=F)

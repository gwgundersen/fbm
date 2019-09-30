# Generate data for the linear regression model used as an example in the 
# technical report on "Annealed importance sampling" by R. M. Neal.

set.seed(1)

K <- 10
N <- 100

rho <- 0.9
sigma <- 1

b <- c(1,0.5,-0.5,0,0,0,0,0,0,0)

C <- matrix(rho,K,K)
diag(C) <- 1
U <- chol(C)

X <- matrix(rnorm(K*N),N,K) %*% U
y <- X %*% b + rnorm(N,0,sigma^2)

X.test <- matrix(rnorm(K*10000),10000,K) %*% U
y.test <- X.test %*% b + rnorm(10000,0,sigma^2)

lmod <- lm(y~X-1)

write.table(round(cbind(X,y),6),"sdata"," ")
write.table(round(cbind(X.test,y.test),6),"sdata.test"," ")
write.table(diag(K),"sdiag"," ")

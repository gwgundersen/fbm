# Generate data for the random effects example.

set.seed(2)

u <- 25
v1 <- 12^2
v2 <- 6^2
K <- 18

i <- c(1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,4)

t0 <- rep(0,K)
t1 <- rep(0,K)

for (j in 1:K)
{ x <- rnorm (i[j], rnorm(1,u,sqrt(v1)), sqrt(v2))
  t0[j] <- mean(x)
  if (i[j]==1) t1[j] <- 0 else t1[j] <- var(x)
}

t0 <- round(t0,3)
t1 <- round(t1,4)

write.table(cbind(i,t0,t1),"edata",dimnames.write=F,sep=" ")

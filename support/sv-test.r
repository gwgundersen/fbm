# Test different approaches to implementing sample-values for gp module.


# SET UP TEST.

set.seed(1)

n <- 5				# Number of training cases

X <- matrix(rnorm(n*n),n,n)		# Covariance of latent values
V <- X %*% t(X)

C <- V + diag(rep(0.1,n))		# Covariance of targets

t <- as.vector(chol(C) %*% rnorm(n))	# Target values


# METHOD 1

method1 <- function()
{
  U <- chol(C)
  Ci <- chol2inv(U)
  VCi <- V %*% Ci
  sm <- as.vector (VCi %*% t)
  sU <- chol(V-VCi%*%V)
  list(m=sm,U=sU)
}

cat ("method1:",system.time (m1 <- method1()),"\n")


# METHOD 2

method2 <- function()
{ 
  U <- chol(C)
  LiV <- forwardsolve(t(U),V)
  sm <- as.vector (forwardsolve(t(U),t) %*% LiV)
  sU <- chol(V-t(LiV) %*% LiV)
  list(m=sm,U=sU)
}

cat ("method2:",system.time (m2 <- method2()),"\n")


# METHOD 3

method3 <- function()
{ 
  U <- chol(C)
  sm <- as.vector (V %*% backsolve(U,forwardsolve(t(U),t)))
   
  Vi <- chol2inv(chol(V))
  diag(Vi) <- diag(Vi) + 10
  sU <- chol(chol2inv(chol(Vi)))

  list(m=sm,U=sU)
}

cat ("method3:",system.time (m3 <- method3()),"\n")


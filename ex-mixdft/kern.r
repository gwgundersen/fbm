# KERN.R - Kernel density estimation. */

kern.pr <- function (x, y, h)
{
  k <- ncol(x)
  h <- rep(h,length=k)

  n <- nrow(x)
  m <- nrow(y)

  p <- rep(0,m)

  for (j in 1:m)
  { for (i in 1:n)
    { p[j] <- p[j] + exp( - 0.5*sum(((x[i,]-y[j,])/h)^2))
    }
  }

  p <- p / sqrt(prod(2*pi*h^2))
  p <- p / n

  log(p)
}

kern.sample <- function (x, m, h)
{
  k <- ncol(x)
  h <- rep(h,length=k)

  n <- nrow(x)

  y <- matrix(0,m,k)

  for (j in 1:m)
  { i <- sample(n,1)
    y[j,] <- rnorm(k,x[i,],h)
  }

  y
}

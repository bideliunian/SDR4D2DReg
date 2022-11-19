
## l2 distance between two beta distribution
dist4beta_l2 <- function(a1, b1, a2, b2, ndSup=101){
  # Param
  #   dist1: Beta(a1, b1)
  #   dist2: Beta(a2, b2)
  #   ndSup: number of grids on support
  # Return 
  #   l2 distance between two beta distribution
  dSup <- seq(0,1,length.out = ndSup)[2:(ndSup-1)]
  d1 <- dbeta(dSup,a1,b1)
  d2 <- dbeta(dSup,a2,b2)
  return(sqrt(pracma::trapz(dSup, (d1 - d2)^2)))
}

## gram matrix with l2 distance gaussian kernel
gram_beta_l2 <- function(a1, b1, a2, b2, complexity){
  # Param
  #   a1, b1, a2, b2, length n vector 
  #   dist1: Beta(a1, b1)
  #   dist2: Beta(a2, b2)
  #   complexity controls the tuning parameter in kernel
  # Return 
  #   n by m gram matrix
  n <- length(a1)
  m <- length(a2)
  a1 <- as.vector(a1)
  b1 <- as.vector(b1)
  a2 <- as.vector(a2)
  b2 <- as.vector(b2)
  k <- outer(
    1:n, 1:m,
    Vectorize(function(i,k) dist4beta_l2(a1[i],b1[i],a2[k],b2[k]))
  )
  sigma <- sum(k^2)/choose(n,2)
  gamma <- complexity/(sigma)
  return(exp(-gamma*(k^2)))
}

## wasserstein distance between two beta distribution
dist4beta_w <- function(a1,b1,a2,b2,nqSup=101){
  d1 <- list()
  d2 <- list()
  qSup <- seq(0,1,length.out = nqSup)
  d1$y <- qbeta(qSup,a1,b1)
  d1$x <- qSup
  d2$y <- qbeta(qSup,a2,b2)
  d2$x <- qSup
  return(dist4den(d1, d2, fctn_type = 'quantile'))
}

## gram matrix with wasserstein2 distance gaussian kernel
gram_beta_w <- function(a1,b1,a2,b2,complexity){
  n <- length(a1)
  m <- length(a2)
  a1 <- as.vector(a1)
  b1 <- as.vector(b1)
  a2 <- as.vector(a2)
  b2 <- as.vector(b2)
  k <- outer(
    1:n, 1:m,
    Vectorize(function(i,k) dist4beta_w(a1[i], b1[i], a2[k], b2[k]))
  )
  sigma <- sum(k^2)/choose(n,2)
  gamma <- complexity/(sigma)
  return(exp(-gamma*(k^2)))
}


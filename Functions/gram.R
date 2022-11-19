##########################################################################
######## function to calculate the gaussian gram matrix of distributions
#########################################################################

## function calculate the gram matrix of euclidean random vectors

gram_gauss_eucl <- function(x, x_new, complexity){
  # Param: 
  #   x n by p matrix
  #   x_new m by p matrix
  #   complexity controls the tuning parameter in kernel
  # Return
  #   n by m gram matrix
  x <- as.matrix(x)
  x_new <- as.matrix(x_new)
  n <- dim(x)[1]
  m <- dim(x_new)[1]
  k2 <- x%*%t(x)
  k1 <- t(matrix(diag(k2),n,n))
  k3 <- t(k1)
  k <- k1 - 2*k2 + k3
  sigma <- sum(sqrt(k)) / (2*choose(n,2))
  gamma <- complexity / (2*sigma^2)
  k.new.1 <- matrix(diag(x%*%t(x)), n, m)
  k.new.2 <- x%*%t(x_new)
  k.new.3 <- matrix(diag(x_new%*%t(x_new)), m, n)
  return(exp(-gamma*(k.new.1 - 2*k.new.2 + t(k.new.3))))
}

## function for warsserstein gram matrix between two distributions
gram_gauss_w <- function(x, x_new, complexity){
  # Param: 
  #   x n by p matrix, each row is one observation of univariate distribution
  #   x_new m by p matrix, each row is one observation of univariate distribution
  #   complexity controls the tuning parameter in kernel
  # Return
  #   n by m gram matrix
  x <- as.matrix(x)
  x_new <- as.matrix(x_new)
  n <- dim(x)[1]
  m <- dim(x_new)[1]
  k <- outer( 
    1:n, 1:m,
    Vectorize(function(i,k) wasserstein1d(x[i,],x_new[k,]))
  )
  sigma2 <- sum(k^2)/choose(n,2)
  gamma <- complexity/(sigma2)
  return(exp(-gamma*k^2))
}

## function for sliced warsserstein gram matrix between two multi-dim distributions
gram_gauss_sw <- function(x, x_new, complexity, n_proj = 50, n_seed=1){
  # Param: 
  #   x is n by p by d array, each row is one observation of d-dim distribution
  #   x_new is m by p by d array, each row is one observation of d-dim distribution
  #   complexity controls the tuning parameter in kernel
  #   n_proj is the number of projections when computing sliced wass distance
  #   n_seed is the random seed
  # Return
  #   n by m gram matrix
  n <- dim(x)[1]
  m <- dim(x_new)[1]
  k <- outer( 
    1:n, 1:m,
    Vectorize(function(i,k) sliced_wass_distance(x[i,,],x_new[k,,],
                                                 n_projections = n_proj,
                                                 n_seed = n_seed))
  )
  sigma2 <- sum(k^2) / choose(n,2)
  gamma <- complexity / (sigma2)
  return(exp(-gamma*k^2))
}

## function for l2 gram matrix between two distributions
gram_gauss <- function(x, x_new, complexity, metric){
  # Param: 
  #   x is n by p by d array, each row is one observation of a 1-dim distribution
  #   x_new is m by p by d array, each row is one observation of a 1-dim distribution
  #   complexity controls the tuning parameter in kernel
  #   n_seed is the random seed
  # Return
  #   n by m gram matrix
  n <- dim(x)[1]
  m <- dim(x_new)[1]
  if (metric == 'l2'){
    k <- outer( 
      1:n, 1:m,
      Vectorize(function(i,k) l2density(x[i,],x_new[k,]))
    )  
  }
  else if (metric == 'l1'){
    k <- outer( 
      1:n, 1:m,
      Vectorize(function(i,k) l1density(x[i,],x_new[k,]))
    )
  }
  sigma2 <- sum(k^2) / choose(n,2) / 2
  gamma <- complexity / (sigma2)
  return(exp(-gamma*k^2))
}


## sliced warsserstein kernel gram matrix between two multi-dim distributions(faster)
gram_sw <- function(x, complexity, n_proj = 50, n_seed=1){
  # Param
  #   x is n by m by d array, each row is one observation of a d-dim distribution
  #   complexity controls the tuning parameter in kernel
  #   n_proj is the number of projections when computing sliced wass distance
  #   n_seed is the random seed
  # Return
  #   n by n gram matrix
  n <- dim(x)[1]
  ##upper triangular index
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  dist <- function(i,k){
    return(sliced_wass_distance(x[i,,], x[k,,], n_projections = n_proj, n_seed = n_seed))
  }
  kupper <- mapply(dist, i=ind[,1], k=ind[,2])
  k <- diag(n)
  k[upper.tri(k,diag = FALSE)] <- kupper^2
  k <- k+t(k)-2*diag(n)
  sigma2 <- sum(kupper^2)/choose(n,2)
  gamma <- complexity/(2*sigma2)
  return(exp(-gamma*k))
}

## warsserstein kernel gram matrix between two 1-d distributions(faster)
gram_w <- function(x, complexity){
  # Param
  #   x is n by m matrix, each row is one observation of a d-dim distribution
  #   complexity controls the tuning parameter in kernel
  #   n_proj is the number of projections when computing sliced wass distance
  #   n_seed is the random seed
  # Return
  #   n by n gram matrix
  x <- as.matrix(x)
  n <- dim(x)[1]
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  dist <- function(i,k){
    return(wasserstein1d(x[i,],x[k,]))
  }
  kupper <- mapply(dist, i=ind[,1], k=ind[,2])
  k <- diag(n)
  k[upper.tri(k,diag <- FALSE)]<-kupper^2
  k <- k+t(k)-2*diag(n)
  sigma2 <- sum(kupper^2)/choose(n,2)
  gamma <- complexity/(2*sigma2)
  return(exp(-gamma*k))
}

## l2 kernel gram matrix between two 1-d distributions(faster)
gram <- function(x, complexity, metric){
  # Param
  #   x is n by p matrix, each row is one observation of 1-dim distribution
  #   complexity controls the tuning parameter in kernel
  # Return
  #   n by m gram matrix
  x <- as.matrix(x)
  n <- dim(x)[1]
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  
  if (metric == 'l2'){
    dist <- function(i, k){
      return(l2density(x[i,],x[k,]))
    }  
  }
  else if(metric == 'l1'){
    dist <- function(i, k){
      return(l1density(x[i,],x[k,]))
    }
  }
  
  kupper <- mapply(dist, i=ind[,1], k=ind[,2])
  k <- diag(n)
  k[upper.tri(k,diag = FALSE)] <- kupper^2
  k <- k+t(k)-2*diag(n)
  sigma2 <- sum(kupper^2)/choose(n,2)
  gamma <- complexity/(2*sigma2)
  return(exp(-gamma*k))
}
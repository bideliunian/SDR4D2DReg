
##################
## distance cor 
## input:
##      X: n by n distance matrix
##      Y: n by n distance matrix
## output: 
##      distance correaltion
##############################################
dcor <- function(X, Y, centered=TRUE){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  m <- nrow(Y)
  if (n != m) stop("Sample sizes must agree")
  
  ## input distance matrix of Y
  if (!centered){X <- center(X); Y <- center(Y)}
  
  dCov <- sqrt(mean(X * Y))
  dVarX <- sqrt(mean(X * X))
  dVarY <- sqrt(mean(Y * Y))
  V <- sqrt(dVarX * dVarY)
  if (V > 0)
    dCor <- dCov / V else dCor <- 0
  return(dCor)
}


## centering a distance matrix
center <- function(x) {
  d <- as.matrix(x)
  m <- rowMeans(d)
  M <- mean(d)
  a <- sweep(d, 1, m)
  b <- sweep(a, 2, m)
  return(b + M)
}

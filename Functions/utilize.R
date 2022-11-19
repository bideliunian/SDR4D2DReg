##############################
## trace
tr <- function(x) return(sum(diag(x)))## trace of matrix

## operator norm
onorm <- function(a) return(eigen(round((a+t(a))/2,8))$values[1])

## symmetry a asymmetric matrix
sym <- function(a) return(round((a+t(a))/2,9))

## matrix power
matpower <- function(a,alpha){
  a <- (a+t(a))/2;
  tmp <- eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))
}

# trace
trace <- function(x){
  x <- as.matrix(x)
  return(sum(diag(x)))
}

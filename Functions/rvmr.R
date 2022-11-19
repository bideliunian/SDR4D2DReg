##################################
## function to calculate RVMR ###
####################################

sign <- function(x){
  if(sum(abs(x))==0) return(x)
  else return(x/sqrt(sum(x^2)))
}

rvmr <- function(x, y){
  # Param
  #   x: n by r matrix
  #   y: m by s matrix
  # Return
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  n <- dim(x)[1]
  m <- dim(y)[1]
  r <- dim(x)[2]
  s <- dim(y)[2]
  
  tx <- matrix(0, nrow = n, ncol = r)
  ty <- matrix(0, nrow = m, ncol = s)
  
  fx <- function(i,j) {sign(x[i,]-x[j,])}
  fy <- function(i,j) {sign(y[i,]-y[j,])}
  
  # compute multivariate rank of x
  for(k in 1:n){
    if(r==1) {
      tx[k,] <- sum(mapply(i=c(1:n)[-k], fx, j=k)) / n
    }
    else {
      tx[k,] <- rowSums(mapply(i=c(1:n)[-k], fx, j=k)) / n 
    }
  }
  # compute multivariate rank of y
  for (k in 1:m) {
    if(s==1) {
      ty[k,] <- sum(mapply(i=c(1:m)[-k], fy, j=k)) / m
    }
    else {
      ty[k,] <- rowSums(mapply(i=c(1:m)[-k], fy, j=k)) / m 
    }
  }
  
  rv <- tr(cov(tx,ty)%*%cov(ty,tx)) / 
    sqrt(tr(var(tx)%*%var(tx)) * tr(var(ty)%*%var(ty)))
  
  return(rv)
}




########################################
## generalized CV to tune epsilon_x and epsilon_y
######################################
gcv <- function(x, y, eps, which, ytype, xtype, complex_x, complex_y){
  n <- dim(x)[1]
  if(xtype=="scalar") Kx <- gram_gauss_eucl(x,x,complex_x)
  else if(xtype=='distribution') Kx <- gram_w(x,complex_x)
  else if(xtype=='mult_distribution') Kx <- gram_sw(x,complex_x)
  else if(xtype=='gram') Kx <- x
  
  if(ytype=="scalar") Ky <- gram_gauss_eucl(y,y,complex_y)
  else if(ytype=="distribution") Ky <- gram_w(y,complex_y)
  else if(ytype=="mult_distribution") Ky <- gram_sw(y,complex_y)
  else if(ytype=='gram') Ky <- y
  
  if(which=="ey") {
    G1 <- Kx
    G2 <- Ky
  }
  else if(which=="ex") {
    G1 <- Ky
    G2 <- Kx
  }
  
  G2inv <- solve(G2 + eps*onorm(G2)*diag(n))
  nu <- sum((G1 - G2%*%G2inv%*%G1)^2)
  de <- (1 - sum(diag(G2inv%*%G2))/n)^2
  
  return(nu / de)
}


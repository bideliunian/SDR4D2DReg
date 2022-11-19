###############################
## wasserstein generalized sir
###############################
wgsir <- function(x, x_new, y, xtype, ytype, atype, ex, ey, complex_x, complex_y, r){
  ## Param
  ##  x: the n*m matrix of the predictors; 
  ##  x_new: can be set to x or the predictor matrix for the testing set; 
  ##  y: the response object;
  ##  ytype: string{"distribution", "scalar"}, indicating the type of response;
  ##  atype: string "identity" or "Gyinv", indicating the type of matrix A;
  ##  ex,ey: control the regularization magnitude;
  ##  complex_x and complex_y: measure the complexity of HX and HY ;
  ##  r: is the number of leading GSIR predictors in the output.
  ## Return:
  ##  list of estimated sufficient predictors and GSIR candidate matrix
  if(r == 0){
    warning("The dimension is 0", call. = FALSE)
  }
  
  n <- dim(x)[1]
  Q <- diag(n)-rep(1,n)%*%t(rep(1,n))/n
  if(xtype=="scalar") Kx <- gram_gauss_eucl(x, x, complex_x)
  else if(xtype=='distribution') Kx <- gram_w(x, complex_x)
  else if(xtype=='mult_distribution') Kx <- gram_sw(x, complex_x)
  else if(xtype=='gram') Kx <- x
  
  if(ytype=="scalar") Ky <- gram_gauss_eucl(y, y, complex_y)
  else if(ytype=='distribution') Ky <- gram_w(y, complex_y)
  else if(ytype=='mult_distribution') Ky <- gram_sw(y, complex_y)
  else if(ytype=='gram') Ky <- y
  
  Gx <- Q%*%Kx%*%Q
  Gy <- Q%*%Ky%*%Q
  Gxinv <- matpower(sym(Gx+ex*(eigen(Gx)$values[1])*diag(n)),-1)
  Gyinv <- matpower(sym(Gy+ey*onorm(Gy)*diag(n)),-1)
  a1 <- Gxinv%*%Gx
  
  if (!(atype %in% c('identity', 'Gyinv'))) {
    stop("A matrix type must be one of 'identity' 'Gyinv'")
  }
  if(atype=="identity") a2 <- Gy
  else a2 <- Gy%*%Gyinv
  
  gsir <- a1%*%a2%*%t(a1)
  gsir <- sym(gsir)
  v <- eigen(gsir)$vectors[,1:r]
  
  if(xtype=="scalar") Kx_new <- gram_gauss_eucl(x,x_new,complex_x)
  else if(xtype=='distribution') Kx_new <- gram_gauss_w(x,x_new,complex_x)
  else if(xtype=='mult_distribution') Kx_new <- gram_gauss_sw(x,x_new,complex_x)
  else if(xtype=='gram') Kx_new <- x_new
  
  pred.est <- t(t(v)%*%Gxinv%*%Q%*%Kx_new)
  
  return(list(pred.est=pred.est, mat=gsir))
}

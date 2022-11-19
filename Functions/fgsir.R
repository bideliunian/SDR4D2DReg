## fgsir fucntional generalized sir

#################################
## fgsir
####################################
fgsir <- function(x,x_new,y,xtype,ytype,atype,ex,ey,complex_x,complex_y,r,
                  Kx=NULL,Ky=NULL,Kx_new=NULL){
  # Param
  #   x is the n*m matrix of the predictors; 
  #   x_new can be set to x or the predictor matrix for the testing set; y is the response variable;
  #   ytype is a variable indicating the type of response: it can be "distribution" or "scalar";
  #   atype is a variable indicating the type of matrix A: it can be "identity" or "Gyinv";
  #   ex,ey control the regularization magnitude;
  #   complex_x and complex_y measure the complexity of HX and HY ;
  #   r is the number of leading GSIR predictors in the output.
  #   function for l2 gram matrix between two densities
  n <- dim(x)[1]
  Q <- diag(n)-rep(1,n)%*%t(rep(1,n))/n
  
  if(is.null(Kx)){
    if(xtype=="scalar") Kx <- gram_gauss_eucl(x,x,complex_x)
    else if(xtype=='distribution') Kx <- gram_gauss_l2(x,x,complex_x)
  }
  
  if(is.null(Ky)){
    if(ytype=="scalar") Ky <- gram_gauss_eucl(y,y,complex_y)
    else if(ytype=="distribution") Ky <- gram_gauss_l2(y,y,complex_y)
  }
  
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
  v <- eigen(sym(gsir))$vectors[,1:r]
  
  if(is.null(Kx_new)){
    if(xtype=="scalar") Kx_new <- gram_gauss_eucl(x,x_new,complex_x)
    else if(xtype=='distribution') Kx_new <- gram_gauss_l2(x,x_new,complex_x)
  }
  
  pred.est <- t(t(v)%*%Gxinv%*%Q%*%Kx_new)
  
  return(list(pred.est=pred.est,mat=gsir))
}


###############
## generalized CV to tune epsilon_x and epsilon_y
gcv_l2 <- function(x, y, eps, which, ytype, xtype, 
                   complex_x, complex_y, Kx=NULL, Ky=NULL){
  n <- dim(x)[1]
  if(is.null(Kx)){
    if(xtype=="scalar") Kx <- gram_gauss_eucl(x, x, complex_x)
    else if(xtype=='distribution') Kx <- gram_l2(x, complex_x)
  }
  if(is.null(Ky)){
    if(ytype=="scalar") Ky <- gram_gauss_eucl(y, y, complex_y)
    else if(ytype=="distribution") Ky <- gram_l2(y, complex_y)
  }
  if(which=="ey") {G1 <- Kx;G2 <- Ky}
  else if(which=="ex") {G1 <- Ky;G2 <- Kx}
  
  G2inv <- solve(G2+eps*onorm(G2)*diag(n))
  nu <- sum((G1-G2%*%G2inv%*%G1)^2)
  de <- (1-sum(diag(G2inv%*%G2))/n)^2
  
  return(nu/de)
}

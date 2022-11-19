#################################
## BIC for order determination ##
#################################

bic <- function(x, y, ytype, xtype, atype, ex, ey, complex_x, 
                complex_y, criterion, c0=2){
  ## Param
  ##  x: the n*m matrix of the predictors; 
  ##  y: the response object;
  ##  ytype: string{"distribution", "scalar"}, indicating the type of response;
  ##  atype: string "identity" or "Gyinv", indicating the type of matrix A;
  ##  ex,ey: control the regularization magnitude;
  ##  complex_x and complex_y: measure the complexity of HX and HY;
  ##  criterion: string {'lal' or 'zmp'}
  ## Return:
  ##  list of estimated sufficient predictors and GSIR candidate matrix
  n <- dim(x)[1]
  candmat <- wgsir(x=x,x_new=x,y=y,xtype=xtype,ytype=ytype,atype=atype,
                   ex=ex,ey=ey,complex_x=complex_x,complex_y=complex_y,r=n)$mat
  
  candmat <- (candmat + t(candmat))/2
  out <- eigen(candmat)
  lam <- out$values
  if(criterion=="lal"){
    gn <- c(0)
    for(k in 1:n){
        gn <- c(gn, sum(lam[1:k]) - (c0*lam[1])*n^(-1/2)*(log(n))^(1/2)*k)}
  }
  if(criterion=="zmp"){
    gn <- numeric()
    for(k in 0:(n-1)){
      c1 = (lam[1]/3)*(0.5* log(n)+0.1* n^(1/3))/(2*n)
      c2 = k*(2*n-k+1)
      gn = c(gn, sum(log(lam[(k+1):n]+1)-lam[(k+1):n])-c1*c2)
    }
    gn = c(gn,-c1*p*(2*p-p+1))
  }
  
  return(list(rhat = which.max(gn)-1, rcurve = gn))
}




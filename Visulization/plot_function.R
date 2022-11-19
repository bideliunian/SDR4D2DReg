############################################################# 
##  plotting function for SDR in s2d or d2d regression ######
###############################################################

plot_wgsir <- function(n, m, p=NULL, model, complex.x, complex.y){
  
  result <- list()## result container
  
  ### generate data
  data <- gendata(n=n, m =m, p=p, model = model)
  
  n <- data$n; m <- data$m; d0 <- data$d0
  ytype <- data$ytype #'scalar', 'distribution'
  xtype <- data$xtype
  pred <- as.matrix(data$pred)[1:floor(n/2),]%>%as.matrix()
  pred.test <- as.matrix(data$pred)[(floor(n/2)+1):n,]%>%as.matrix()
  
  if(ytype=='scalar') {
    ytrain <- data$y[1:floor(n/2)]
    ytest <- data$y[(floor(n/2)+1):n]
    Ky <- gram_gauss_eucl(ytrain, ytrain, complex.y)
  }
  
  else if(ytype=='distribution') {
    ytrain <- data$y[1:floor(n/2),];ytest <- data$y[(floor(n/2)+1):n,]
    Ky <- gram_w(ytrain, complex.y)
  } 
  
  else if(ytype=='mult_distribution') {
    ytrain <- data$y[1:floor(n/2),,];ytest <- data$y[(floor(n/2)+1):n,,]
    Ky <- gram_sw(ytrain, complex.y)
  }
  
  if(xtype=="scalar") {
    xtrain <- data$x[1:floor(n/2)];xtest <- data$x[(floor(n/2)+1):n]
    Kx <- gram_gauss_eucl(xtrain,xtrain,complex.x)
    Kx.new <- gram_gauss_eucl(xtrain,xtest,complex.x)
  }
  
  else if(xtype=='distribution') {
    xtrain <- data$x[1:floor(n/2),];xtest <- data$x[(floor(n/2)+1):n,]
    Kx <- gram_w(xtrain,complex.x)
    Kx.new <- gram_gauss_w(xtrain,xtest,complex.x)
  }
  
  else if(xtype=='mult_distribution') {
    xtrain <- data$x[1:floor(n/2),,];xtest <- data$x[(floor(n/2)+1):n,,]
    Kx <- gram_sw(xtrain,complex.x)
    Kx.new <- gram_gauss_sw(xtrain,xtest,complex.x)
  }
  
  eps.grid <- 10^(-6:0)
  epsx <- lapply(eps.grid, FUN = gcv, x=Kx, y=Ky, which='ex', ytype='gram', xtype='gram',
                 complex_x=complex.x,complex_y=complex.y)
  epsx.min <- eps.grid[which.min(epsx)]
  
  epsy <- lapply(eps.grid, FUN = gcv, x=Kx, y=Ky, which='ey', ytype='gram', xtype='gram',
                 complex_x=complex.x,complex_y=complex.y)
  epsy.min <- eps.grid[which.min(epsy)]
  
  # using true dimension
  r0_iden <- d0
  r0_yinv <- d0
  
  result$estpred.sir1 <- wgsir(x=Kx, x_new=Kx.new, y=Ky, xtype='gram', ytype='gram', atype='identity',
                               ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
                               r=r0_iden)$pred.est
  result$estpred.sir2 <- wgsir(x=Kx, x_new=Kx.new, y=Ky, xtype='gram', ytype='gram', atype='Gyinv',
                               ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y,
                               r=r0_yinv)$pred.est
  
  if(ytype=='distribution') {
    result$density <- apply(X=ytest, MARGIN = 1, FUN = density)
  }
  else if(xtype=='mult_distribution') {
    ytest.list <- list()
    for (i in 1:floor(n/2)) {
      ytest.list[[i]] = ytest[i,,]
    }
    result$ytest = ytest.list
    result$density <- lapply(X=ytest.list, FUN = function(x) return(kde2d(x[,1], x[,2], h = rep(0.1, 2), n = 25))) 
  }
  
  result$truepred <- pred.test
  
  return(result)
}

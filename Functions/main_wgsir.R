#################################################
# main function to implement wgsir ###############
###################################################
main_wgsir <- function(times, n, m, p=NULL, model, complex.x, complex.y, order.method, seed){
  # Param
  #   times: times of repeatment of experiment
  #   n: sample size
  #   m: the number of iid observations from a distribution
  #   p: dimension of predictors if X is vector
  #   model: string in 's2d_ex1-6', 'd2d_ex1-6', 'd2s_ex1-3', 's2s_ex1', 'md2md_ex1-5'
  #   order_method: string in 'true', 'bic', 'ladle'
  # Return
  #   list of estimation error and estimated dimension for wgsir1 and wgsir2
  print(paste("-----------------------n=", n, "; m=", m , "; model=", model,"----------------------"))
  result <- list(wgsir1.rvmr=vector(),wgsir2.rvmr=vector(),
                 wgsir1.dcor=vector(),wgsir2.dcor=vector(),
                 wgsir1.dim=vector(),wgsir2.dim=vector(),
                 wgsir1.correct=vector(),wgsir2.correct=vector()) # result container
  set.seed(seed)
  for (i in 1:times) {
    cat(sprintf("-----------------%i th experiment------------", i))
    tic('Time of Experiment')
    
    ## generate data
    data <- gendata(n=n, m =m, p=p, model = model)
    n <- data$n 
    m <- data$m 
    d0 <- data$d0
    ytype <- data$ytype #'scalar', 'distribution', 'mult_distribution'
    xtype <- data$xtype # 'scalar', distribution', 'mult_distribution'
    pred <- as.matrix(as.matrix(data$pred)[1:floor(n/2),])
    pred.test <- as.matrix(as.matrix(data$pred)[(floor(n/2)+1):n,])
    
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
    
    if (order.method == 'bic'){
      r0_iden <- bic(x=Kx, y=Ky, xtype='gram', ytype='gram', atype='identity',
                     ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
                     criterion='lal', c0=2)$rhat
      r0_yinv <- bic(x=Kx, y=Ky, xtype='gram', ytype='gram', atype='Gyinv',
                     ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
                     criterion='lal', c0=4)$rhat  
    }
    else if (order.method == 'ladle'){
      r0_iden <- ladle(Kx=Kx, Ky=Ky, nboot = 100, atype='identity', ex=epsx.min, ey=epsy.min)$rhat
      r0_yinv <- ladle(Kx=Kx, Ky=Ky, nboot = 100, atype='Gyinv', ex=epsx.min, ey=epsy.min)$rhat  
    }
    else if (order.method == 'true'){
      r0_iden <- d0
      r0_yinv <- d0
    }
    
    estpred.sir1 <- wgsir(x=Kx, x_new=Kx.new, y=Ky, xtype='gram', ytype='gram', atype='identity',
                          ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
                          r=r0_iden)$pred.est
    estpred.sir2 <- wgsir(x=Kx, x_new=Kx.new, y=Ky, xtype='gram', ytype='gram', atype='Gyinv',
                          ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y,
                          r=r0_yinv)$pred.est
    
    result$wgsir1.rvmr[i] <- rvmr(estpred.sir1, pred.test)
    result$wgsir2.rvmr[i] <- rvmr(estpred.sir2, pred.test)
    
    estpred.sir1.dist <- dist(estpred.sir1)
    estpred.sir2.dist <- dist(estpred.sir2)
    pred.test.dist <- dist(pred.test)
    result$wgsir1.dcor[i] <- dcor(estpred.sir1.dist, pred.test.dist)
    result$wgsir2.dcor[i] <- dcor(estpred.sir2.dist, pred.test.dist)
    
    result$wgsir1.dim[i] <- r0_iden
    result$wgsir2.dim[i] <- r0_yinv
    
    result$wgsir1.correct[i] <- as.numeric(d0 == r0_iden)
    result$wgsir2.correct[i] <- as.numeric(d0 == r0_yinv)
    toc()
  }
  return(result)
}

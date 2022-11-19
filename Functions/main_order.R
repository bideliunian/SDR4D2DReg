#################################################################
########## main function to implement order determination #######
#################################################################
main_order <- function(times, n, m, p=NULL, model, complex.x, complex.y, order.method){
  # Param
  #   times: times of repeatment of experiment
  #   n: sample size
  #   m: the number of iid observations from a distribution
  #   p: dimension of predictors if X is vector
  #   model: string in 's2d_ex1-6', 'd2d_ex1-6', 'd2s_ex1-3', 's2s_ex1', 'md2md_ex1-5'
  #   order_method: string in 'bic', 'ladle'
  # Return
  #   list of estimation error and estimated dimension for wgsir1 and wgsir2
  print(paste("-----------------------n=", n, "; m=", m , "; model=", model,"----------------------"))
  result <- list(wgsir1.dim=vector(),wgsir2.dim=vector(),
                 wgsir1.correct=vector(),wgsir2.correct=vector()) # result container
  for (i in 1:times) {
    cat(sprintf("-----------------%i th experiment------------", i))
    tic('Time of Experiment')
    
    ## generate data
    data <- gendata(n=n, m =m, p=p, model = model)
    y <- data$y
    x <- data$x
    n <- data$n 
    m <- data$m 
    d0 <- data$d0
    ytype <- data$ytype #'scalar', 'distribution', 'mult_distribution'
    xtype <- data$xtype # 'scalar', distribution', 'mult_distribution'
    pred <- data$pred
    
    if(ytype=='scalar') {
      Ky <- gram_gauss_eucl(y, y, complex.y)
    }
    else if(ytype=='distribution') {
      Ky <- gram_w(y, complex.y)
    } 
    else if(ytype=='mult_distribution') {
      Ky <- gram_sw(y, complex.y)
    }
    
    if(xtype=="scalar") {
      Kx <- gram_gauss_eucl(x,x,complex.x)
    }
    else if(xtype=='distribution') {
      Kx <- gram_w(x,complex.x)
    }
    else if(xtype=='mult_distribution') {
      Kx <- gram_sw(x, complex.x)
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
    
    result$wgsir1.dim[i] <- r0_iden
    result$wgsir2.dim[i] <- r0_yinv
    
    result$wgsir1.correct[i] <- as.numeric(d0 == r0_iden)
    result$wgsir2.correct[i] <- as.numeric(d0 == r0_yinv)
    toc()
  }
  return(result)
}

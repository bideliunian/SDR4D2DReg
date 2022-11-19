## generate sample from truncated gamma
rgammat <- function(n, range=c(0.2,2), shape, scale = 1) {
  F.a <- pgamma(min(range), shape = shape, scale = scale)
  F.b <- pgamma(max(range), shape = shape, scale = scale)
  
  u <- runif(n, min = F.a, max = F.b)
  return(qgamma(u, shape = shape, scale = scale))
}

rnormalt <- function(n, range=c(1,3), mu, sigma){
  F.a <- pnorm(min(range), mean = mu, sd = sigma)
  F.b <- pnorm(max(range), mean = mu, sd = sigma)
  
  u <- runif(n, min = F.a, max = F.b)
  return(qnorm(u, mean = mu, sd = sigma))
}

## s2d---scalar to distribution model
## d2d---distribution to distribution model
## md2md--- multidim distribution to multidim_distribution model

gendata <- function(n, m, p=NULL, model){
  ## Input: n is sample_size; m is # of discrete observations from a x/y;
  ##        p is the dimension of predictors in d2s model
  ##        model is 's2d_ex1-6', 'd2d_ex1-6', 'd2s_ex1-3', 's2s_ex1', 'md2md_ex1-5'
  ##        
  ## Output: d2d models:data$x- n*m matrix; data$y- n*m matrix for distributional data;
  ##         s2d models:data$x- n*m matrix; data$y- n*1 matrix.
  ##         md2md models: data$x- n*m*d array; data$y- n*m*dim array
  data <- list()
  #################
  ## scalar on distribution regression
  ##################################
  if(model=='s2d_ex1'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    dw1 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=1)
    dw2 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=3)
    data$pred <- exp(dw1^2)+exp(dw2^2)
    data$y <- rnorm(n, mean = data$pred, sd=0.1)
    data$d0 <- 1
    data$ytype <- 'scalar'
    data$xtype <- 'distribution'
  }
  else if(model=='s2d_ex2'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    dw1 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=1)
    dw2 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=3)
    data$pred <- exp(dw1)/(1+exp(dw2))
    data$y <- rnorm(n, mean = data$pred, sd=0.1)
    data$d0 <- 1
    data$ytype <- 'scalar'
    data$xtype <- 'distribution'
  }
  else if(model=='s2d_ex3'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    dw1 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=1)
    dw2 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=3)
    data$pred <- sin(pi*(exp(dw1+dw2))/2)
    data$y <- rnorm(n, mean = data$pred, sd=0.1)
    data$d0 <- 1
    data$ytype <- 'scalar'
    data$xtype <- 'distribution'
  }
  else if(model=='s2d_ex4'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    d1 <- mapply(hellinger_dist4beta, a1=a, b1=b, a2=2, b2=1)
    d2 <- mapply(hellinger_dist4beta, a1=a, b1=b, a2=2, b2=3)
    data$pred <- exp(d1)+exp(d2)
    data$y <- rnorm(n, mean = data$pred, sd=0.1)
    data$d0 <- 1
    data$ytype <- 'scalar'
    data$xtype <- 'distribution'
  }
  else if(model=='s2d_ex5'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    d1 <- mapply(kl_dist4beta, a1=a, b1=b, a2=2, b2=1)
    d2 <- mapply(kl_dist4beta, a1=2, b1=3, a2=a, b2=b)
    data$pred <- exp(d2)
    data$y <- rnorm(n, mean = data$pred, sd=0.1)
    data$d0 <- 1
    data$ytype <- 'scalar'
    data$xtype <- 'distribution'
  }
  else if(model=='s2d_ex6'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    d1 <- a/(a+b)
    d2 <- (a+1)*a/(a+b+1)/(a+b)
    data$pred <- (d1+d2)
    data$y <- rnorm(n, mean = data$pred, sd=0.1)
    data$d0 <- 1
    data$ytype <- 'scalar'
    data$xtype <- 'distribution'
  }
  
  ###########################
  ## univariate distribution on uni-distribution regression
  ###############################################
  else if(model=='d2d_ex1'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    # a <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    # b <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    dw1 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=1)
    dw2 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=3)
    data$pred <- exp(dw1^2)+exp(dw2^2)
    mu <- rnorm(n, mean = data$pred, sd = 0.2)
    sigma <- 1
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$d0 <- 1
    data$ytype <- 'distribution'
    data$xtype <- 'distribution'
  }
  else if(model=='d2d_ex2'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    # a <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    # b <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    dw1 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=1)
    dw2 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=3)
    data$pred <- cbind(dw1,dw2)
    mu <- rnorm(n, mean = exp(dw1^2), sd = 0.2)
    sigma <- exp(dw2^2)
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$d0 <- 2
    data$ytype <- 'distribution'
    data$xtype <- 'distribution'
  }
  else if(model=='d2d_ex3'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    # a <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    # b <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    dw1 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=1)
    dw2 <- mapply(dist4beta_w, a1=a, b1=b, a2=2, b2=3)
    data$pred <- cbind(dw1,dw2)
    mu <- rnorm(n, mean = exp(dw1^2), sd = 0.2)
    sigma <- rgamma(n,shape=dw2^2,rate=dw2)
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$d0 <- 2
    data$ytype <- 'distribution'
    data$xtype <- 'distribution'
  }
  else if(model=='d2d_ex4'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    # a <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    # b <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    dw1 <- mapply(hellinger_dist4beta, a1=a, b1=b, a2=2, b2=1)
    dw2 <- mapply(hellinger_dist4beta, a1=a, b1=b, a2=2, b2=3)
    data$pred <- cbind(dw1,dw2)
    mu <- rnorm(n, mean = exp(dw1), sd = 0.2)
    sigma <- exp(dw2)
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$d0 <- 2
    data$ytype <- 'distribution'
    data$xtype <- 'distribution'
  }
  else if(model=='d2d_ex5'){
    a <- rgamma(n, shape = 2, rate = 1)
    b <- rgamma(n, shape = 2, rate = 3)
    # a <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    # b <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    d1 <- mapply(kl_dist4beta, a1=a, b1=b, a2=2, b2=1)
    d2 <- mapply(kl_dist4beta, a1=2, b1=3, a2=a, b2=b)
    data$pred <- d2
    mu <- rnorm(n, mean = sqrt(exp(d2)), sd = 0.2)
    sigma <- 1
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$d0 <- 1
    data$ytype <- 'distribution'
    data$xtype <- 'distribution'
  }
  else if(model=='d2d_ex6'){
    # a <- rgamma(n, shape = 2, rate = 1)
    # b <- rgamma(n, shape = 2, rate = 3)
    a <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    b <- rnormalt(n, mu = 2, sigma = 1, range = c(1, 3))
    data$x <- t(mapply(rbeta, n=m, shape1=a, shape2=b))
    d1 <- a/(a+b)
    d2 <- (a+1)*a/(a+b+1)/(a+b)
    data$pred <- cbind(d1,d2)
    mu <- rnorm(n, mean = d1, sd = 0.2)
    sigma <- rgamma(n,shape=d2,rate=sqrt(d2))
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$d0 <- 2
    data$ytype <- 'distribution'
    data$xtype <- 'distribution'
  }
  
  ###############################
  ## uni-distribution on scalar regression/ frechet regression
  ######################################
  else if(model=='d2v_ex1'){
    data$x <- matrix(rnorm(n*p, 0, 1), nrow = n)
    data$pred <- data$x[,1]^2+data$x[,2]^2#*log(data$x[,1]^2+data$x[,2]^2)
    mu <- rnorm(n, mean = data$pred, sd = 0.2)
    sigma <- 1
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$d0 <- 1
    data$ytype <- 'distribution'
    data$xtype <- 'scalar'
  }
  else if(model=='d2v_ex2'){
    data$x <- matrix(rnorm(n*p, 0, 1), nrow = n)
    d1 <- data$x[,1]/(1+exp(data$x[,2]))
    mu <- rnorm(n, mean = d1, sd = 0.2)
    d2<- sqrt(data$x[,3]^2+data$x[,4]^2)
    sigma <- d2
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$pred <- cbind(d1,d2)
    data$d0 <- 2
    data$ytype <- 'distribution'
    data$xtype <- 'scalar'
  }
  else if(model=='d2v_ex3'){
    data$x <- matrix(rnorm(n*p, 0, 1), nrow = n)
    d1 <- data$x[,1]^2+data$x[,2]^2
    mu <- rnorm(n, mean = d1, sd = 0.2)
    d2 <- sqrt(data$x[,3]^2+data$x[,4]^2)
    sigma <- rgamma(n,shape=d2^2, rate=d2)
    data$y <- t(mapply(rnorm, n=m, mean=mu, sd=sigma))
    data$pred <- cbind(d1,d2)
    data$d0 <- 2
    data$ytype <- 'distribution'
    data$xtype <- 'scalar'
  }
  
  ####################
  ## scalar on scalar regression
  ##############################
  else if(model=='s2s_ex1'){
    data$x <- matrix(rnorm(n*p, 0, 1), nrow = n)
    d1 <- sqrt(data$x[,1]^2+data$x[,2]^2)*log(sqrt(data$x[,1]^2+data$x[,2]^2))
    data$y <- rnorm(n, mean = d1, sd = 0.5)
    data$pred <- d1
    data$d0 <- 1
    data$ytype <- 'scalar'
    data$xtype <- 'scalar'
  }
  
  #########################
  ## multivariate distribution data
  ######################################
  
  else if(model=='md2md_ex1'){
    a <- rnorm(n, 0.5, 0.5)
    b <- rbeta(n, 2, 3)
    data$x <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    data$x[,,1] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    data$x[,,2] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    d1 <- rep(0,n); d2 <- rep(0,n)
    for(i in 1:n){
      d1[i] <- dist4normal_md(m1=a[i]*c(1,1), sigma1=b[i]*diag(2), m2=c(-1,0), sigma2=diag(c(1,1/2)))
      d2[i] <- dist4normal_md(m1=a[i]*c(1,1), sigma1=b[i]*diag(2), m2=c(0,1), sigma2=diag(c(1/2,1)))
    } 
    data$pred <- d1
    mu1 <- rnorm(n, (d1), 1)
    mu2 <- rnorm(n, (d1), 1)
    data$y <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    for(i in 1:n){
      data$y[i,,] <- mvrnorm(m, mu = c(mu1[i],mu2[i]), Sigma = diag(2))
    }
    data$d0 <- 1
    data$ytype <- 'mult_distribution'
    data$xtype <- 'mult_distribution'
  }
  
  else if(model=='md2md_ex2'){
    a <- rnorm(n, 0.5, 0.5)
    b <- rbeta(n, 2, 3)
    data$x <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    data$x[,,1] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    data$x[,,2] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    d1 <- rep(0,n); d2 <- rep(0,n)
    for(i in 1:n){
      d1[i] <- dist4normal_md(m1=a[i]*c(1,1), sigma1=b[i]*diag(2), m2=c(-1,0), sigma2=diag(c(1,1/2)))
      d2[i] <- dist4normal_md(m1=a[i]*c(1,1), sigma1=b[i]*diag(2), m2=c(0,1), sigma2=diag(c(1/2,1)))
    } 
    data$pred <- cbind(d1,d2)
    Gam <- matrix(c(1,1,-1,1), nrow = 2, byrow = TRUE)/sqrt(2)
    # lambda1 <- rgamma(n, shape=d2^2, rate=d2)
    # lambda2 <- rgamma(n, shape=d2^2, rate=d2)
    lambda1 <- abs(rnorm(n, mean = d2, sd = 0.5))
    lambda2 <- abs(rnorm(n, mean = d2, sd = 0.5))
    mu <- c(d1, d1)
    data$y <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    for(i in 1:n){
      data$y[i,,] <- mvrnorm(m, mu = c(sqrt(d1[i]),sqrt(d1[i])), Sigma = Gam%*%diag(c(lambda1[i], lambda2[i]))%*%t(Gam))
    }
    data$d0 <- 2
    data$ytype <- 'mult_distribution'
    data$xtype <- 'mult_distribution'
  }
  
  else if(model=='md2md_ex3'){
    a <- rnorm(n, 0.5, 0.5)
    b <- rbeta(n, 2, 3)
    data$x <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    data$x[,,1] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    data$x[,,2] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    d1 <- rep(0,n); d2 <- rep(0,n)
    for(i in 1:n){
      d1[i] <- dist4normal_md(m1=a[i]*c(1,1), sigma1=b[i]*diag(2), m2=c(-1,0), sigma2=diag(c(1,1/2)))
      d2[i] <- dist4normal_md(m1=a[i]*c(1,1), sigma1=b[i]*diag(2), m2=c(0,1), sigma2=diag(c(1/2,1)))
    } 
    data$pred <- cbind(d1,d2)
    Gam <- matrix(c(1,1,-1,1), nrow = 2, byrow = TRUE)/sqrt(2)
    lambda1 <- rgammat(n, shape=d2^2, scale=1/d2)
    lambda2 <- rgammat(n, shape=d2^2, scale=1/d2)
    mu1 <- rnorm(n, d1, 1)
    mu2 <- rnorm(n, d1, 1)
    data$y <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    for(i in 1:n){
      data$y[i,,] <- mvrnorm(m, mu = c(mu1[i],mu2[i]), Sigma = Gam%*%diag(c(lambda1[i], lambda2[i]))%*%t(Gam))
    }
    data$d0 <- 2
    data$ytype <- 'mult_distribution'
    data$xtype <- 'mult_distribution'
  }
  else if(model=='md2md_ex4'){
    a <- rnorm(n, 0, 1)
    b <- rbeta(n, 2, 3)
    data$x <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    data$x[,,1] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    data$x[,,2] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    d1 <- rep(0,n); d2 <- rep(0,n)
    for(i in 1:n){
      d1[i] <- hellinger_dist4normal(mu1=a[i]*c(1,1), sigma1=b[i]*diag(2), mu2=c(-1,0), sigma2=diag(c(1,1/2)))
      d2[i] <- hellinger_dist4normal(mu1=a[i]*c(1,1), sigma1=b[i]*diag(2), mu2=c(0,1), sigma2=diag(c(1/2,1)))
    }
    data$pred <- cbind(d1,d2)
    Gam <- matrix(c(1,1,1,-1), nrow = 2, byrow = TRUE)/sqrt(2)
    lambda1 <- rgammat(n, shape=d2^2, scale=1/d2)
    lambda2 <- rgammat(n, shape=d2^2, scale=1/d2)
    mu1 <- rnorm(n, d1, 1)
    mu2 <- rnorm(n, d1, 1)
    data$y <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    for(i in 1:n){
      data$y[i,,] <- mvrnorm(m, mu = c(mu1[i],mu2[i]), Sigma = Gam%*%diag(c(lambda1[i], lambda2[i]))%*%t(Gam))
    }
    data$d0 <- 2
    data$ytype <- 'mult_distribution'
    data$xtype <- 'mult_distribution'
  }
  else if(model=='md2md_ex5'){
    a <- rnorm(n, 0, 1)
    b <- rbeta(n, 2, 3)
    data$x <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    data$x[,,1] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    data$x[,,2] <- t(mapply(rnorm, n=m, mean=a, sd=b))
    d1 <- rep(0,n); d2 <- rep(0,n)
    for(i in 1:n){
      d1[i] <- dist4normal_md(m1=a[i]*c(1,1), sigma1=b[i]*diag(2), m2=c(-1,0), sigma2=diag(c(1,1/2)))
      d2[i] <- dist4normal_md(m1=a[i]*c(1,1), sigma1=b[i]*diag(2), m2=c(0,1), sigma2=diag(c(1/2,1)))
    } 
    data$pred <- cbind(d1,d2)
    Gam <- matrix(c(1,1,-1,1), nrow = 2, byrow = TRUE)/sqrt(2)
    lambda1 <- abs(rnorm(n, mean = d2, sd = 0.1))
    lambda2 <- abs(rnorm(n, mean = d2, sd = 0.1))
    mu1 <- rnorm(n, d1, 1)
    mu2 <- rnorm(n, d1, 1)
    data$y <- array(rep(NA, n*m*2), dim = c(n, m ,2))
    for(i in 1:n){
      data$y[i,,] <- mvrnorm(m, mu = c(mu1[i],mu2[i]), Sigma = Gam%*%diag(c(lambda1[i], lambda2[i]))%*%t(Gam))
    }
    data$d0 <- 2
    data$ytype <- 'mult_distribution'
    data$xtype <- 'mult_distribution'
  }
  data$n <- n
  data$m <- m
  return(data)
}

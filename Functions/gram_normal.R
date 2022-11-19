
## l2 distance between two normal distribution
dist4normal_l2 <- function(mu1, sigma1, mu2, sigma2){
  # Param
  #   dist1: N(mu1, sigma1)
  #   dist2: N(mu2, sigma2)
  # Return 
  #   l2 distance between two normal distribution
  d <- sqrt(abs((1/(2*sqrt(pi)))*(1/sigma1+1/sigma2)
                -2*dnorm(x = mu1,mean = mu2,sd = sqrt(sigma1^2+sigma2^2))))
  return(d)
}

## gram matrix with l2 distance gaussian kernel between two group of normal
gram_normal_l2 <- function(mu1, sigma1, mu2, sigma2, complexity){
  n <- length(mu1)
  m <- length(mu2)
  k <- outer(
    1:n, 1:m,
    Vectorize(function(i,k) dist4normal_l2(mu1[i], sigma1[i], mu2[k], sigma2[k]))
  )
  sigma <- sum(k^2)/choose(n,2)
  gamma <- complexity/(sigma)
  return(exp(-gamma*(k^2)))
}

## wasserstein distance between two gaussian measures
dist4normal_w <- function(a1, b1, a2, b2){
  return((a1 - a2)^2+(b1 - b2)^2)
}

## gram matrix with wasserstein distance gaussian kernel between two group of normal
gram_normal_w <- function(mu1, sigma1, mu2, sigma2, complexity){
  n <- length(mu1)
  m <- length(mu2)
  k <- outer(
    1:n, 1:m,
    Vectorize(function(i,k) dist4normal_w(mu1[i], sigma1[i], mu2[k], sigma2[k]))
  )
  sigma <- sum(k^2)/choose(n,2)
  gamma <- complexity/(sigma)
  return(exp(-gamma*(k^2)))
}


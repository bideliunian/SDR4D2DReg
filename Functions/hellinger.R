## calculate hellinger distance between two beta distributions
hellinger_dist4beta <- function(a1,b1,a2,b2){
  return(1-beta((a1+a2)/2,(b1+b2)/2)/sqrt(beta(a1,b1)*beta(a2,b2)))
}

## calculate hellinger distance between two normal distributions
hellinger_dist4normal <- function(mu1,sigma1,mu2,sigma2){
  lambda1 = eigen(sigma1)$values
  lambda2 = eigen(sigma2)$values
  d1 = det(sigma1); d2 = det(sigma2); d3 = det(sigma1/2+sigma2/2)
  return(as.numeric(1-(d1*d2)^{1/4}/sqrt(d3)*exp(-t(mu1-mu2)%*%solve(sigma1/2+sigma2/2)%*%(mu1-mu2))))
}

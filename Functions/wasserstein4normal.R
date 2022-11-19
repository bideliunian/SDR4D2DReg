## calculate wass distance between two multi-dim gaussian measures
dist4normal_md <- function(m1,sigma1,m2,sigma2){
  # Param
  #   wasserstein distance between N(m1, sigma1) and N(m2, sigma2)
  sigma1_sqrt <- matpower(sigma1, 1/2)
  sigma2_sqrt <- matpower(sigma2, 1/2)
  return(sum(m1-m2)^2+trace(sigma1+sigma2-2*matpower(sigma1_sqrt%*%sigma2%*%sigma1_sqrt, 1/2)))
}
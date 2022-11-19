
## calculate kl divergence between two beta distributions
kl_dist4beta <- function(a1,b1,a2,b2){
  dgamma1 <- digamma(a1+b1)
  return(log(beta(a2,b2)/beta(a1,b1))+(a1-a2)*(digamma(a1)-dgamma1)+(b1-b2)*(digamma(b1)-dgamma1))
}

## calculate wass distance between two gaussian measures
dist4normal <- function(a1,b1,a2,b2){
  return((a1-a2)^2+(b1-b2)^2)
}
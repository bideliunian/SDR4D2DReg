###################
## This file define functions computing sliced wasserstein distance
#######################

############## 
## get random projections
#############
get_random_projections <- function(n_projections, d, n_seed){
  ## Generates n_projections samples from the uniform on the unit sphere of dimension d-1
  ## Parameter:
  ##  n_projections: number of projections
  ##  d: dimension
  ##  seed: seed_number
  ## Return: 
  ##  n_projections * d matrix, each row is a random projection direction
  set.seed(seed = n_seed)
  projections <- matrix(rnorm(n = n_projections*d, mean = 0, sd = 1), nrow = n_projections, ncol = d)
  projections <- apply(projections, MARGIN = 2, FUN = function(x) x/sqrt(sum(x^2)))
  return(projections)
}


sliced_wass_distance <- function(X_s, X_t, n_projections=50, n_seed = 1){
  # Computes a Monte-Carlo approximation of the 2-Sliced Wasserstein distance
  ###############
  ## Parameters
  ##  X_s : n_samples_a * dim
  ##  X_t : n_samples_b * dim
  ##  n_projections : Number of projections used for the Monte-Carlo approximation
  ## Returns
  ##  res: Sliced Wasserstein distance
  X_s <- as.matrix(X_s)
  X_t <- as.matrix(X_t)
  X_s <- na.omit(X_s)
  X_t <- na.omit(X_t)
  
  n <- dim(X_s)[1]
  m <- dim(X_t)[1]
  
  if (dim(X_s)[2] != dim(X_t)[2]){
    stop("X_s and X_t must have the same number of dimensions")
  }
  d <- dim(X_s)[2]
  
  projections <- get_random_projections(n_projections, d, n_seed = n_seed)
  
  X_s_projections <- X_s%*%t(projections)
  X_t_projections <- X_t%*%t(projections)
  
  projected_emd <- NULL
  
  res = 0
  
  for (i in 1:n_projections){
    emd <- wasserstein1d(X_s_projections[,i], X_t_projections[,i])
    if (is.null(projected_emd)){
      projected_emd[i] <- emd^2 
    }
    res <- emd^2 + res
  }
  
  res <- sqrt(res/n_projections)
  return(res)
}


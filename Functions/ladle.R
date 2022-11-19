#############################################################################
#################  Ladle estimator for order determination #############
###########################################################################

ladle <- function(Kx, Ky, nboot, p=NULL, atype, ex, ey){
  # Para:
  #   Kx: n by n gram matrix for predictor
  #   Ky: n by n gram matrix for response
  #   nboot: number of times for bootstrap
  #   p: dimension of predictor for distribution on vector regression
  n <- dim(Kx)[1]
  if (is.null(p)) {p <- n}
  # kmax is the upper bound of candidate orders
  if(p <= 10) {
    kmax <- p-2
  }
  else {
    kmax <- round(p/log(p)) 
  }

  candidate_matrix <- wgsir(x=Kx, x_new=Kx, y=Ky, xtype='gram', ytype='gram', 
               atype=atype, ex=ex, ey=ey, r=1)$mat
  eigen_val_old <- eigen(candidate_matrix)$values  # eigenvalues
  eigen_vec_old <- eigen(candidate_matrix)$vectors  # eigenvectors
  
  # compute phi_n
  eps = 1 # adding 1 for robustness;
  phi <- function(kmax, eigen_val){
    eigen_all <- eps + sum(eigen_val[1:(kmax+1)]) 
    # warning: adding eps may affect the relative magnitude between phi and f parts, 
    # need more careful design
    return(eigen_val[1:(kmax+1)]/eigen_all)
  }
  phi_n <- phi(kmax, eigen_val_old)
  
  f <- function(kmax, evec1, evec2){
    out <- c(0)
    for(k in 1:kmax){
      if(k == 1) out <- c(out, 1 - abs(t(evec1[,1])%*%evec2[,1]))
      if(k!=0 & k!=1) out <- c(out, 1 - abs(det(t(evec1[,1:k])%*%evec2[,1:k])))
    }
    return(out)
  }
  fn0 <- 0
  for(iboot in 1:nboot){
    bootindex <- sample(1:n, n, replace=TRUE)
    Kxs <- Kx[bootindex, bootindex]
    Kys <- Ky[bootindex, bootindex]
    mat <- wgsir(x=Kxs, x_new=Kxs, y=Kys, xtype='gram', ytype='gram', 
                 atype=atype, ex=ex, ey=ey, r=1)$mat
    eigen_val <- eigen(mat)$values
    eigen_vec <- eigen(mat)$vectors
    fn0 <- fn0 + f(kmax, eigen_vec_old, eigen_vec) / nboot
  }
  f_n <- fn0 / (eps + sum(fn0)) # again, here adding 1 for robustness
  
  g_n <- phi_n + f_n
  rhat <- which.min(g_n) - 1
  
  return(list(rset=(0:kmax), gn=g_n, rhat=rhat))
}

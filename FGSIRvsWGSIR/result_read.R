#######################################################
###########   Simulations for Scenario III: d2d  ########
########    including Model II-1-6##########   

#################### Preparation #################################

function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/FGSIRvsWGSIR"
save_path <- "~/work/DR4D2D/FGSIRvsWGSIR/Results"

require(dplyr)

# set global parameters

#####  global parameters
n <- 200
m <- 100
Models <- c('s2d_ex1','s2d_ex2','s2d_ex3','d2d_ex1','d2d_ex2', 'd2d_ex3', 'd2d_ex4')
complex.x <- 1
complex.y <- 1
order.method <- 'bic'

gsir.l1.rvmr <- numeric();gsir.l2.rvmr <- numeric()
gsir.l1.dim <- numeric();gsir.l2.dim <- numeric()
gsir.l1.dcor <- numeric();gsir.l2.dcor <- numeric()
gsir.l1.correct <- numeric();gsir.l2.correct <- numeric()

############################### read data ##############################
for (i in 1:100) {
  load(file = paste(save_path,"/l1gsir.",i,".RData",sep=""))
  load(file = paste(save_path,"/l2gsir.",i,".RData",sep=""))
  
  result_gsir.l1.rvmr <- result_gsir_l1[1, ]%>%unlist()
  result_gsir.l2.rvmr <- result_gsir_l2[1, ]%>%unlist()
  result_gsir.l1.dcor <- result_gsir_l1[2, ]%>%unlist()
  result_gsir.l2.dcor <- result_gsir_l2[2, ]%>%unlist()
  result_gsir.l1.dim <- result_gsir_l1[3, ]%>%unlist()
  result_gsir.l2.dim <- result_gsir_l2[3, ]%>%unlist()
  result_gsir.l1.correct <- result_gsir_l1[4, ]%>%unlist()
  result_gsir.l2.correct <- result_gsir_l2[4, ]%>%unlist()
  
  gsir.l1.rvmr <- cbind(gsir.l1.rvmr, result_gsir.l1.rvmr)
  gsir.l2.rvmr <- cbind(gsir.l2.rvmr, result_gsir.l2.rvmr)
  gsir.l1.dcor <- cbind(gsir.l1.dcor, result_gsir.l1.dcor)
  gsir.l2.dcor <- cbind(gsir.l2.dcor, result_gsir.l2.dcor)
  gsir.l1.dim <- cbind(gsir.l1.dim, result_gsir.l1.dim)
  gsir.l2.dim <- cbind(gsir.l2.dim, result_gsir.l2.dim)
  gsir.l1.correct <- cbind(gsir.l1.correct, result_gsir.l1.correct)
  gsir.l2.correct <- cbind(gsir.l2.correct, result_gsir.l2.correct)
}

grid.result <- expand.grid(n = n, m=m, complex.x=complex.x, complex.y=complex.y,
                           model = Models, stringsAsFactors=FALSE)

gsir.l1.rvmr_mean <- apply(gsir.l1.rvmr, MARGIN = 1, mean)%>%round(3)
gsir.l1.rvmr_sd <- apply(gsir.l1.rvmr, MARGIN = 1, sd)%>%round(3)
gsir.l2.rvmr_mean <- apply(gsir.l2.rvmr, MARGIN = 1, mean)%>%round(3)
gsir.l2.rvmr_sd <- apply(gsir.l2.rvmr, MARGIN = 1, sd)%>%round(3)

gsir.l1.dcor_mean <- apply(gsir.l1.dcor, MARGIN = 1, mean)%>%round(3)
gsir.l1.dcor_sd <- apply(gsir.l1.dcor, MARGIN = 1, sd)%>%round(3)
gsir.l2.dcor_mean <- apply(gsir.l2.dcor, MARGIN = 1, mean)%>%round(3)
gsir.l2.dcor_sd <- apply(gsir.l2.dcor, MARGIN = 1, sd)%>%round(3)

gsir.l1.correct <- apply(gsir.l1.correct, MARGIN = 1, mean)%>%round(3)
gsir.l2.correct <- apply(gsir.l2.correct, MARGIN = 1, mean)%>%round(3)

gsir <- cbind(grid.result, gsir.l1.rvmr_mean, gsir.l1.rvmr_sd, gsir.l1.dcor_mean, 
               gsir.l1.dcor_sd, gsir.l1.correct, gsir.l2.rvmr_mean, gsir.l2.rvmr_sd, 
               gsir.l2.dcor_mean, gsir.l2.dcor_sd, gsir.l2.correct)

#save(wgsir, file = paste(save_path,"/d2d_wgsir_summary.Rdata"))

################################ write in latex ################################
gsir.l1.rvmr <- gsir.l1.rvmr_mean
gsir.l2.rvmr <- gsir.l2.rvmr_mean
gsir.l1.dcor <- gsir.l1.dcor_mean
gsir.l2.dcor <- gsir.l2.dcor_mean
for (i in 1:nrow(grid.result)) {
  gsir.l1.rvmr[i] <- paste(gsir.l1.rvmr_mean[i],"(",gsir.l1.rvmr_sd[i],")")
  gsir.l2.rvmr[i] <- paste(gsir.l2.rvmr_mean[i],"(",gsir.l2.rvmr_sd[i],")")
  gsir.l1.dcor[i] <- paste(gsir.l1.dcor_mean[i],"(",gsir.l1.dcor_sd[i],")")
  gsir.l2.dcor[i] <- paste(gsir.l2.dcor_mean[i],"(",gsir.l2.dcor_sd[i],")")
}


library(xtable)
xtable(cbind(gsir.l1.rvmr, gsir.l2.rvmr), digits = 3)
xtable(cbind(gsir.l1.dcor, gsir.l2.dcor), digits = 3)

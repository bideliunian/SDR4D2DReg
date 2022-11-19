################################################
####   Simulations for Scenario II: d2v  ########
####    including Model II-1, II-2, II-3   ########

#################### PART 1: Preparation #################################

function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/d2v"
save_path <- "~/work/DR4D2D/d2v/Results"

require(dplyr)

# set global parameters
n_list <- c(200, 400)
p <- 10
m <- 100
Models <- c('d2v_ex1','d2v_ex2','d2v_ex3')
complex.x <- 1
complex.y <- 1

wgsir1.rvmr <- numeric();wgsir2.rvmr <- numeric()
wgsir1.dim <- numeric();wgsir2.dim <- numeric()
wgsir1.dcor <- numeric();wgsir2.dcor <- numeric()
wgsir1.correct <- numeric();wgsir2.correct <- numeric()

############################### read data ##############################
for (i in 1:100) {
  load(file = paste(save_path,"/wgsir.d2v.",i,".RData",sep=""))
  
  result_wgsir1.rvmr <- result_wgsir[1, ]%>%unlist()
  result_wgsir2.rvmr <- result_wgsir[2, ]%>%unlist()
  result_wgsir1.dcor <- result_wgsir[3, ]%>%unlist()
  result_wgsir2.dcor <- result_wgsir[4, ]%>%unlist()
  result_wgsir1.dim <- result_wgsir[5, ]%>%unlist()
  result_wgsir2.dim <- result_wgsir[6, ]%>%unlist()
  result_wgsir1.correct <- result_wgsir[7, ]%>%unlist()
  result_wgsir2.correct <- result_wgsir[8, ]%>%unlist()
  
  wgsir1.rvmr <- cbind(wgsir1.rvmr, result_wgsir1.rvmr)
  wgsir2.rvmr <- cbind(wgsir2.rvmr, result_wgsir2.rvmr)
  wgsir1.dcor <- cbind(wgsir1.dcor, result_wgsir1.dcor)
  wgsir2.dcor <- cbind(wgsir2.dcor, result_wgsir2.dcor)
  wgsir1.dim <- cbind(wgsir1.dim, result_wgsir1.dim)
  wgsir2.dim <- cbind(wgsir2.dim, result_wgsir2.dim)
  wgsir1.correct <- cbind(wgsir1.correct, result_wgsir1.correct)
  wgsir2.correct <- cbind(wgsir2.correct, result_wgsir2.correct)
}

grid.result <- expand.grid(n = n_list, m=m, complex.x=complex.x, complex.y=complex.y,
                           model = Models, stringsAsFactors=FALSE)

wgsir1.rvmr_mean <- apply(wgsir1.rvmr, MARGIN = 1, mean)%>%round(3)
wgsir1.rvmr_sd <- apply(wgsir1.rvmr, MARGIN = 1, sd)%>%round(3)
wgsir2.rvmr_mean <- apply(wgsir2.rvmr, MARGIN = 1, mean)%>%round(3)
wgsir2.rvmr_sd <- apply(wgsir2.rvmr, MARGIN = 1, sd)%>%round(3)

wgsir1.dcor_mean <- apply(wgsir1.dcor, MARGIN = 1, mean)%>%round(3)
wgsir1.dcor_sd <- apply(wgsir1.dcor, MARGIN = 1, sd)%>%round(3)
wgsir2.dcor_mean <- apply(wgsir2.dcor, MARGIN = 1, mean)%>%round(3)
wgsir2.dcor_sd <- apply(wgsir2.dcor, MARGIN = 1, sd)%>%round(3)

wgsir1.correct <- apply(wgsir1.correct, MARGIN = 1, mean)%>%round(3)
wgsir2.correct <- apply(wgsir2.correct, MARGIN = 1, mean)%>%round(3)

wgsir <- cbind(grid.result, wgsir1.rvmr_mean, wgsir1.rvmr_sd, wgsir1.dcor_mean, 
               wgsir1.dcor_sd, wgsir1.correct, wgsir2.rvmr_mean, wgsir2.rvmr_sd, 
               wgsir2.dcor_mean, wgsir2.dcor_sd, wgsir2.correct)

#save(wgsir, file = paste(save_path,"/d2v_wgsir_summary.Rdata"))

################################ write in latex ################################
wgsir1.rvmr <- wgsir1.rvmr_mean
wgsir2.rvmr <- wgsir2.rvmr_mean
wgsir1.dcor <- wgsir1.dcor_mean
wgsir2.dcor <- wgsir2.dcor_mean
for (i in 1:nrow(grid.result)) {
  wgsir1.rvmr[i] <- paste(wgsir1.rvmr_mean[i],"(",wgsir1.rvmr_sd[i],")")
  wgsir2.rvmr[i] <- paste(wgsir2.rvmr_mean[i],"(",wgsir2.rvmr_sd[i],")")
  wgsir1.dcor[i] <- paste(wgsir1.dcor_mean[i],"(",wgsir1.dcor_sd[i],")")
  wgsir2.dcor[i] <- paste(wgsir2.dcor_mean[i],"(",wgsir2.dcor_sd[i],")")
}

library(xtable)
xtable(cbind(wgsir1.rvmr, wgsir1.dcor, wgsir1.correct, wgsir2.rvmr, wgsir2.dcor, wgsir2.correct)
       , digits = 3)

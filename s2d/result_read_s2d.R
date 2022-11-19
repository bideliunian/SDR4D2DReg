#######################################################
###########   Simulations for Scenario I: s2d  ########
##    including Model I-1, I-2, I-3, I-4, I-5   #######

#################### Preparation #################################

function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/s2d"
save_path <- "~/work/DR4D2D/s2d/Results"

require(dplyr)

# set global parameters
n_list <- c(200, 400)
m <- c(50, 100, 200)
Models <- c('s2d_ex1','s2d_ex2','s2d_ex3','s2d_ex4','s2d_ex5')
complex.x <- 1
complex.y <- 1
wgsir1.rvmr <- numeric();wgsir2.rvmr <- numeric()
wgsir1.dim <- numeric();wgsir2.dim <- numeric()
wgsir1.dcor <- numeric();wgsir2.dcor <- numeric()
wgsir1.correct <- numeric();wgsir2.correct <- numeric()

############################### read data ##############################
for (i in 1:100) {
  load(file = paste(save_path,"/wgsir.s2d.",i,".RData",sep=""))
  
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

#save(wgsir, file = paste(save_path,"/s2d_wgsir_summary.Rdata"))

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

m50.index <- which(grid.result$m == 50)
m100.index <- which(grid.result$m == 100)

wgsir1.rvmr.m50 <- wgsir1.rvmr[m50.index]
wgsir2.rvmr.m50 <- wgsir2.rvmr[m50.index]
wgsir1.dcor.m50 <- wgsir1.dcor[m50.index]
wgsir2.dcor.m50 <- wgsir2.dcor[m50.index]

wgsir1.rvmr.m100 <- wgsir1.rvmr[m100.index]
wgsir2.rvmr.m100 <- wgsir2.rvmr[m100.index]
wgsir1.dcor.m100 <- wgsir1.dcor[m100.index]
wgsir2.dcor.m100 <- wgsir2.dcor[m100.index]

library(xtable)
xtable(cbind(rep(c(100, 200), 5), wgsir1.rvmr.m50, wgsir1.rvmr.m100, wgsir2.rvmr.m50, wgsir2.rvmr.m100)
       , digits = 3)

xtable(cbind(rep(c(100, 200), 5), wgsir1.dcor.m50, wgsir1.dcor.m100, wgsir2.dcor.m50, wgsir2.dcor.m100)
       , digits = 3)

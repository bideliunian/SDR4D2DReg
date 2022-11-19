########################################################
####   Simulations to compare GSIR using l1 distance with WGSIR  ########
########################################################

#################### PART 1: Preparation #################################

function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/FGSIRvsWGSIR"
save_path <- "~/work/DR4D2D/FGSIRvsWGSIR/Results"

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"/main_gsir.R",sep=""))

require(dplyr)
require(tictoc)
require(pracma)

#####  global parameters
n <- 200
m <- 100
Models <- c('s2d_ex1','s2d_ex2','s2d_ex3','d2d_ex1','d2d_ex2', 'd2d_ex3', 'd2d_ex4')
complex.x <- 1
complex.y <- 1
order.method <- 'bic'
metric <- 'l1'

############## PART 2: implement the wgsir on s2d models #############

args <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(2021*args)

grid.result <- expand.grid(n = n, m=m, complex.x=complex.x, complex.y=complex.y,
                           model = Models, stringsAsFactors=FALSE)
result_gsir_l1 <- mapply(main_gsir, times = 1, n=grid.result$n, m=grid.result$m,
                       model=grid.result$model, complex.x=grid.result$complex.x,
                       complex.y=grid.result$complex.y, order.method=order.method, 
                       metric=metric, seed = 2021*args)

save(result_gsir_l1, file = paste(save_path, "/l1gsir.",args,".RData",sep=""))


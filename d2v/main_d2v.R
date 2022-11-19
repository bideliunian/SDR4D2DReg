################################################
####   Simulations for Scenario II: d2v  ########
####    including Model II-1, II-2, II-3   ########

#################### PART 1: Preparation #################################

function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/d2v"
save_path <- "~/work/DR4D2D/d2v/Results"

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)

require(dplyr)
require(tictoc)
# require(pracma)

# set global parameters
n_list <- c(200, 400)
p <- 10
m <- 100
Models <- c('d2v_ex1','d2v_ex2','d2v_ex3')
complex.x <- 1
complex.y <- 1
order.method <- 'bic'

############## PART 2: implement the wgsir on d2v models #####################

args <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(2021*args)

grid.result <- expand.grid(n=n_list, m=m, p=p, complex.x=complex.x, complex.y=complex.y,
                           model=Models, stringsAsFactors=FALSE)
result_wgsir <- mapply(main_wgsir, times = 1, n=grid.result$n, m=grid.result$m,
                       p=grid.result$p, model=grid.result$model,
                       complex.x=grid.result$complex.x, complex.y=grid.result$complex.y,
                       order.method=order.method)

save(result_wgsir, file = paste(save_path, "/wgsir.d2v.",args,".RData",sep=""))


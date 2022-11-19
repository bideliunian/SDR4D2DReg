#######################################################
####   Order simulations for Scenario III: d2d  #######
####    including Model III-(1-6)   ###################

#################### PART 1: Preparation #################################

function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/d2d"
save_path <- "~/work/DR4D2D/d2d/Results"

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)

require(dplyr)
require(tictoc)
# require(pracma)

# set global parameters
n_list <- c(100) #, 200
m <- 100
Models <- c('d2d_ex1','d2d_ex2','d2d_ex3','d2d_ex4','d2d_ex5','d2d_ex6')
complex_x <- 1
complex_y <- 1
order_method <- 'bic'

############## PART 2: order determination on d2d models #############

args <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(2021*args)

grid.result <- expand.grid(n = n_list, m=m, complex.x=complex_x, complex.y=complex_y,
                           model = Models, stringsAsFactors=FALSE)
result_order <- mapply(main_order, times = 1, n=grid.result$n, m=grid.result$m, p=10,
                       model=grid.result$model,complex.x=grid.result$complex.x,
                       complex.y=grid.result$complex.y, order.method=order_method)

save(result_order, file = paste(save_path, "/wgsir.d2d.",args,".RData",sep=""))


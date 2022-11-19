################################################
####   Simulations for Scenario IV: md2md  ########
####    including Model IV-1, IV-2, IV-3   ########

#################### PART 1: Preparation #################################

function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/md2md"
save_path <- "~/work/DR4D2D/md2md/Results"

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)

require(dplyr)
require(tictoc)
require(MASS)
# require(pracma)

# set global parameters
n_list <- c(200, 400)
m <- c(50, 100)
Models <- c('md2md_ex1','md2md_ex2','md2md_ex3','md2md_ex4','md2md_ex5')
complex.x <- 1
complex.y <- 1
order.method <- 'bic'


############## PART 2: implement the wgsir on d2v models #############

args <- as.numeric(commandArgs(trailingOnly=TRUE))

grid.result <- expand.grid(n = n_list, m=m, complex.x=complex.x, complex.y=complex.y,
                           model = Models, stringsAsFactors=FALSE)
result_wgsir <- mapply(main_wgsir, times = 1, n=grid.result$n, m=grid.result$m,
                       model=grid.result$model, complex.x=grid.result$complex.x,
                       complex.y=grid.result$complex.y, order.method=order.method,
                       seed = 2021*args)

save(result_wgsir, file = paste(save_path, "/wgsir.md2md.",args,".RData",sep=""))


###################################################
########## PART 1: Preparation #############################
##################################################

function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/CalgaryTempData"
save_path <- "~/work/DR4D2D/FGSIRvsWGSIR/CalgaryTempData"

function_path <- "D:/Research/DR4D2D/Codes/Functions"
working_path <- "D:/Research/DR4D2D/Codes/CalgaryTempData"
save_path <- "D:/Research/DR4D2D/Codes/CalgaryTempData"


library(pracma)
library(dplyr)
library(ggplot2)
library(plotly)
library(RColorBrewer)

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)

##########################################################
################# PART 2: Processing Data ##############
#######################################################
spring <- read.csv(file = paste(working_path, '/spring_Calgary.csv', sep=''))
summer <- read.csv(file = paste(working_path, '/summer_Calgary.csv', sep=''))
fall <- read.csv(file = paste(working_path, '/fall_Calgary.csv', sep=''))
winter <- read.csv(file = paste(working_path, '/winter_Calgary.csv', sep=''))

spring <- na.omit(spring)
summer <- na.omit(summer)
fall <- na.omit(fall)
winter <- na.omit(winter)


### using summer to predict fall
## predictor is 136 * 92 * 2
n = dim(summer)[2]/2
m_spring = dim(spring)[1]
m_summer = dim(summer)[1]
m_fall = dim(fall)[1]
m_winter = dim(fall)[1] 
d = 2
spring_arr = array(rep(NA, n*m_spring*2), dim = c(n,m_spring,d))
summer_arr = array(rep(NA, n*m_summer*2), dim = c(n,m_summer,d))
fall_arr = array(rep(NA, n*m_fall*2), dim = c(n,m_fall,d))
winter_arr = array(rep(NA, n*m_winter*2), dim = c(n,m_winter,d))
for(i in 1:n){
  spring_arr[i,,] <- as.matrix(spring[,(2*i-1):(2*i)])
  summer_arr[i,,] <- as.matrix(summer[,(2*i-1):(2*i)])
  fall_arr[i,,] <- as.matrix(fall[,(2*i-1):(2*i)])
  winter_arr[i,,] <- as.matrix(winter[,(2*i-1):(2*i)])
}


########################################
##### PART 3: training #################
###################################


## paras
complex.x <- 1
complex.y <- 1


## spring vs summer
K_spring <- gram_sw(spring_arr, complexity = complex.x, n_seed = 2021)
K_summer <- gram_sw(summer_arr, complexity = complex.y, n_seed = 2021)
# K_fall <- gram_sw(fall_arr, complexity = complex.x, n_seed = 2021)
# K_winter <- gram_sw(winter_arr, complexity = complex.x, n_seed = 2021)


eps.grid <- 10^(-6:0)
epsx <- lapply(eps.grid, FUN = gcv, x=K_spring, y=K_summer, which='ex', ytype='gram', xtype='gram',
               complex_x=complex.x,complex_y=complex.y)
epsx.min <- eps.grid[which.min(epsx)]

epsy <- lapply(eps.grid, FUN = gcv, x=K_spring, y=K_summer, which='ey', ytype='gram', xtype='gram',
               complex_x=complex.x,complex_y=complex.y)
epsy.min <- eps.grid[which.min(epsy)]

r0_iden <- bic(x=K_spring, y=K_summer, xtype='gram', ytype='gram', atype='identity',
               ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
               criterion='lal', c0=2)$rhat
r0_yinv <- bic(x=K_spring, y=K_summer, xtype='gram', ytype='gram', atype='Gyinv',
               ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
               criterion='lal', c0=4)$rhat  

estpred.wgsir1 <- wgsir(x=K_spring, x_new=K_spring, y=K_summer, xtype='gram', ytype='gram', atype='identity',
                        ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
                        r=r0_iden)$pred.est
estpred.wgsir2 <- wgsir(x=K_spring, x_new=K_spring, y=K_summer, xtype='gram', ytype='gram', atype='Gyinv',
                        ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y,
                        r=r0_iden)$pred.est


#################################################
###### PART 4: generating figures ############
####################################################
x_cord = y_cord = c()
for (i in 1:n) {
  x_cord = c(x_cord, summer_arr[i,,1])
  y_cord = c(y_cord, summer_arr[i,,2])
}

plotdata<- data.frame(
  x = x_cord,
  y = y_cord,
  estpred1 = rep(estpred.wgsir1[,1], each = m_summer),
  estpred2 = rep(estpred.wgsir2[,1], each = m_summer)
)


plot_data_est = function (data, i) {
  plotdata_md_ex = data
  ggplot(plotdata[dense_rank(plotdata$estpred2)==(13*i),], aes(x=x, y=y)) +  geom_point()+
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(density)),
      contour = FALSE
    ) + 
    scale_fill_viridis_c() + 
    labs(x = "min temperature", y = 'range' ) + 
    scale_fill_distiller(palette= "Spectral") + 
    xlim(-10, 20)+ylim(0, 40) +
    theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

myplots_est <- lapply(1:10, plot_data_est, data = plotdata)

library(ggpubr) 
figure1 <- ggarrange(myplots_est[[1]], myplots_est[[2]], myplots_est[[3]], myplots_est[[4]], myplots_est[[5]],
                     myplots_est[[6]], myplots_est[[7]], myplots_est[[8]], myplots_est[[9]], myplots_est[[10]],
                     ncol = 10)
figure1 <- ggarrange(myplots_est[[1]], myplots_est[[3]], myplots_est[[5]], myplots_est[[7]], myplots_est[[9]],
                     ncol = 5)


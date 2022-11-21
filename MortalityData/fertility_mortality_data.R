###################################################
########## PART 1: Preparation #############################
##################################################
function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/MortalityData"
save_path <- "~/work/DR4D2D/FGSIRvsWGSIR/MortalityData"

library(frechet)
#library(pracma)
library("readxl")
#library(dplyr)
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
# read file
my_data <- list()

deaths_age <- read_excel(paste(working_path, "WPP2019_MORT_F17_1_ABRIDGED_LIFE_TABLE_BOTH_SEXES.xlsx", sep = '/'),
                         sheet = "ESTIMATES", na = '...')%>%
  filter(Type=='Country/Area'&Period=='2015-2020'&`Age (x)`!=0)%>%
  select(`Region, subregion, country or area *`,`Country code`,`Age (x)`,`Number of deaths d(x,n)`)%>%
  rename(country.or.area=`Region, subregion, country or area *`, Country.code=`Country code`)
birth_age <- read_excel(paste(working_path, "WPP2019_FERT_F06_BIRTHS_BY_AGE_OF_MOTHER.xlsx", sep = '/'),
                        sheet = "ESTIMATES", na = '...')%>%
   filter(Type=='Country/Area'&`Period`=='2015-2020')%>%
   select(`Region, subregion, country or area *`,`Country code`,`15-19`:`45-49`)%>%
   rename(country.or.area=`Region, subregion, country or area *`, Country.code=`Country code`)

my_data <- list(deaths_age, birth_age)
var_names <- c('deaths_age','birth_age')
names(my_data) <- var_names
country_code <- lapply(my_data, `[[`, "Country.code")
common_country <- Reduce(intersect, country_code)
common_country <- sort(common_country)
my_data_common <- lapply(my_data, filter, Country.code%in%common_country)%>% lapply(arrange,Country.code)

############# create death density #####################
n <- length(common_country)
m <- nrow(my_data_common$deaths_age)/n
y <- as.matrix(my_data_common$deaths_age[["Age (x)"]][1:m])
for (i in common_country) {
  y <- cbind(y,filter(my_data_common$deaths_age, Country.code==i)[["Number of deaths d(x,n)"]])
}
rownames(y) <- y[,1]
y <- y[,-1]
ylist <- split(y,col(y)) 
colnames(y) <- common_country
## using frechet packege
x0 <-seq(0,100,length.out=101)#outputGrid
biny <- c(0,1,seq(5,100, by=5))
death_density <- lapply(ylist, FUN = CreateDensity, bin=biny, optns=list(outputGrid=x0),y = NULL, histogram = NULL)
death_density <- lapply(death_density, function(x) x <- x[-1])

########### create birth_age density ############
m_birth <- ncol(my_data_common$birth_age)-2
x <- t(as.data.frame(my_data_common$birth_age[,-c(1:2)]))
colnames(x) <- common_country
xlist <- split(x,col(x)) 
## using frechet packege
x0_birth <-seq(15,50,length.out=71)#outputGrid
binx <- c(seq(15,50, by=5))
birth_density <- lapply(xlist, FUN = CreateDensity, bin=binx, optns=list(outputGrid=x0_birth),y = NULL, histogram = NULL)
birth_density <- lapply(birth_density, function(x) x <- x[-1])

############# save data ################
save(death_density, file = paste(save_path, "/mort_density.RData", sep=""))
save(birth_density, file = paste(save_path, "/fert_density.RData", sep=""))

############ plot the densities ######################
p_death <- ggplot()
plot_death <- function(l){
  df <- data.frame(x=l$x,y=l$y)
  p_death <<- p_death + geom_line(data=df, aes(x,y),alpha=0.1)
}
lapply(death_density, plot_death)
p_death+labs(title="",x="Age(0-100)", y = "Density of age at death") + theme_bw()
p_birth <- ggplot()
plot_birth <- function(l){
  df <- data.frame(x=l$x,y=l$y)
  p_birth <<- p_birth + geom_line(data=df, aes(x,y),alpha=0.1)
}
lapply(birth_density, plot_birth)
p_birth+labs(title="",x="Age(15-49)", y = "Density of mother's age at birth") + theme_bw()



###########################################################
############ PART 3: SDR using WGSIR ########################
###########################################################

## paras

complex.x <- 1
complex.y <- 1
## read data ##############
load(paste(working_path, 'mort_density.RData', sep="/"))
load(paste(working_path, 'fert_density.RData', sep="/"))


gram_wass_smooth <- function(x, complexity){
  # Param
  #   x is a list of densities
  #   complexity controls the tuning parameter in kernel
  # Return
  #   n by n gram matrix
  n <- length(x)
  
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  
  dist <- function(i,k){
    return(frechet::dist4den(d1=x[[i]],d2=x[[k]]))
  }
  
  kupper <- mapply(dist, i=ind[,1], k=ind[,2])
  k <- diag(n)
  k[upper.tri(k,diag <- FALSE)] <- kupper^2
  k <- k+t(k)-2*diag(n)
  sigma2 <- sum(kupper^2)/choose(n,2)
  gamma <- complexity/(2*sigma2)
  
  return(exp(-gamma*k))
}


Kx <- gram_wass_smooth(birth_density, complexity = complex.x)
Ky <- gram_wass_smooth(death_density, complexity = complex.y)

##### if exculde outliers #########
outlier=c(121,97,84,101,62,148,185)
birth_density <- birth_density[-outlier]
death_density <- death_density[-outlier]
Kx <- gram_wass_smooth(birth_density, complexity = complex.x)
Ky <- gram_wass_smooth(death_density, complexity = complex.y)


############# wgsir #######################

eps.grid <- 10^(-6:0)
epsx <- lapply(eps.grid, FUN = gcv, x=Kx, y=Ky, which='ex', ytype='gram', xtype='gram',
               complex_x=complex.x,complex_y=complex.y)
epsx.min <- eps.grid[which.min(epsx)]

epsy <- lapply(eps.grid, FUN = gcv, x=Kx, y=Ky, which='ey', ytype='gram', xtype='gram',
               complex_x=complex.x,complex_y=complex.y)
epsy.min <- eps.grid[which.min(epsy)]

r0_iden <- bic(x=Kx, y=Ky, xtype='gram', ytype='gram', atype='identity',
                 ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
                 criterion='lal', c0=2)$rhat
r0_yinv <- bic(x=Kx, y=Ky, xtype='gram', ytype='gram', atype='Gyinv',
                 ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
                 criterion='lal', c0=4)$rhat  

estpred.wgsir1 <- wgsir(x=Kx, x_new=Kx, y=Ky, xtype='gram', ytype='gram', atype='identity',
                      ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y, 
                      r=r0_iden)$pred.est
estpred.wgsir2 <- wgsir(x=Kx, x_new=Kx, y=Ky, xtype='gram', ytype='gram', atype='Gyinv',
                      ex=epsx.min, ey=epsy.min, complex_x=complex.x, complex_y=complex.y,
                      r=r0_iden)$pred.est
n <- length(estpred.wgsir1)
nSup <- length(death_density[[1]]$x)

########################################################
############### PART 4: generating plots #################
############################################################
palette <- colorRampPalette(c('#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                              '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))
plotdata <- data.frame(
  x = unlist(lapply(death_density, "[[","x")),
  z = unlist(lapply(death_density, "[[", "y")),
  estpred = rep(estpred.wgsir1, each = nSup),
  estpred2 = rep(estpred.wgsir2, each = nSup),
  y = rep(runif(n,0,1), each = nSup))

################### response vs estimated pred ###################

## wgsir1
fig_estpred <- plot_ly(plotdata, x = ~x, y = ~estpred, z = ~z, type = 'scatter3d', 
                       mode = 'lines', color = ~z, split = ~estpred2, alpha=1, 
                       colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "Age-at-death"),
      yaxis = list(title = "1st WGSIR1 predictor"),
      zaxis = list(title = "density")
    ),showlegend = FALSE)
fig_estpred%>% hide_colorbar()

## wgsir2
fig_estpred2 <- plot_ly(plotdata, x = ~x, y = ~estpred2, z = ~z, type = 'scatter3d', 
                        mode = 'lines', color = ~z, split = ~estpred2, alpha=1, 
                        colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "Age-at-death"),
      yaxis = list(title = "1st WGSIR2 predictor"),
      zaxis = list(title = "density")
    ),showlegend = FALSE)
fig_estpred2%>% hide_colorbar()

## chaos
fig_chaos <- plot_ly(plotdata, x = ~x, y = ~y, z = ~z, type = 'scatter3d',
                     mode = 'lines', color = ~z, split = ~y, alpha=1,colors = palette(100))%>% 
  layout(scene = list(
    xaxis = list(title = "Age-at-death"),
    yaxis = list(title = ""),
    zaxis = list(title = "density")),showlegend = FALSE)
fig_chaos%>% hide_colorbar()

################ summary stat plot ########################
plotdata_sum <- data.frame()
plotdata_sum <- data.frame(
  wgsir1_estpred = estpred.wgsir1,
  wgsir2_estpred = estpred.wgsir2,
  mean = sapply(X = death_density,FUN = function(l) sum(l[["y"]]*l[["x"]]))%>%as.vector(),
  mode = sapply(X = lapply(death_density,"[[","y"), FUN = function(l) (which.max(l)-1))%>%as.vector(),
  mom2 = sapply(X = death_density, FUN = function(l) sum(l[["y"]]*(l[["x"]])^2)),
  mom3 = sapply(X = death_density, FUN = function(l) sum(l[["y"]]*(l[["x"]])^3))
)

plotdata_sum$var <- plotdata_sum$mom2-plotdata_sum$mean^2
plotdata_sum$sd <- sqrt(plotdata_sum$var)
plotdata_sum$skew <- (plotdata_sum$mom3-3*plotdata_sum$mean*plotdata_sum$var-plotdata_sum$mean^3)/sqrt(plotdata_sum$var^{3})

fig_mean <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir1_estpred, y =plotdata_sum$mean, 
                    name = 'Mean', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR1 predictor"))
fig_mode <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir1_estpred, y =plotdata_sum$mode, 
                    name = 'Mode', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR1 predictor"))
fig_var <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir1_estpred, y =plotdata_sum$var,  
                   name = 'Variance', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR1 predictor"))
fig_skew <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir1_estpred, y =plotdata_sum$skew,  
                    name = 'Skewness', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR1 predictor"))
fig_sd <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir1_estpred, y =plotdata_sum$sd,  
                  name = 'Sd', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR1 predictor"))
fig_sum <- subplot(fig_mean, fig_mode, fig_sd, fig_skew, nrows = 2, shareX = TRUE, titleX = TRUE)
fig_sum

### wgsir2
fig_mean2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir2_estpred, y =plotdata_sum$mean, 
                     name = 'Mean', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR2 predictor"))
fig_mode2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir2_estpred, y =plotdata_sum$mode, 
                     name = 'Mode', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR2 predictor"))
fig_var2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir2_estpred, y =plotdata_sum$var,  
                    name = 'Variance', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR2 predictor"))
fig_skew2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir2_estpred, y =plotdata_sum$skew,  
                     name = 'Skewness', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR2 predictor"))
fig_sd2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$wgsir2_estpred, y =plotdata_sum$sd,  
                   name = 'Sd', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = "1st WGSIR2 predictor"))
fig_sum2 <- subplot(fig_mean2,fig_mode2,fig_sd2,fig_skew2,nrows = 2, shareX = TRUE, titleX = TRUE)
fig_sum2



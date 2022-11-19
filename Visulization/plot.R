###################################################
########## PART 1: Preparation #############################
##################################################
function_path <- "~/work/DR4D2D/Functions"
working_path <- "~/work/DR4D2D/Visulization"
save_path <- "~/work/DR4D2D/FGSIRvsWGSIR/Visulization"


function_path <- "D:/Research/DR4D2D/Codes/Functions"
working_path <- "D:/Research/DR4D2D/Codes/Visulization"
save_path <- "D:/Research/DR4D2D/Codes/Visulization"


# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path, '/plot_function.R', sep=''))

require(dplyr)
# require(pracma)
require(MASS)
require(ggplot2)
require(plotly)
require(RColorBrewer)
library(ggpubr)

########################################################
## PART 2: visulizing  #####################
######################################################

palette <- colorRampPalette(c('#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                           '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))

################################
##### d2d ex1
#################################

result_ex1 <- plot_wgsir(n=400, m=100, model='d2d_ex1',complex.x=1,complex.y=1)
nSup <- length(result_ex1$density[[1]]$x)
plotdata_ex1 <- data.frame(
  x = unlist(lapply(result_ex1$density, "[[","x")),
  z = unlist(lapply(result_ex1$density, "[[", "y")),
  estpred = rep(result_ex1$estpred.sir1, each = nSup),
  truepred = rep(result_ex1$truepred, each = nSup))

### response vs estpred
fig_ex1_estpred <- plot_ly(plotdata_ex1, x = ~x, y = ~estpred, z = ~z, type = 'scatter3d', 
                        mode = 'lines', color = ~z, split = ~estpred, alpha=1, colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "support",range = c(-3, 7)),
      yaxis = list(title = "1st WGSIR predictor"),
      zaxis = list(title = "density")
    ),showlegend = FALSE)

fig_ex1_estpred%>% hide_colorbar()

### response vs truepred
fig_ex1_truepred <- plot_ly(plotdata_ex1, x = ~x, y = ~truepred, z = ~z, type = 'scatter3d', 
                         mode = 'lines', color = ~z, split = ~truepred, alpha=1, colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "support",range = c(-3, 7)),
      yaxis = list(title = "True predictor"),
      zaxis = list(title = "density")
    ),showlegend = FALSE)

fig_ex1_truepred%>% hide_colorbar()

### scatter plot estpred vs truepred
pred_ex1 <- data.frame(truepred=result_ex1$truepred,
                   estpred=result_ex1$estpred.sir1)
ggplot(pred_ex1, aes(x=truepred, y=estpred)) + geom_point() + 
  labs(y = "Sufficient Predictor", x = "True predictor") +  theme(text = element_text(size = 15)) + theme_bw() 



########################################
###### d2d ex2
########################################

set.seed(1)
result_ex2 <- plot_wgsir(n=400, m=100, model='d2d_ex2',complex.x=1,complex.y=1)
nSup <- length(result_ex2$density[[1]]$x)
plotdata_ex2 <- data.frame(
  x = unlist(lapply(result_ex2$density, "[[","x")),
  z = unlist(lapply(result_ex2$density, "[[", "y")),
  estpred1 = rep(result_ex2$estpred.sir1[,1], each = nSup),
  estpred2 = rep(result_ex2$estpred.sir1[,2], each = nSup),
  truepred1 = rep(result_ex2$truepred[,1], each = nSup),
  truepred2 = rep(result_ex2$truepred[,2], each = nSup)
  )

### response vs estpred
fig_estpred1 <- plot_ly(plotdata_ex2, x = ~x, y = ~estpred1, z = ~z, type = 'scatter3d', 
                   mode = 'lines', color = ~z, split = ~estpred1, alpha=1, colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "support",range = c(-4, 6)),
      yaxis = list(title = "1st WGSIR predictor"),
      zaxis = list(title = "density")
    ),showlegend = FALSE)
fig_estpred2 <- plot_ly(plotdata_ex2, x = ~x, y = ~estpred2, z = ~z, type = 'scatter3d', 
                        mode = 'lines', color = ~z, split = ~estpred2, alpha=1, colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "support",range = c(-4, 6)),
      yaxis = list(title = "2nd WGSIR predictor"),
      zaxis = list(title = "density")
    ),showlegend = FALSE)

fig_estpred1

### response vs truepred
fig_truepred1 <- plot_ly(plotdata, x = ~x, y = ~truepred1, z = ~z, type = 'scatter3d', 
                       mode = 'lines', color = ~z, split = ~truepred1, alpha=1, colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "support",range = c(-4, 6)),
      yaxis = list(title = "1st true predictor"),
      zaxis = list(title = "density")
    ),showlegend = FALSE)
fig_truepred2 <- plot_ly(plotdata, x = ~x, y = ~truepred2, z = ~z, type = 'scatter3d', 
                         mode = 'lines', color = ~z, split = ~truepred2, alpha=1, colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "support",range = c(-4, 6)),
      yaxis = list(title = "2nd true predictor"),
      zaxis = list(title = "density")
    ),showlegend = FALSE)
fig_truepred1
# fig_truepred2

### scatter plot estpred vs truepred
pred <- data.frame(truepred1=result_ex2$truepred[,1],truepred2=result_ex2$truepred[,2],
                   estpred1=result_ex2$estpred.sir1[,1],estpred2=result_ex2$estpred.sir1[,2])

ggplot(pred, aes(x=estpred1, y=estpred2)) + geom_point(aes(colour=truepred1)) + 
  scale_color_gradient(low="grey", high="red") + labs(y = "2nd Sufficient Predictor", x = "1st Sufficient predictor") + 
  theme(text = element_text(size = 15)) + theme_bw() 
ggplot(pred, aes(x=estpred1, y=estpred2)) + geom_point(aes(colour=truepred2)) + 
  scale_color_gradient(low="grey", high="blue") + labs(y = "2nd Sufficient Predictor", x = "1st Sufficient predictor") + 
  theme(text = element_text(size = 15)) + theme_bw() 


#############################
## multidim ex2
#################################

set.seed(2021)
result_md <- plot_wgsir(n=200, m=100, model='md2md_ex2',complex.x=0.1,complex.y=0.1)
x_cord = y_cord = c()
nSup = dim(result_md$estpred.sir2)[1]

for (i in 1:n) {
  x_cord = c(x_cord, result_md[["ytest"]][[i]][,1])
  y_cord = c(y_cord, result_md[["ytest"]][[i]][,2])
}
plotdata_md <- data.frame(
  x = x_cord,
  y = y_cord,
  estpred1 = rep(result_md$estpred.sir2[,1], each = nSup),
  estpred2 = rep(result_md$estpred.sir2[,2], each = nSup),
  truepred1 = rep(result_md$truepred[,1], each = nSup),
  truepred2 = rep(result_md$truepred[,2], each = nSup)
)


plot_data_est = function (data, i) {
  ggplot(data[dense_rank(data$estpred1)==(10*i),], aes(x=x, y=y)) +
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(density)),
      contour = FALSE
    ) + scale_fill_viridis_c() + scale_fill_distiller(palette= "Spectral")+
    labs(x = "y1", y = 'y2' )+# , subtitle = paste0(10*i, '%')
    xlim(-5, 8) +ylim(-5,8) + theme_bw() +
    theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
}


plot_data_true = function (data, i) {
  ggplot(data[dense_rank(data$truepred1)==(10*i),], aes(x=x, y=y)) +
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(density)),
      contour = FALSE
    ) + scale_fill_viridis_c()+ scale_fill_distiller(palette= "Spectral")+
    labs(x = "y1", y = 'y2' )+#subtitle = paste0(10*i, '%')
    xlim(-5, 8) +ylim(-5, 8) + theme_bw() + 
    theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

myplots_est <- lapply(1:10, plot_data_est, data = plotdata_md_ex2)
myplots_true <- lapply(1:10, plot_data_true, data = plotdata_md_ex2)


figure2 <- ggarrange(myplots_est[[1]], myplots_est[[3]], myplots_est[[5]], myplots_est[[7]],myplots_est[[9]],
                    myplots_true[[1]], myplots_true[[3]], myplots_true[[5]], myplots_true[[7]],myplots_true[[9]],
                    ncol = 5, nrow = 2)
figure2



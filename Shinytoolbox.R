#################################################-
## Toolbox 
#################################################-
## Preliminaries ----
#################################################-
library(shiny)
#library(rgl)
library(scatterplot3d)
library(pracma)
library(mvtnorm)
library(ggplot2)
library(igraph)
library(scatterplot3d)
library(deSolve)
library(latex2exp)
#################################################-
## Functions: calculate metrics ----
## Reference paper
#################################################-
mu <- c(1,1)
A <- matrix(c(1, 0.2, 0.2, 1), 2,2)

pars <- list(mu = mu, A = A)
pars$sign_A11 <- sign(pars$A[1,1])

compute_NFD <- function(pars){
  # computes all parameters for plotting
  
  # equilibrium densities
  pars$N_star <- 1/diag(pars$A)
  pars$c_ij <- sqrt(abs(pars$A[2,2]*pars$A[1,2]/
                          pars$A[2,1]/pars$A[1,1]))
  
  pars$mu_i <- pars$mu[1]
  pars$r_i <- pars$mu_i - pars$A[1,2]*pars$N_star[2]
  pars$eta_i <- pars$mu_i - pars$A[1,1]*pars$c_ij*pars$N_star[2]
  
  pars$ND <- round((pars$r_i - pars$eta_i)/(pars$mu_i - pars$eta_i),3)
  pars$FD <- round((0 - pars$eta_i)/(pars$mu_i - pars$eta_i),3)
  
  pars$N_j_dagger <- pars$mu_i/pars$A[1,2]
  pars$N_i_star <- pars$mu_i/pars$A[1,1]
  
  # compute plotting densities
  pars$N_i_star_plot <- pars$N_i_star/(pars$c_ij*pars$N_star[2])
  pars$N_j_dagger_plot <- pars$N_j_dagger/pars$N_star[2]
  pars$c_ij_N_j_plot <- 1
  pars$N_j_star_plot <- 1
  
  pars$N_i_hash <- pars$N_i_star_plot*(pars$N_j_star_plot-pars$N_j_dagger_plot)/(0-pars$N_j_dagger_plot)
  
  return(pars)
}

compute_LV_from_dens <- function(N_i_star, N_j_dagger, pars){
  # if densities are changed computed LV parameters and then recomputes NFD
  pars$A[1,1] <- pars$sign_A11*abs(pars$mu_i/N_i_star)
  pars$A[1,2] <- pars$mu_i/N_j_dagger
  pars$mu[1] <- sign(N_i_star/pars$A[1,1])
  
  # re-compute niche and fitness differences
  pars <- compute_NFD(pars)
  return(pars)
}

compute_LV_from_NFD <- function(ND, FD, pars){
  # if ND or FD are changed computes LV parameters
  abs_a12 <- abs(pars$mu[2]/pars$A[2,2]*(1-ND)/(1-FD))
  abs_a11 <- abs(pars$A[2,1]*pars$A[2,2]/pars$mu[2]/(1-ND)/(1-FD))
  
  # check whether information aligns with intraspecific facilitation
  pars$A[1,1] <- round(abs_a11*pars$sign_A11,3)
  pars$A[1,2] <- round(abs_a12*sign((1-ND)*pars$A[1,1]),3)
  pars$mu[1] <- sign(pars$A[1,1]*(1-FD))
  
  pars <- compute_NFD(pars)
  return(pars)
}


pars <- compute_NFD(pars)
ND <- 0.8
FD <- 0.9
pars <- compute_LV_from_NFD(ND, FD, pars)

print(c(pars$ND, ND))
print(c(pars$FD, FD))

print(pars$A)
print(pars$mu)

pars <- compute_LV_from_dens(0.6, 0.8, pars)
print(pars$N_i_star)
print(pars$N_j_dagger)

#################################################-
## Functions: make plots ----
#################################################-

mk_graph_2sp <- function(alphamat, rs){
  alphamat.pos <- unname(abs(alphamat))
  g <- igraph::graph_from_adjacency_matrix(alphamat.pos  > 0)
  igraph::E(g)$weight <- as.numeric(alphamat.pos)
  widths <- igraph::E(g)$weight *4
  color.vertex <- c()
  if(rs[1] >= 0){
    color.vertex[1] <- "grey"}else{color.vertex[1] <- "red"}
  if(rs[2] >= 0){
    color.vertex[2] <- "grey"}else{color.vertex[2] <- "red"}
  
  color.mat <- as.numeric(alphamat)
  color.mat[which(color.mat > 0)] <- 1
  color.mat[which(color.mat < 0)] <- 3
  igraph::E(g)$lty <- color.mat
  #widths[widths > 1] <- sqrt(widths)
  plot(g,
       margin = c(0.3,0 ,0, 0),
       #xlim = c(-1.25, 1.25), 
       ylim = c(-1.2, -0.5),
       vertex.label = c("species i","species j"),
       vertex.label.family="Helvetica",   
       vertex.label.cex = 2,
       vertex.label.dist= -2,
       vertex.label.degree = c(pi/2,pi/2),
       vertex.label.color = "black",
       vertex.size = 10 * abs(rs),
       vertex.color = color.vertex,
       vertex.frame.color = "transparent",
       edge.curved = TRUE,
       edge.width = 2*widths,
       edge.arrow.size = 2,
       edge.arrow.mode = c(2, 1,
                           1, 2),
       edge.color = "black",
       edge.loop.angle = c(pi/2,0,0,pi/2),
       layout = matrix(c(-1,0.5, 1, 0.5), 
                       ncol = 2,
                       byrow = TRUE))
}

#mk_graph_2sp(A,mu)

library(purrr)
library(pracma) 
graph_growth_plane <- function(pars){
  normal_vec = pracma::cross(c(-pars$N_i_star_plot,0, pars$mu_i),
                             c(0, -pars$N_j_dagger_plot, pars$mu_i)) # I am not sure this does what you want
  x <- seq(0, 1.2*max(c(pars$N_i_star_plot,
                        pars$N_j_dagger_plot,
                        pars$N_j_star_plot,
                        pars$c_ij_N_j_plot)), length.out = 21)
  pars$y <- rep(x, length(x))
  pars$x <- rep(x, each = length(x))
  offset = normal_vec[3]*pars$mu_i
  pars$z = (offset - normal_vec[1]*pars$x - normal_vec[2]*pars$y)/normal_vec[3]
  return(pars)
}

graph_NFD_plan <- function(pars){
  
  # Generate some point data
  
  plot.df <- data.frame(point =c("Fi"),x=0,
                        y=0,z=0,color.pt="darkgreen")
  plot.df[2,] <-c("Ni",0,0,pars$r_i,color.pt="darkgreen")
  plot.df[3,] <- c("0",0,0,pars$eta_i,color.pt="darkgreen")
  plot.df[4,] <- c("Ni-Fi",pars$N_i_hash,pars$N_j_star_plot,
                    0,color.pt="darkgreen")
  #plot.df[5,] <- c("1",NA,0,pars$mu_i,color.pt="darkgreen")
  plot.df[5,] <- c("Ni#",pars$N_i_hash,0,0,color.pt="red")
  plot.df[6,] <- c("Ni*",pars$N_i_star_plot,0,0,color.pt="red")
  plot.df[7,] <- c("Nj*",0,pars$N_j_star_plot,0,color.pt="red")
  plot.df[8,] <- c("Njdag",0,pars$N_j_dagger_plot,0,color.pt="red")
  plot.df[9,] <- c("CijNj*",pars$c_ij_N_j_plot,0,0,color.pt="red")
  plot.df[10,] <- c("r_i",0,pars$N_j_star_plot,pars$r_i,color.pt="blue")
  plot.df[11,] <- c("mu_i-1",0,0,pars$mu_i,color.pt="blue")
  plot.df[12,] <- c("eta_i",pars$c_ij_N_j_plot,0,pars$eta_i,color.pt="blue")
  
  
  # Ensure columns of plane_coords are numeric
  plot.df[,c(2,3,4)] <- as.data.frame(lapply(plot.df[,c(2,3,4)], as.numeric))
  m <- list(
    l = 0,
    r = 0,
    b = 0,
    t = 0
  )
  # generate lines data
  # z axis
  line1 <- data.frame(x = 0,
                      y = 0,
                      z = c(min(plot.df$z)-1,max(plot.df$z)+1))
  # x axis
  line2 <- data.frame(x = c(min(plot.df$x)-1,max(plot.df$x)+1),
                      y = 0,
                      z = 0)
  # y axis
  line3 <- data.frame(x = 0,
                      y = c(min(plot.df$y)-1,max(plot.df$y)+1),
                      z = 0)
  # eta_i -> 0
  orline1 <- data.frame(x = c(pars$c_ij_N_j_plot, pars$c_ij_N_j_plot, 0),
                        y = c(0,0,0),
                        z = c(0, pars$eta_i, pars$eta_i))
  # r_i -> ND
  orline2 <- data.frame(x = c(0,0,0),
                        y = c(pars$N_j_star_plot, pars$N_j_star_plot, 0),
                        z = c(0, pars$r_i, pars$r_i))
  # N_j_star -> N_i_hash
  orline3 <- data.frame(x = c(pars$N_i_hash,pars$N_i_hash,0),
                        y = c(0, pars$N_j_star_plot, pars$N_j_star_plot),
                        z = c(0, 0, 0))
  
  
  blueline <- data.frame(x = c(pars$N_i_star_plot,0, pars$N_i_hash),
                         y = c(0,pars$N_j_dagger_plot, pars$N_j_star_plot),
                         z = c(0,0,0)) # blue line 
  
  # generate plane data
  plane.data <- data.frame(x = graph_growth_plane(pars)$x,
                           y = graph_growth_plane(pars)$y,
                           z = graph_growth_plane(pars)$z) %>%
    dplyr::filter(z < max(plot.df$z)+0.1) %>%
    dplyr::filter(y < max(plot.df$y)+0.1) %>%
    dplyr::filter(x < max(plot.df$x)+0.1) %>%
    dplyr::filter(z > min(plot.df$z)-0.5) %>%
    dplyr::filter(y > min(plot.df$y)-0.5) %>%
    dplyr::filter(x > min(plot.df$x)-0.5)
  
  
  # Create a 3D scatter plot
  fig <- plotly::plot_ly() %>%
    plotly::add_trace(data=plot.df,
                      x=~x,
                      y=~y,
                      z=~z,
                      marker=list( color = plot.df$color.pt),
                      type="scatter3d", mode="markers") %>%
    plotly::add_trace(data=line1,x = ~x, y =~y,
                      z = ~z,
                      type = 'scatter3d', mode = 'lines', 
                      line = list(color = 'black', width = 4)) %>%
    plotly::add_trace(data=line2,x = ~x, y =~y,
                      z = ~z,
                      type = 'scatter3d', mode = 'lines', 
                      line = list(color = 'black', width = 4)) %>%
    plotly::add_trace(data=line3,x = ~x, y =~y,
                      z = ~z,
                      type = 'scatter3d', mode = 'lines', 
                      line = list(color = 'black', width = 4)) %>%
    # orange lines
    plotly::add_trace(data=orline1,x = ~x, y =~y,
                      z = ~z,
                      type = 'scatter3d', mode = 'lines', 
                      line = list(color = 'orange',
                                  width = 6,dash = 'dash', width = 6)) %>%
    plotly::add_trace(data=orline2,x = ~x, y =~y,
                      z = ~z,
                      type = 'scatter3d', mode = 'lines', 
                      line = list(color = 'orange',
                                  width = 6,dash = 'dash', width = 6)) %>%
    plotly::add_trace(data=orline3,x = ~x, y =~y,
                      z = ~z,
                      type = 'scatter3d', mode = 'lines', 
                      line = list(color = 'orange',
                                  width = 6,dash = 'dash', width = 6)) %>%
    # blue line
    plotly::add_trace(data=blueline,x = ~x, y =~y,
                      z = ~z,
                      type = 'scatter3d', mode = 'lines', 
                      line = list(color = 'blue',
                                  width = 6,dash = 'dash', width = 6)) %>%
    plotly::layout(showlegend = F, 
                   autosize = F, 
                   margin = m,
                   layout.title = list(text =  'c. Growth rate landscape of species i',
                                       family = "Courier New",
                                       size = 16,
                                       color = "black"), 
                   scene = list(xaxis = list(showgrid = F,
                                             showline = T,
                                             zerolinecolor = 'white',
                                             mirror = T, showticklabels = F),
                                yaxis = list(showgrid = F, 
                                             showline = T,
                                             zerolinecolor = 'white',
                                             mirror = T, showticklabels = F),
                                zaxis = list( showgrid = F, 
                                              showline = T,
                                              zerolinecolor = 'white',
                                              mirror = T, showticklabels = F)))
  
  fig <- fig %>% plotly::add_text(data=plot.df,
                                  x=~x,
                                  y=~y,
                                  z=~z,
                                  text = ~ point,
                                  textfont = list(color = '#000000',size = 16),
                                  textposition = "top",showlegend = FALSE)
  fig <- fig %>% plotly::add_mesh(data=plane.data,
                                   z=~z,
                                   x=~x,
                                   y=~y, 
                                   #type = 'scatter3d',
                                   #mode="surface",
                                   #color=~as.factor(col),
                                   #colors=c("lightblue","lightblue","lightblue"),
                                   #surfacecolor = list(color="lightblue"),
                                   opacity = 0.8)
  
  return(fig)
}

#graph_NFD_plan(pars)

graph_density_lv <- function(pars,time){
  
  Nt <- data.frame(Ni=0.01,Nj=pars$N_j_star,time=0)
  
  N0 <- c(0.01, pars$N_j_star)
  time_ode <- seq(0,time, 0.1)
  
  dN_dt <- function(t, N, pars){
    return(list(N*(pars$mu - pars$A%*%N)))
  }
  
  Nt <- ode(N0, time_ode, dN_dt, pars)
  colnames(Nt) <- c("time", "Ni", "Nj")
  
  # this was your previous version, I've changed to below
  # i didn't understand what you were attempting to do
  # also from my understanding you don't need to add these columns, but maybe I'm wrong...
  #Nt <-  Nt %>%
  #  as.data.frame() %>%
  #  mutate(mu_i = log10(unlist(abs(pars$mu_i[1]))+1),
  #         eta_i = log10(unlist(pars$eta_i)+1),
  #         r_i = log10(unlist(pars$r_i)+1))
  
  delta_time <- 0.1*time
  ggplot(Nt,aes(x=time))+
    geom_line(aes(y=Ni),color="red") +
    geom_line(aes(y=Nj),color="black") +
    theme_bw() +
    geom_segment(aes(y =Nt[1,"Ni"],yend = Nt[1, "Ni"] * exp(delta_time*pars$r_i))
                 ,x = 0, xend = delta_time, 
                 color="blue",
                 arrow = arrow(angle = 30, length = unit(0.25, "inches"),
                               ends = "last", type = "open")) +
    geom_text(aes(y = Nt[1, "Ni"] * exp(delta_time*pars$r_i)), label= "r_i",
              x=delta_time*0.8,size=6,color="darkblue") +
    geom_segment(aes(y =Nt[1,"Ni"],yend = Nt[1, "Ni"] * exp(delta_time*pars$mu_i))
                 ,x = 0, xend = delta_time, 
                 color="blue",
                 arrow = arrow(angle = 30, length = unit(0.25, "inches"),
                               ends = "last", type = "open")) +
    geom_text(aes(y = Nt[1, "Ni"] * exp(delta_time*pars$mu_i)), label="mu_i",
              x=delta_time*0.8,size=6,color="darkblue") +
    geom_segment(aes(y =Nt[1,"Ni"],yend = Nt[1, "Ni"] * exp(delta_time*pars$eta_i))
                 ,x = 0, xend = delta_time, 
                 color="blue",
                 arrow = arrow(angle = 30, length = unit(0.25, "inches"),
                               ends = "last", type = "open")) +
    geom_text(aes(y = Nt[1, "Ni"] * exp(delta_time*pars$eta_i)), label="eta_i",
              x=delta_time*0.8,size=6,color="darkblue") +
    geom_text(aes(y =Nj[1]), 
              label="Nj*",x=0, size=8,color="red") +
    geom_text(aes(y =Ni[1]),
              label="Ni",x=0,size=8,color="red") +
    scale_y_log10() +
    labs(#title="b. Abundance over time when species i invades species j \nand for species i in monoculturure",
      y="Abundance of species (log10)",x="time") +
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "bottom",
           strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=28),
           legend.text=element_text(size=20),
           legend.title=element_text(size=20),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_text(size=20, angle=66, hjust=1),
           axis.text.y= element_text(size=20),
           axis.title.x= element_text(size=24),
           axis.title.y= element_text(size=24),
           title=element_text(size=16))
}

#graph_density_lv(pars,20)

graph_NFD_2D <- function(pars,time){
  
  if(pars$ND > 2 | pars$ND < -1){
    ND.lim <- c(-1,2,pars$ND)
    ND.lim <- c(min(ND.lim),max(ND.lim))
  }else{ND.lim <- c(-1,2)}
  if(pars$FD > 2 | pars$FD < -1){
    FD.lim <- c(-1,2,pars$FD)
    FD.lim <- c(min(FD.lim),max(FD.lim))
  }else{FD.lim <- c(-1,2)}
  
  print(dev.size())
  ggplot() +
    geom_point(aes(x=pars$ND,y=pars$FD), color="red",size=4)+
    geom_abline(slope=1, intercept=c(0,0)) +
    geom_hline(yintercept=0) +
    geom_hline(yintercept=1) +
    geom_vline(xintercept=1) +
    geom_vline(xintercept=0) +
    scale_y_continuous(limits = FD.lim) +
    scale_x_continuous(limits = ND.lim) +
    annotate("text", x=c(-0.5,0.5,1.5,-1.0,-1,-1),
             y=c(-1.0,-1.0,-1.0,-0.5,0.5,1.5),
             size=5,angle=c(0,0,0,90,90,90),
             label=c("Stronger Inter\nthan Intra",
                     "Weaker Inter\nthan Intra",
                     "Inter and Intra\ndiffer in sign",
                     "Higher equil.\nabundance",
                     "Lower equil.\nabundance",
                     "Negative equil.\nabundance"),
             color="blue") +
    annotate("text", x=c(0.5,-0,1),
             y=c(0.5,1.5,1.5),
             size=5,angle=c(atan(dev.size()[2]/dev.size()[1])*180/pi,90,90),
             label=c(expression(paste(Ni^{'#'},'= 0')),
                     "No Frequency \ndependence",
                     "No species \ninteractions"),
             color="black") +
    labs(#title="d. Graph of niche and fitness difference MAP",
      y="Fitness difference",x="Niche differences") +
    theme_bw()+
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "bottom",
           strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=28),
           legend.text=element_text(size=20),
           legend.title=element_text(size=20),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_text(size=20, angle=66, hjust=1),
           axis.text.y= element_text(size=20),
           axis.title.x= element_text(size=24),
           axis.title.y= element_text(size=24),
           title=element_text(size=16))
}
#graph_NFD_2D(pars)

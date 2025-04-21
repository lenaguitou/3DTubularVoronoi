library(deldir)
library(ggplot2)
#library(ggvoronoi)
#library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)


# rad_coef = parameters$ratio_rad
# steps = parameters$n_steps
# Radius = parameters$apical_rad #5/(2*pi) #params$radiusA
# Radius2 <- rad_coef*Radius
# cyl_thickness <- Radius2-Radius
# cyl_length = parameters$cyl_length
# cyl_width = Radius*(2*pi)
# rec <- list()
# rad <- list()
# L = parameters$n_layers
# n = parameters$n_cells
# A0 <- ((Radius2+Radius)*pi*cyl_length)/n
# xmin <- 0
# ymin <- 0
# xmax <- cyl_width
# ymax <- cyl_length
# for(k in 1:L){
#   rad[[k]]<- Radius+(k-1)*(cyl_thickness/(L-1)) #the radius of the layer k
#   rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
# }
# lamad = parameters$lambda
# gamad = parameters$gamma



#FUNCTIONS INVOLVED IN THE ALGORITHM

#to see the functions in detail, see the N_cylinder_algorithm.R script
#The functions here are the same, excepting that we need more parameters, since
#we do the function call from another function (and so we don't define global
#parameters)
#
move_points <-function(pt,wid,len,rc,n,ind){ 
  theta <- runif(1,0,2*pi)  # random angle
  ptinx <- pt$x[ind] + rc*cos(theta)
  ptiny <- pt$y[ind] + rc*sin(theta)
  while((ptinx<0 || ptinx>wid)||(ptiny<0 || ptiny>len)){
    theta <- runif(1,0,2*pi)  # random angle
    ptinx <- pt$x[ind] + rc*cos(theta)
    ptiny <- pt$y[ind] + rc*sin(theta)}
  pt$x[c(ind,ind+n,ind+2*n)]<-c(ptinx,ptinx+wid,ptinx+2*wid)
  pt$y[c(ind,ind+n,ind+2*n)]<-ptiny
  return(pt)
}

move_points_2 <-function(pt,wid,len,rc,n){
  ind <- sample(1:n,1)
  ptinx <- pt$x[ind]+rnorm(1,mean=0,sd=rc)
  ptiny <- pt$y[ind]+rnorm(1,mean=0,sd=rc)
  while((ptinx<0 || ptinx>wid)||(ptiny<0 || ptiny>len)){
    ptinx <- pt$x[ind]+rnorm(1,mean=0,sd=rc)
    ptiny <- pt$y[ind]+rnorm(1,mean=0,sd=rc)
  }
  pt$x[c(ind,ind+n,ind+2*n)]<-c(ptinx,ptinx+wid,ptinx+2*wid)
  pt$y[c(ind,ind+n,ind+2*n)]<-ptiny
  return(pt)
}

tesellation_energy_N<-function(xt, yt, A0, rec, rad, gamad, lamad, n, Layer, s0){
  tesener<-sapply(1:Layer,function(i) {
    tesel<-deldir(xt*(rad[[i]]/rad[[1]]), yt, rw = rec[[i]])
    tilest<-tile.list(tesel)[(n+1):(2*n)]
    perims<-(tilePerim(tilest)$perimeters)/sqrt(A0)
    areas<-sapply(tilest,function(x){x$area})/A0
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/s0)
    sum((areas-1)^2+(gam/2)*(perims^2)+lamad*perims)
  })
  return(sum(tesener)/Layer)
}

tesellation_energy_N_parallel <- function(xt, yt, A0, rec, rad, gamad, lamad, n, Layer, s0) {
  tesener <- unlist(mclapply(1:Layer, function(i) {
    tesel <- deldir(xt * (rad[[i]] / rad[[1]]), yt, rw = rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    perims <- (tilePerim(tilest)$perimeters) / sqrt(A0)
    areas <- sapply(tilest, function(x) {x$area}) / A0
    gam <- gamad * exp((1 - (rad[[i]]/rad[[1]])) / s0)
    sum((areas-1)^2 + (gam/2) * (perims^2) + lamad * perims)
  }, mc.cores = detectCores() - 1))  # Use all but 1 core
  return(sum(tesener) / Layer)
}

elastic_contractile_adhesion_energy_parallel <- function(xt, yt, A0, rec, rad, gamad, lamad, n_cells, Layer, s0_ratio) {
energy_components <- do.call(rbind, mclapply(1:Layer, function(i) {
  tesel <- deldir(xt * (rad[[i]] / rad[[1]]), yt, rw = rec[[i]])
  tilest <- tile.list(tesel)[(n_cells + 1):(2 * n_cells)]
  perims <- (tilePerim(tilest)$perimeters) / sqrt(A0)
  areas <- sapply(tilest, function(x) { x$area }) / A0
  gam <- gamad * exp((1 - (rad[[i]] / rad[[1]])) / s0_ratio)
  # Compute energy components
  elastic_energy <- sum((areas - 1)^2) / Layer
  contractile_energy <- sum((gam / 2) * (perims^2)) / Layer
  adhesion_energy <- sum(lamad * perims) / Layer
  return(c(elastic_energy, contractile_energy, adhesion_energy))
}, mc.cores = detectCores() - 1))  # Use multiple cores
# Sum the components
return(list(
  ElasticEnergy = sum(energy_components[, 1]),
  ContractileEnergy = sum(energy_components[, 2]),
  AdhesionEnergy = sum(energy_components[, 3])
))}

elastic_contractile_adhesion_energy <- function(xt, yt, A0, rec, rad, gamad, lamad, n_cells, Layer, s0_ratio) {
  energy_components <- sapply(1:Layer, function(i){
    tesel <- deldir(xt * (rad[[i]] / rad[[1]]), yt, rw = rec[[i]])
    tilest <- tile.list(tesel)[(n_cells + 1):(2 * n_cells)]
    perims <- (tilePerim(tilest)$perimeters) / sqrt(A0)
    areas <- sapply(tilest, function(x) { x$area }) / A0
    gam <- gamad * exp((1 - (rad[[i]] / rad[[1]])) / s0_ratio)
    # Calculate energy components
    elastic_energy <- sum((areas - 1)^2)/Layer
    contractile_energy <- sum((gam / 2) * (perims^2))/Layer
    adhesion_energy <- sum(lamad * perims)/Layer
    # Return the energy components
    return(c(elastic_energy, contractile_energy, adhesion_energy))
  })
  
  # Return the sum of energy components for each layer
  return(list(
    ElasticEnergy = sum(energy_components[1, ]),
    ContractileEnergy = sum(energy_components[2, ]),
    AdhesionEnergy = sum(energy_components[3, ])
  ))
}

choice_metropolis<-function(delta,beta){
  p<-numeric()
  if(delta<=0){p<-1}else if(-delta*beta==0){p<-0}else{p<-exp(-delta*beta)}
  a=sample(c(0,1), size = 1, replace=TRUE, prob = c(1-p,p))
  return(a)
}

ggplotvor<-function(plotpoints,tit,wid){
  rectangle <- data.frame(x=c(xmin,xmin,xmax+2*wid,xmax+2*wid),y=c(ymin,ymax,ymax,ymin))
  pl <- ggplot(plotpoints,aes(x,y)) +
    geom_voronoi(aes(fill=as.factor(y)),size=.125, outline = rectangle,show.legend = FALSE) +
    geom_vline(xintercept = xmax,color = 'white',linetype='solid',size=1) +
    geom_vline(xintercept = xmax+wid,color = 'white',linetype='solid',size=1) +
    stat_voronoi(geom="path",outline = rectangle) +
    geom_point(size=2) +
    theme(
      panel.grid.major = element_blank() # Remove gridlines (major)
      ,panel.grid.minor = element_blank() # Remove gridlines (minor)
      ,panel.background = element_blank() # Remove grey background
      ,plot.title = element_text(hjust = 0, size = 20, colour = "#323232") # Title size and colour
      ,plot.caption = element_text(vjust = 0.3, size = 11, colour = "#323232") # Caption size and colour
      ,axis.ticks.y = element_blank() # Remove tick marks (Y-Axis)
      ,axis.text.y =  element_blank() # Remove scale marks (Y-Axis)
      ,axis.title.y = element_blank() # Remove axis label (Y-Axis) 
      ,axis.ticks.x = element_blank() # Remove tick marks (X-Axis)
      ,axis.text.x  = element_blank() # Remove axis scale (X-Axis)
      ,axis.title.x = element_blank() # Remove axis label (X-Axis) 
      ,legend.position="bottom"
    ) +
    labs(title = tit # Title text
         ,caption = "Author: Eloy Serrano        ")
  show(pl)
}

areasideplots<-function(px,py,rect,tit1,tit2,A0){
  tsl<-deldir(px,py,rw=rect)
  til<-tile.list(tsl)[(n_cells+1):(2*n_cells)]
  celledgearea<-data.frame(edges=numeric(),area=numeric())
  for (i in 1:length(til)) {
    celledgearea[i,c(1,2)]<-c(length(til[[i]]$x),til[[i]]$area/A0)
  }
  plotareaedges <- ggplot(celledgearea, aes(x = edges, y = area, colour = area))+
    geom_point()+xlab("Number of sides")+ylab("Relative area")+
    ggtitle(tit1)+
    stat_summary(aes(y = area,group=1), fun=mean, colour="#00BFC4", geom="line",group=1)
  show(plotareaedges)
  
  histedges <- ggplot(celledgearea,aes(edges))+geom_histogram(colour="#F8766D" ,fill="#7CAE00",bins=10,alpha=0.6)+
    xlab("Number of edges")+ylab("Frequency")+
    ggtitle(tit2)+
    xlim(2,11)
  show(histedges)
}

plotenergy <- function(en){
  ploten <- ggplot(en,aes(x=iteration,y=energy))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Tesselation energy")+
    scale_y_continuous(trans = "log")+
    ggtitle("Energy relaxation of the tesselation")
  show(ploten)
}


nu_sq <- function(points, rec, RadB= 2.5*5/(2*pi) , n_cells=100){
  Lay <- length(rec)
  teselap <- deldir(points$x, points$y, rw = rec[[1]])
  teselba <- deldir(RadB*points$x, points$y, rw = rec[[Lay]])
  tilap <- tile.list(teselap)[(n_cells+1):(2*n_cells)]
  tilba <- tile.list(teselba)[(n_cells+1):(2*n_cells)]
  
  cellsdf<-data.frame(edgesA=integer(),edgesB=integer())
  for (i in 1:length(tilap)) {
    cellsdf[i,c(1,2)]<-c(length(tilap[[i]]$x),length(tilba[[i]]$x))
  }
  num <- sum((cellsdf[,1]-cellsdf[,2])^2)/n_cells
  den <- 2*(sum(cellsdf[,1])/n_cells)*(sum(cellsdf[,2])/n_cells)
  return(num/den)
}


#Program starts

#We define a function to run the algorithm, 
#this way allows us to run simulations with different parameters at the same
#time, calling the function with different initial values
#The function works exactly as the script in N_cylinder_algorithm.R

metropolisad<-function(seed = 666, n_steps = 250, n_cells = 100, n_layers=5,
                       apical_rad = 5/(2*pi), ratio_rad = 2.5, cyl_length = 20,
                       gamma = 0.15, lambda = 0.04, beta = 100, s0_ratio=1,par=FALSE){
  
  #We define our variables
  
  RadiusB <- ratio_rad*apical_rad
  cyl_width_A <- 2*pi*apical_rad
  cyl_width_B <- 2*pi*RadiusB
  
  cyl_thickness <- RadiusB-apical_rad
  
  #We define the vertices of the plane
  xmin <- 0
  xmax <- cyl_width_A
  ymin <- 0
  ymax <- cyl_length
  
  r <- 1/10 #radius to make the moves dimensionless
  Am <- ((RadiusB+apical_rad)*pi*cyl_length)/n_cells
  
  rec <- vector(mode = "list", length = n_layers)
  rad <- numeric(n_layers)
  for(k in 1:n_layers){
    rad[[k]]<- apical_rad+(k-1)*(cyl_thickness/(n_layers-1)) #the radius of the layer k
    rec[[k]]<- c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
  }
  
  #Start, first iteration
  
  set.seed(seed)
  x1 <- runif(n_cells,xmin,xmax)
  y1 <- runif(n_cells,ymin,ymax)
  x <- c(x1,x1+cyl_width_A,x1+2*cyl_width_A)
  y <- c(y1,y1,y1)

  pointsinit <- data.frame(x=x,y=y)

  energyinit <- tesellation_energy_N(pointsinit$x, pointsinit$y, A0 = Am, rec = rec , rad =rad,
                                    gamad = gamma, lamad = lambda, n = n_cells, Layer = n_layers, s0 = s0_ratio)
  energiesinit <- elastic_contractile_adhesion_energy(pointsinit$x, pointsinit$y, A0 = Am, rec , rad,
                                     gamad = gamma, lamad = lambda, n = n_cells, Layer = n_layers, s0_ratio = s0_ratio)
  points <- pointsinit
  
  #We create the variables to store the results
  energhist <- data.frame(Iteration=numeric(n_steps), Total=numeric(n_steps),Elastic=numeric(n_steps),Contractile=numeric(n_steps),Adhesion=numeric(n_steps))

  energhist[1,c(1,2)] <- c(0,energyinit)
  energhist[1,c(3,4,5)] <- energiesinit
  
  histpts <- data.frame(x=numeric(3*n_cells*n_steps),
                        y=numeric(3*n_cells*n_steps),
                        Iteration = integer(3*n_cells*n_steps))
  
  histpts[1:(3*n_cells),c(1,2)] <- pointsinit
  histpts[1:(3*n_cells),3] <- 0
  
  energytesel <- energyinit


  #Start of the loop
  for (j in 1:n_steps) {
    print(paste0("Iteration ",j))
    if (par == TRUE){
      iter_start_time <- Sys.time()
      elapsed_time <- difftime(Sys.time(), start_time, units="secs")  
      log_msg <- sprintf("Simulation %d - Iteration %d - Time elapsed: %.2f sec  \n",seed,j, elapsed_time)
      cat(log_msg, file=log_file, append=TRUE)}
  
    indexes = sample(1:n_cells) #pick a random cell 
    for(k in indexes) {  
      points2 <- move_points(points,wid = cyl_width_A,len = cyl_length,rc = r,n = n_cells,ind = k)  
      energytesel2 <- tesellation_energy_N(points2$x, points2$y, A0 = Am,rec = rec , rad =rad,
                                         gamad = gamma, lamad = lambda, n = n_cells, Layer = n_layers, s0 = s0_ratio)
      c <- choice_metropolis(energytesel2-energytesel,beta) 
      cond <- c==1
      if(cond){
        points <- points2
        energytesel <- energytesel2
        energiestesel <- elastic_contractile_adhesion_energy(points$x, points$y, A0 = Am, rec = rec , rad = rad,
                                                             gamad = gamma, lamad = lambda, n_cells = n_cells, Layer = n_layers, s0_ratio = s0_ratio)
      }}
    gc()
    histpts[(j*3*n_cells+1):(j*3*n_cells+3*n_cells),c(1,2)] <- points
    histpts[(j*3*n_cells+1):(j*3*n_cells+3*n_cells),3]<-j
    energhist[j+1,c(1,2)] <- c(j,energytesel)
    energhist[j+1,c(3,4,5)] <- energiestesel
  }
  # nu2 <- nu_sq(points = points, rec = rec, n_cells = 100)
  return(list(points_evolution=histpts 
              ,energy_evolution=energhist
              # nu_2=nu2
              ))
}

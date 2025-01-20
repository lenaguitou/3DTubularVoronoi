library(deldir)
library(ggplot2)
#library(ggvoronoi)
library(stats)
library(dplyr)
library(nls.multstart)


funaux <-function(p,rec = c(0,15,0,20), n = 100,  cylen,apical_radius){
  
  A0 <- (2*pi*apical_radius*(cylen))/n
  
  #funaux returns a dataframe which first column is the edges of the cell and the
  # second column is the area of that cell
  
  tsl<-deldir(p$x,p$y,rw=rec)
  til<-tile.list(tsl)
  celledgearea<-data.frame(edges=integer(),area=double())
  for (i in 1:n) {
    celledgearea[i,c(1,2)]<-c(length(til[[i+n]]$x),til[[i+n]]$area/A0)
  }
  return (celledgearea)
}

funaux2sim<-function(ptsord, n = 100, ps = 100, cylen,apical_radius){
  
  A0 <- (2*pi*apical_radius*(cylen))/n
  #ps is the number of simulations done
  #First execute the ord function(for results)
  
  edgear<-data.frame(edges=integer(),area=double(),Frame=integer())
  for (i in 1:ps) {
    a<-length(edgear$edges)
    edgear[(a+1):(a+n),c(1,2)]<-funaux(ptsord[((i-1)*(3*n)+1):(i*3*n),c(1,2)],n, cylen,apical_radius)
    edgear[(a+1):(a+n),3]<-i
  }
  return(edgear)
}
funaux2simDOUBLE<-function(ptsord, n = 100, ps = 100, Ratio = 2.5,  cylen,apical_radius){
  
  A0 <- (2*pi*apical_radius*(cylen))/n
  #ps is the number of simulations done
  #First execute the ord function(for results)
  
  #This function returns a big dataframe, made up of one dataframe for each simulation
  # the df for one single simulation represents the edges and area of each of the cells.
  
  edgearA <- data.frame(edges=integer(), area=double(), Frame=integer())
  for (i in 1:ps) {
    a <- length(edgearA$edges)
    edgearA[(a+1):(a+n),c(1,2)] <- funaux(ptsord[((i-1)*(3*n)+1):(i*3*n),c(1,2)], rec=c(0,15,0,20), n, cylen,apical_radius)
    edgearA[(a+1):(a+n),3] <- i
  }
  
  ptsordB <- ptsord
  ptsordB$x <- ptsordB$x*Ratio
  
  edgearB <- data.frame(edges=integer(),area=double(),Frame=integer())
  for (i in 1:ps) {
    a <- length(edgearB$edges)
    edgearB[(a+1):(a+n),c(1,2)] <- funaux(ptsordB[((i-1)*(3*n)+1):(i*3*n),c(1,2)], rec=c(0,15*Ratio,0,20),n,cylen,apical_radius*Ratio)
    edgearB[(a+1):(a+n),3]<-i
  }
  return(list(edgearA,edgearB))
}

stationarylewisBasal<-function(edgear){
  #First execute ord function (for results) and funaux2 with the data, then this function makes the plots.

  plotareaedges <- ggplot(edgear, aes(x = edges, y = area, colour = area))+
    geom_point() + xlab("Number of sides") + ylab("Relative area")+
    ggtitle("Relative area of the cells by sides for the system in equilibrium. Basal surface")+
    stat_summary(aes(y = area, group = 1), fun = mean, colour = "#00BFC4", geom = "line", group = 1)
  show(plotareaedges)
  
  mined<-min(edgear$edges)
  maxed<-max(edgear$edges)
  quant<-maxed-mined
  
  datahist<-data.frame(edges=numeric(),frec=numeric())
  
  for (i in min(edgear$Frame):max(edgear$Frame)){
    for (j in mined:maxed) {
      dat<-dplyr::filter(edgear, edges == j & Frame == i)
      pos<-(i-1)*quant+j-mined+1
      datahist[pos,c(1,2)]<-c(j, length(dat$edges))
    }
  }
  
  meandat <- data.frame(edges=numeric(), meanfrec=numeric())
  varcoefdat <- data.frame(edges=numeric(),varcoef=numeric())
  
  for (j in mined:maxed) {
    dat2<-dplyr::filter(datahist, edges==j)
    pos<-j-mined+1
    meandat[pos,c(1,2)]<-c(j,mean(dat2$frec))
    varcoefdat[pos,c(1,2)]<-c(j,var(dat2$frec)/mean(dat2$frec))
  }
  
  print(datahist)
  print(meandat)
  print(varcoefdat)
  
  histedges<-ggplot(meandat,aes(edges,meanfrec/sum(meanfrec)))+
    geom_col(colour="#F8766D", fill="#7CAE00", alpha=0.6)+
    xlab("Number of edges")+ylab("Average frequency")+
    ggtitle("Average quantity of sides of the cells for system in equilibrium. Basal surface")+
    xlim(2,11)
  show(histedges)
}

stationarylewisApical<-function(edgear){
  #First execute ord function (for results) and funaux2 with the data, then this function makes the plots.
  
  plotareaedges <- ggplot(edgear, aes(x = edges, y = area, colour = area))+
    geom_point() + xlab("Number of sides") + ylab("Relative area")+
    ggtitle("Relative area of the cells by sides for the system in equilibrium. Apical surface")+
    stat_summary(aes(y = area, group = 1), fun = mean, colour = "#00BFC4", geom = "line", group = 1)
  show(plotareaedges)
  
  mined<-min(edgear$edges)
  maxed<-max(edgear$edges)
  quant<-maxed-mined
  
  datahist<-data.frame(edges=numeric(),frec=numeric())
  
  for (i in min(edgear$Frame):max(edgear$Frame)){
    for (j in mined:maxed) {
      dat<-dplyr::filter(edgear, edges == j & Frame == i)
      pos<-(i-1)*quant+j-mined+1
      datahist[pos,c(1,2)]<-c(j, length(dat$edges))
    }
  }
  
  meandat <- data.frame(edges=numeric(), meanfrec=numeric())
  varcoefdat <- data.frame(edges=numeric(),varcoef=numeric())
  
  for (j in mined:maxed) {
    dat2<-dplyr::filter(datahist, edges==j)
    pos<-j-mined+1
    meandat[pos,c(1,2)]<-c(j,mean(dat2$frec))
    varcoefdat[pos,c(1,2)]<-c(j,var(dat2$frec)/mean(dat2$frec))
  }
  
  print(datahist)
  print(meandat)
  print(varcoefdat)
  
  histedges<-ggplot(meandat,aes(edges,meanfrec/sum(meanfrec)))+
    geom_col(colour="#F8766D", fill="#7CAE00", alpha=0.6)+
    xlab("Number of edges")+ylab("Average frequency")+
    ggtitle("Average quantity of sides of the cells for system in equilibrium. Apical surface")+
    xlim(2,11)
  show(histedges)
}

ord <- function(results, iter = 150, n=100, sim=100){
  #With this function we extract the points for one iteration of every simulation
  #iter is the specific iteration of the algorithm that we extract to make the analysis
  iter = iter+1
  n_points = 3*n
  ptsord<-data.frame(x=numeric(n_points*sim),
                     y=numeric(n_points*sim), 
                     Frame=numeric(n_points*sim),
                     simulation=numeric(n_points*sim))
  for (i in 0:(sim-1)) {
    ptsord[(i*(n_points*iter)+1):((i+1)*(n_points*iter)),c(1,2,3)]<-results[[i+1]][[1]][c(1,2,3)]
    ptsord[(i*(n_points*iter)+1):((i+1)*(n_points*iter)),4]<-i+1
  }
  return(ptsord)
}

regnls<-function(energh){
  
  ##This functions makes a single regression for one iteration with the function
  #nls
  
  x<-unlist(lapply(energh$iteration,as.numeric))
  y<-unlist(lapply(energh$energy,as.numeric))
  m<-nls(y~I(a+b*(1-exp(-x/c))),start = list(a=2,b=0.3,c=0.1))
  #plot(x,y)
  #lines(x,predict(m),col="red",lwd=3)
  summary(m)
}

regnls2<-function(energh, n=100){
  
  #This functions makes a single regression for one iteration with the function
  #nls_multistart (which is much better for all cases)
  
  x<-unlist(lapply(energh$iteration,as.numeric))
  y<-unlist(lapply(energh$energy,as.numeric))
  m<-nls_multstart(y~I(a+b*(1-exp(-x/c))),
                   iter = 500,
                   start_lower = list(a=0,b=-5,c=5),
                   start_upper = list(a=5,b=5,c=30))
  #plot(x,y)
  #lines(x,predict(m),col="red",lwd=3)
  summary(m)
  return(summary(m)[["coefficients"]][c(1,2,3)])
}

adjsim<-function(results,nsim=100,it=150){
  
  #This function makes a regression of the energy of the results for each iteration,
  # and then it makes the average.
  
  coefest<-data.frame(a=double(nsim),b=double(nsim),c=double(nsim))
  for (i in 1:nsim) {
    coefest[i,c(1,2,3)]<-regnls2(results[[i,1]]$energy_evolution,n=nsim)
  }
  a<-mean(coefest$a)
  b<-mean(coefest$b)
  c<-mean(coefest$c)
  adj<-function(x) a+b*(1-exp(-x/c))
  adj_data<-data.frame(x=1:it,y=adj(1:it))
  ploten<-ggplot(data = adj_data, aes(x=x,y=y))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Average energy of the cells")+
    ggtitle("Average energy relaxation of the system")
  show(ploten)
  return(paste("Regression function for the energy: a+b*(1-exp(-x/c) with a=",a,"b=",b,"c=",c))
}

scutoids_analysis_oneiter<-function(pointsAx, pointsAy, pointsBx, pointsBy,ratio,ap_rad, n = 100,cylen){
  
  #We define the vertices of the plane
  xmin <- 0
  xmax <- 2*pi*ap_rad
  ymin <- 0
  ymax <- cylen

  rect1 <- c(xmin,xmin+3*xmax,ymin,ymax)
  rect2 <- c(xmin,xmin+3*(xmax*ratio),ymin,ymax)
  
  tslA<-deldir(pointsAx,pointsAy,rw=rect1)
  tilA<-tile.list(tslA)[(n+1):(2*n)]
  tslB<-deldir(pointsBx,pointsBy,rw=rect2)
  tilB<-tile.list(tslB)[(n+1):(2*n)]
  
  cellsdf<-data.frame(edgesA=integer(),edgesB=integer())
  
  for (i in 1:length(tilA)) {
    cellsdf[i,c(1,2)]<-c(length(tilA[[i]]$x),length(tilB[[i]]$x))
  }
  countdf<- cellsdf %>%
    group_by(edgesA,edgesB) %>%
    summarize(count=n())
  
  total_count <- sum(countdf$count)
  
  # Convert counts to percentages
  countdf$percent <- (countdf$count / total_count) * 100

  # Filter to only keep percentages greater  
  #countdf <- dplyr::filter(countdf, percent >= 1)
  
  scutoidsplot<-ggplot(countdf, aes(x = edgesA, y = edgesB, label=count))+
    geom_count(shape = "square", aes(color= percent,size=10))+
    xlab("Edges on apical surface")+ylab("Edges on basal surface")+
    ggtitle("Polygon class of apical and basal surfaces")+
    guides(colour = "colorbar", size = "none")+
    labs(color = "Percentages")+
    scale_size_area(max_size = 70)+
     #geom_label()+
    #scale_color_gradient(low = "yellow", high = "blue") +
    scale_color_viridis(option = "D",direction = -1, limits = c(0,30)) +
    xlim(4,8) +
    ylim(4,8)
  show(scutoidsplot)
  
  }

scutoids_percr <-function(histpts, n = 100, it=150,ratio,ap_rad,cylen,plotShow=TRUE){
  
  #We define the vertices of the plane
  xmin <- 0
  xmax <- 2*pi*ap_rad
  ymin <- 0
  ymax <- cylen
  
  rect1 <- c(xmin,xmin+3*xmax,ymin,ymax)
  rect2 <- c(xmin,xmin+3*(xmax*ratio),ymin,ymax)
                        
  #Computes the evolution of the percentage of escutoids in every iteration of the algorithm

  perc<-data.frame(it=integer(it), perc_sc=double(it))
  for (i in 1:it) {
    pointsAx<-filter(histpts, Frame == i)$x
    pointsAy<-filter(histpts, Frame == i)$y
    pointsBx <- pointsAx*2.5
    pointsBy <- pointsAy*2.5
    tslA<-deldir(pointsAx,pointsAy,rw=rect1)
    tilA<-tile.list(tslA)[(n+1):(2*n)]
    tslB<-deldir(pointsBx,pointsBy,rw=rect2)
    tilB<-tile.list(tslB)[(n+1):(2*n)]
    cellsdf<-data.frame(edgesA=integer(),edgesB=integer())
    for (j in 1:length(tilA)) {
      cellsdf[j,c(1,2)]<-c(length(tilA[[j]]$x),length(tilB[[j]]$x))
    }
    #percen<-length(filter(cellsdf, edgesA==edgesB)[[1]])
    #perc[i,c(1,2)]<-c(i,percen)
    
    # Count the number of matching edges (where edgesA == edgesB)
    matching_tiles <- length(filter(cellsdf, edgesA == edgesB)[[1]])
    
    # Calculate the total number of tiles in this iteration
    total_tiles <- nrow(cellsdf)
    
    # Calculate the percentage of matching tiles
    perc_sc <- (matching_tiles / total_tiles) * 100
    
    # Store the percentage for this iteration
    perc[i, c(1, 2)] <- c(i, perc_sc)
  }
  
  if(plotShow==TRUE){
    ploten<-ggplot(perc,aes(x=it,y=perc_sc))+
      geom_line(colour="#F8766D")+
      xlab("Iteration of the algorithm")+
      ylab("Percentage of escutoids")+
      ggtitle("Evolution of the percentage of escutoids of the system")
    show(ploten)
  }
  return(perc)
}

  

scutoids_percr_simulations<-function(results, n=100, sim=100, it=150,ap_rad,cylen,ratio){
  xmin <- 0
  xmax <- 2*pi*ap_rad
  ymin <- 0
  ymax <- cylen
  
  rect1 <- c(xmin,xmin+3*xmax,ymin,ymax)
  rect2 <- c(xmin,xmin+3*(xmax*ratio),ymin,ymax)
  
  
  perc<-data.frame(it=1:it, percen=rep(0,it))
  for (j in 1:sim) {
    percen<-scutoids_percr(results[[j,1]]$points_evolution, n, it,ratio,ap_rad,cylen,plotShow = FALSE)$perc_sc
    perc$percen<-perc$percen + percen
  }
  perc$percen<-perc$percen/sim
  ploten<-ggplot(perc,aes(x=it,y=percen))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Percentage of escutoids")+
    ggtitle("Evolution of the average percentage of escutoids")
  show(ploten)
  
  return(perc)
}

scutoids_percr_simulations_par <- function(results, n=100, sim=100, it=150, ap_rad, cylen, ratio) {
  
  # Set up parallel backend (adjust the number of cores as necessary)
  num_cores <- detectCores() - 1  # Leave one core free for other processes
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  xmin <- 0
  xmax <- 2 * pi * ap_rad
  ymin <- 0
  ymax <- cylen
  
  rect1 <- c(xmin, xmin + 3 * xmax, ymin, ymax)
  rect2 <- c(xmin, xmin + 3 * (xmax * ratio), ymin, ymax)
  
  perc <- data.frame(it = 1:it, percen = rep(0, it))
  
  # Use foreach for parallel execution
  foreach(j = 1:50, .combine = 'cbind', .packages = c('deldir', 'dplyr'), .export = c('scutoids_percr')) %dopar% {
    
    # Fetch points for the current simulation
    percen <- scutoids_percr(results[[j, 1]]$points_evolution, n, it, ratio, ap_rad, cylen, plotShow = FALSE)$perc_sc
    
    return(percen)  # Return the perc_sc values for each simulation
  } -> perc_results
  # Stop the parallel cluster
  stopCluster(cl)
  
  # Summing the perc_results and averaging over simulations
  perc$percen <- rowSums(perc_results) / 50
  
  # Plotting the results
  ploten <- ggplot(perc, aes(x = it, y = percen)) +
    geom_line(colour = "#F8766D") +
    xlab("Iteration of the algorithm") +
    ylab("Percentage of escutoids") +
    ggtitle("Evolution of the average percentage of escutoids")
  show(ploten)
  
  
  return(perc_results)
}


scutoids_prep <- function(pointsAx,pointsAy,pointsBx,pointsBy,ratio,ap_rad, n = 100,cylen){

  
  #This function returns a dataframe which counts the number of edges in the 
  # apical and basal surface of each cell, and counts the cells that has the same
  # number of edges in both surfaces.
  #We define the vertices of the plane
  xmin <- 0
  xmax <- 2*pi*ap_rad
  ymin <- 0
  ymax <- cylen
  
  rect1 <- c(xmin,xmin+3*xmax,ymin,ymax)
  rect2 <- c(xmin,xmin+3*(xmax*ratio),ymin,ymax)
  
  
  tslA <- deldir(pointsAx,pointsAy,rw=rect1)
  tilA <- tile.list(tslA)[(n+1):(2*n)]
  tslB <- deldir(pointsBx,pointsBy,rw=rect2)
  tilB <- tile.list(tslB)[(n+1):(2*n)]
  cellsdf <- data.frame(edgesA=integer(),edgesB=integer())
  for (i in 1:n) {
    cellsdf[i,c(1,2)] <- c(length(tilA[[i]]$x),length(tilB[[i]]$x))
  }
  countdf <- cellsdf %>%
    count(edgesA,edgesB)
  return(countdf)
}

scutoids_analysis_stationary <- function(histpts, rect1, rect2, n = 100){
  lon <- 50 #how many iterations we want to have
  histdf_count <- data.frame(edgesA=integer(lon*n),edgesB=integer(lon*n),count=integer(lon*n))
  for (i in 1:lon) {
    histdf_count[((i-1)*n+1):(i*n)]<-
      scutoids_prep(dplyr::filter(histpts, Frame==250+i-1)$x,
                    dplyr::filter(histpts, Frame==250+i-1)$y,
                    (Radius2/Radius)*dplyr::filter(histpts, Frame==250+i-1)$x,
                    dplyr::filter(histpts, Frame==250+i-1)$y,
                    rect1,rect2)
  }
  histdf_avgcount <- histdf_count %>%
    group_by(edgesA,edgesB) %>%
    summarize(avg_count=mean(count))
  
  scutoidsplot<-ggplot(histdf_avgcount, aes(x = edgesA, y = edgesB, label=avg-count))+
    geom_count(shape = "square", aes(color= avg_count))+
    xlab("Average of edges on apical surface")+ylab("Average of edges on basal surface")+
    ggtitle("Average polygon class of apical and basal surfaces")+
    guides(colour = "colorbar", size = "none")+
    scale_size_area(max_size = 30)+
    geom_label()+
    scale_fill_gradient(low = "light blue", high = "deepskyblue")+
    xlim(3.3,7.5)+
    ylim(3.5,7.5)
  show(scutoidsplot)
  return(histdf_avgcount)
}

scutoids_analysis_simulations <- function(results, Ratio = 2.5, n = 100, sim = 100, it = 150, cylen,ap_rad){
  #We define the vertices of the plane
  xmin <- 0
  xmax <- 2*pi*ap_rad
  ymin <- 0
  ymax <- cylen
  
  rect1 <- c(xmin,xmin+3*xmax,ymin,ymax)
  rect2 <- c(xmin,xmin+3*(xmax*Ratio),ymin,ymax)
  
  #sim is how many simulations we have
  #it is the iteration that we have
  
  
  histdf_count <- data.frame(edgesA=integer(),
                             edgesB=integer(),
                             count=integer())
  
  for (i in 1:sim) {
    #each simulation have a different distribution of sides, so we have to adapt
    
    a <- length(histdf_count[[1]])
    ptsx <- dplyr::filter(results[[i]]$points_evolution, Frame==it)$x
    ptsy <- dplyr::filter(results[[i]]$points_evolution, Frame==it)$y
    df<-scutoids_prep(ptsx, ptsy,
                      Ratio*ptsx, ptsy,Ratio,ap_rad,n,cylen)
    len<-length(df[[1]])
    histdf_count[(a+1):(a+len),c(1,2,3)]<-df
  }
  
  #Now, to compute the average, we have to take into account the 1000 simulations
  
  histdf_avcount <- data.frame(edgesA = double(),
                               edgesB = double(),
                               avg_count = double())

  for (edA in min(histdf_count$edgesA):max(histdf_count$edgesA)) {
    for (edB in min(histdf_count$edgesB):max(histdf_count$edgesB)) {
      a <- length(histdf_avcount$edgesA)
      histdf_avcount[a+1,c(1,2,3)] <-
        c(edA, edB, (sum(dplyr::filter(histdf_count, edgesA == edA & edgesB == edB)$count)*100/(sim*100)))
    }
  }
  

  #histdf_avcount <- dplyr::filter(histdf_avcount, avg_count > 1)


  scutoidsplot<-ggplot(histdf_avcount, aes(x = edgesA, y = edgesB, label =avg_count))+
    geom_count(shape = "square", aes(color= avg_count,size=10))+
    xlab("Average of edges on apical surface")+ylab("Average of edges on basal surface")+
    ggtitle("Average polygon class of apical and basal surfaces")+
    guides(colour = "colorbar", size = "none")+
    labs(color = "Percentages")+
    scale_size_area(max_size = 70)+
    #geom_label()+
    #scale_fill_gradient(low = "light blue", high = "deepskyblue")+
    scale_color_viridis(option = "D",direction = -1, limits = c(0,35)) +
    xlim(3,9)+
    ylim(3,9)
  # xlim(min(histdf_avcount$edgesA)-0.5, max(histdf_avcount$edgesA)-0.5)+
  # ylim(min(histdf_avcount$edgesB)-0.5, max(histdf_avcount$edgesB)-0.5)
  show(scutoidsplot)
}









# ADDED 

energy_analysis <- function(points){
  
  
  tesener <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                         Elastic_energy = double(Lay),
                         Tension_energy = double(Lay),
                         Contractile_energy = double(Lay),
                         Bending_energy = double(Lay),
                         Total_energy=double(Lay))
  
  for (i in 1:Lay) {
    tesel <- deldir(points[[i]]$x,points[[i]]$y,rw=rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    
    perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
    areas <- sapply(tilest,function(x){x$area/A0})
    
    elener <- (sum((areas-1)^2)/n)
    tenener <- sum(lamad*perims)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contener <- sum((gam/2)*(perims^2))/n
    tesener[i,c(2,3,4,6)] <- c(elener/Lay,
                               tenener/lay,
                               contener/Lay,
                               sum((areas-1)^2+(gam/2)*(perims^2)+
                                     lamad*perims)/(n*Lay))
    
  }
  
  tesener$Bending_energy[c(1,Lay)]<-0
  
  tesener$Bending_energy[2:(Lay-1)] <- sapply(2:(Lay-1), function(i){
    
    angles <- sapply(1:n, function(j){
      ptcentral <- c(points[[i]]$y[j],
                     rad[[i]]*cos((1/rad[[i]])*points[[i]]$x[j]),
                     rad[[i]]*sin((1/rad[[i]])*points[[i]]$x[j]))
      
      ptinf <- c(points[[i-1]]$y[j],
                 rad[[i-1]]*cos((1/rad[[i-1]])*points[[i-1]]$x[j]),
                 rad[[i-1]]*sin((1/rad[[i-1]])*points[[i-1]]$x[j]))
      
      ptsup <- c(points[[i+1]]$y[j],
                 rad[[i+1]]*cos((1/rad[[i+1]])*points[[i+1]]$x[j]),
                 rad[[i+1]]*sin((1/rad[[i+1]])*points[[i+1]]$x[j]))
      
      vec1 <- ptsup - ptcentral
      vec2 <- ptinf - ptcentral
      
      #We use pmin and pmax to avoid errors in the arc-cosine computation
      v <- pmin(pmax(((vec1%*%vec2)[1,1])/
                       (norm(vec1,type = "2")*norm(vec2,type = "2")),-1.0),1.0)
      ang <- acos(v)        
      return(ang)
    })
    
    return(sum(alpha*((angles-pi)^2))/(Lay*n))
  })
  
  tesener$Total_energy <- tesener$Total_energy+tesener$Bending_energy
  return(tesener)
}

energy_analisis_1sim <- function(histpts, it = 150, lay = 10, n =100){
  histener<-data.frame(it = 1:(it+1), elen = double(it+1), tenen = double(it+1),
                       conten = double(it+1), toten = double(it+1))
  #print(head(histpts))
  for (i in 1:(it+1)) {
    pts <- filter(histpts, Frame==i)
    histener[i,c(2,3,4,5)] <- energy_iteration(pts$x,pts$y, n = n, Lay = lay)
  }
  return(histener)
}


energy_analisis_averages <- function(results, it = 150, lay = 10, n = 100) {
  # Number of simulations
  N_SIM <- nrow(results)
  
  # Initialize a data frame to store averaged energy values
  average_energy <- data.frame(it = 1:(it + 1), elen = double(it + 1), 
                               tenen = double(it + 1), conten = double(it + 1), 
                               toten = double(it + 1))
  
  # Loop through all simulations
  for (i in 1:N_SIM) {
    # Get the historical points from `points_evolution`
    histpts <- results[[i, 1]]$energy_evolution
    
    # Analyze energy for the current simulation
    histener <- energy_analisis_1sim(histpts, it = it, lay = lay, n = n)
    
    # If it's the first simulation, initialize average_energy
    if (i == 1) {
      average_energy <- histener
    } else {
      # Update averages by summing the energies
      average_energy[, 2:5] <- average_energy[, 2:5] + histener[, 2:5]
    }
  }
  
  # Compute the final averages
  average_energy[, 2:5] <- average_energy[, 2:5] / N_SIM
  
  
  p<-ggplot(average_energy, aes(x = it))+
    geom_line(aes(y = toten, colour = "Total energy"))+
    geom_line(aes(y= tenen, colour = "Adhesion"))+
    geom_line(aes(y = elen, colour = "Elastic"))+
    geom_line(aes(y = conten, colour = "Contractility"))+
    # geom_line(aes(x = Layer, y = Bending_energy, colour = "Bending energy"))+
    labs(title = "Decomposition of system energies",
         x = "Iteration", y = "Average energy per cell",
         color = "Energy type") +
    scale_colour_manual("",
                        breaks = c("Total energy",
                                   "Adhesion",
                                   "Elastic",
                                   "Contractility"),
                        # "Bending energy"),
                        values = c("red",
                                   "blue",
                                   "darkgreen",
                                   "purple"))
  # "orange",
  
  
  show(p)
  
  return(average_energy)
}



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
  
  # Calculate the areas of all the Voronoi cells
  all_areas <- sapply(til[(n+1):(2*n)], function(cell) cell$area) # Extract areas for the n cells of interest
  average_area <- mean(all_areas)  # Compute the average area
  
  celledgearea<-data.frame(edges=integer(),area=double())
  for (i in 1:n) {
    #celledgearea[i,c(1,2)]<-c(length(til[[i+n]]$x),til[[i+n]]$area/A0)
    celledgearea[i,c(1,2)]<-c(length(til[[i+n]]$x),til[[i+n]]$area/average_area)
  }
  return (celledgearea)
}

funaux2sim<-function(ptsord, n = 100, ps = 100, cylen,apical_radius,it){
  
  A0 <- (2*pi*apical_radius*(cylen))/n
  #ps is the number of simulations done
  #First execute the ord function(for results)
  
  edgear<-data.frame(edges=integer(),area=double(),Iteration=integer())
  for(j in 1:ps){
    for (i in 1:(it+1)) {
      a<-length(edgear$edges)
      edgear[(a+1):(a+n),c(1,2)] <- funaux(ptsord[((i-1)*(3*n) + 1 + n*(j-1)):(i*3*n + n*(j-1)),c(1,2)],n = n, cylen = cylen,apical_radius = apical_radius)
      edgear[(a+1):(a+n),3]<-i-1
  }
  }
  return(edgear)
}

funaux2simDOUBLE<-function(ptsord, n = 100, ps = 100, Ratio = 2.5,cylen,apical_radius,it){
  
  A0 <- (2*pi*apical_radius*(cylen))/n
  #ps is the number of simulations done
  #First execute the ord function(for results)
  
  #This function returns a big dataframe, made up of one dataframe for each simulation
  # the df for one single simulation represents the edges and area of each of the cells.
  
  edgearA <- data.frame(edges=integer(), area=double(), Iteration=integer())
  for(j in 1:ps){
    for (i in 1:(it+1)) {
    a <- length(edgearA$edges)
    edgearA[(a+1):(a+n),c(1,2)] <- funaux(ptsord[((i-1)*(3*n) + 1 + n*(j-1)):(i*3*n + n*(j-1)),c(1,2)], rec=c(0,0+3*(2*pi*apical_radius),0,cylen), n=n, cylen = cylen,apical_radius=apical_radius)
    edgearA[(a+1):(a+n),3] <- i -1
  }}
  
  ptsordB <- ptsord
  ptsordB$x <- ptsordB$x*Ratio
  
  edgearB <- data.frame(edges=integer(),area=double(),Iteration=integer())
  for(j in 1:ps){
    for (i in 1:(it+1)) {
      a <- length(edgearB$edges)
      edgearB[(a+1):(a+n),c(1,2)] <- funaux(ptsordB[((i-1)*(3*n) + 1 + n*(j-1)):(i*3*n + n*(j-1)),c(1,2)], rec=c(0,0+3*(2*pi*apical_radius)*Ratio,0,cylen),n=n,cylen=cylen,apical_radius=apical_radius*Ratio)
      edgearB[(a+1):(a+n),3]<-i-1
    }
  }
  return(list(edgearA,edgearB))
}

stationarylewisBasal<-function(edgear,it,nsim=1){
  #First execute ord function (for results) and funaux2 with the data, then this function makes the plots.

  # plotareaedges <- ggplot(edgear, aes(x = edges, y = area, colour = area))+
  #   geom_point() + xlab("Number of sides") + ylab("A/<A>")+xlim(3,9)+ylim(0,2)+
  #   ggtitle("Relative area of the cells by sides for the system in equilibrium. Basal surface")+
  #   stat_summary(aes(y = area, group = 1), fun = mean, colour = "#00BFC4", geom = "line", group = 1)
  
  plotareaedges <- ggplot(edgear, aes(x = edges, y = area, colour = area)) +
    geom_point(alpha = 0.6, size = 3) +  # Add transparency and adjust point size
    stat_summary(aes(y = area, group = 1), fun = mean, colour = "#0072B2", 
                 geom = "line", group = 1, linewidth = 1.2) +  # Line for mean values
    scale_colour_gradient(low = "#56B1F7", high = "#132B43", name = "Area") +  # Gradient color for points
    labs(
      title = paste("Cells area at iteration",it),
      subtitle = "Basal Surface",
      x = "Number of Sides",
      y = "A/<A>"
    ) +
    xlim(3, 9) + 
    ylim(0, 2) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "right",
      legend.key.size = unit(0.8, "lines")
    )
  
  show(plotareaedges)
  
  mined<-min(edgear$edges)
  maxed<-max(edgear$edges)
  quant<-maxed-mined
  
  datahist<-data.frame(edges=numeric(),frec=numeric())
  
  for (i in min(edgear$Iteration):max(edgear$Iteration)){
    for (j in mined:maxed) {
      dat<-dplyr::filter(edgear, edges == j & Iteration == i)
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
  
  #print(datahist)
  #print(meandat)
  #print(varcoefdat)
  
  # histedges<-ggplot(meandat,aes(edges,meanfrec/sum(meanfrec)))+
  #   geom_col(colour="#F8766D", fill="#7CAE00", alpha=0.6)+
  #   xlab("Number of edges")+ylab("Average frequency")+
  #   ggtitle("Average quantity of sides of the cells for system in equilibrium. Basal surface")+
  #   xlim(2,11)
  
  histedges <- ggplot(meandat, aes(x = edges, y = meanfrec / sum(meanfrec))) +
    geom_col(
      colour = "#2C3E50",  # Dark border for bars
      fill = "#1ABC9C",    # Teal fill color
      alpha = 0.8         # Slight transparency for bars
    ) +
    labs(
      title = paste("Number of neighbors at iteration",it),
      subtitle = "Basal Surface",
      x = "Number of Edges",
      y = "Normalized Frequency"
    ) +
    scale_x_continuous(limits = c(3, 9), breaks = seq(3, 9, by = 1)) +  # Fine-tuned x-axis
    scale_y_continuous(limits = c(0,0.75),labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05))) +  # Y-axis as percentage
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),  # Subtle grid lines
      panel.grid.minor = element_blank(),                # No minor grid lines
      legend.position = "none"                           # Remove legend
    )
  
  show(histedges)
  
  if(nsim!=1){
    output_file1 <- paste("results/figures/Number of neighbors basal iteration",it, "_Average",nsim,"sim.")
    output_file2 <- paste("results/figures/Cells area basal iteration",it, "_Average",nsim,"sim.")
  }else{
    output_file1 <- paste("results/figures/Number of neighbors basal iteration",it)
    output_file2 <- paste("results/figures/Cells area basal iteration",it)
  }
  
  # Save the plot
  ggsave(
    filename = output_file1,    # File name
    plot = histedges,                  # ggplot object
    device = "svg",            # File format
    width = 8,                 # Width in inches
    height = 6,                # Height in inches
    dpi = 300                  # Resolution (optional for SVG)
  )
  
  ggsave(
    filename = output_file2,  # File name
    plot = plotareaedges,      # ggplot object
    device = "svg",            # File format
    width = 8,                 # Width in inches
    height = 6,                # Height in inches
    dpi = 300                  # Resolution (optional for SVG)
  )
  
  # Confirmation
  print(paste0("Plot saved:", output_file1))
  print(paste0("Plot saved:", output_file2))
  
}

stationarylewisApical<-function(edgear,it,nsim=1){
  #First execute ord function (for results) and funaux2 with the data, then this function makes the plots.
  
  # plotareaedges <- ggplot(edgear, aes(x = edges, y = area, colour = area))+
  #   geom_point() + xlab("Number of sides") + ylab("A/<A>")+ xlim(3,9)+ylim(0,2)+
  #   ggtitle("Relative area of the cells by sides for the system in equilibrium. Apical surface")+
  #   stat_summary(aes(y = area, group = 1), fun = mean, colour = "#00BFC4", geom = "line", group = 1)
  plotareaedges <- ggplot(edgear, aes(x = edges, y = area, colour = area)) +
    geom_point(alpha = 0.6, size = 3) +  # Add transparency and adjust point size
    stat_summary(aes(y = area, group = 1), fun = mean, colour = "#0072B2", 
                 geom = "line", group = 1, linewidth = 1.2) +  # Line for mean values
    scale_colour_gradient(low = "#56B1F7", high = "#132B43", name = "Area") +  # Gradient color for points
    labs(
      title = paste("Cells area at iteration",it),
      subtitle = "Apical Surface",
      x = "Number of Sides",
      y = "A/<A>"
    ) +
    xlim(3, 9) + 
    ylim(0, 2) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "right",
      legend.key.size = unit(0.8, "lines")
    )
  
  show(plotareaedges)
  
  mined<-min(edgear$edges)
  maxed<-max(edgear$edges)
  quant<-maxed-mined
  
  datahist<-data.frame(edges=numeric(),frec=numeric())
  
  for (i in min(edgear$Iteration):max(edgear$Iteration)){
    print(i)
    for (j in mined:maxed) {
      dat<-dplyr::filter(edgear, edges == j & Iteration == i)
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
  
  #print(datahist)
  #print(meandat)
  #print(varcoefdat)
  
  # histedges<-ggplot(meandat,aes(edges,meanfrec/sum(meanfrec)))+
  #   geom_col(colour="#F8766D", fill="#7CAE00", alpha=0.6)+
  #   xlab("Number of edges")+ylab("Average frequency")+
  #   ggtitle("Average quantity of sides of the cells for system in equilibrium. Apical surface")+
  #   xlim(2,11)
  
  histedges <- ggplot(meandat, aes(x = edges, y = meanfrec / sum(meanfrec))) +
    geom_col(
      colour = "#2C3E50",  # Dark border for bars
      fill = "#1ABC9C",    # Teal fill color
      alpha = 0.8         # Slight transparency for bars
    ) +
    labs(
      title = paste("Number of neighbors at iteration",it),
      subtitle = "Apical Surface",
      x = "Number of Edges",
      y = "Normalized Frequency"
    ) + 
    scale_x_continuous(limits = c(3, 9), breaks = seq(3, 9, by = 1)) +  # Fine-tuned x-axis
    scale_y_continuous(limits = c(0,0.75),labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05))) +  # Y-axis as percentage
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),  # Subtle grid lines
      panel.grid.minor = element_blank(),                # No minor grid lines
      legend.position = "none"                           # Remove legend
    )
  
  
  show(histedges)
  
  if(nsim!=1){
    output_file1 <- paste("results/figures/Number of neighbors apical iteration",it, "_Average",nsim,"sim.")
    output_file2 <- paste("results/figures/Cells area apical iteration",it, "_Average",nsim,"sim.")
  }else{
    output_file1 <- paste("results/figures/Number of neighbors apical iteration",it)
    output_file2 <- paste("results/figures/Cells area apical iteration",it)
  }
 
  
  # Save the plot
  ggsave(
    filename = output_file1,    # File name
    plot = histedges,                  # ggplot object
    device = "svg",            # File format
    width = 8,                 # Width in inches
    height = 6,                # Height in inches
    dpi = 300                  # Resolution (optional for SVG)
  )
  
  ggsave(
    filename = output_file2,  # File name
    plot = plotareaedges,      # ggplot object
    device = "svg",            # File format
    width = 8,                 # Width in inches
    height = 6,                # Height in inches
    dpi = 300                  # Resolution (optional for SVG)
  )
  
  # Confirmation
  print(paste0("Plot saved:", output_file1))
  print(paste0("Plot saved:", output_file2))
  
}

ord <- function(results, iter = 150, n=100, sim=100){
  #With this function we extract the points for one iteration of every simulation
  #iter is the specific iteration of the algorithm that we extract to make the analysis
  iter = iter+1
  n_points = 3*n
  ptsord<-data.frame(x=numeric(n_points*sim),
                     y=numeric(n_points*sim), 
                     Iteration=numeric(n_points*sim),
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

regnls2<-function(energh,it){
  
  #This functions makes a single regression for one iteration with the function
  #nls_multistart (which is much better for all cases)
  
  x<-unlist(lapply(energh$Iteration,as.numeric))
  y<-unlist(lapply(energh$energy,as.numeric))
  m<-nls_multstart(y~I(a+b*(1-exp(-x/c))),
                   iter = it,
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
    coefest[i,c(1,2,3)]<-regnls2(results[[i,1]]$energy_evolution$Total,it=it)
  }
  a<-mean(coefest$a)
  b<-mean(coefest$b)
  c<-mean(coefest$c)
  adj<-function(x) a+b*(1-exp(-x/c))
  adj_data<-data.frame(x=0:it,y=adj(0:it))
  ploten<-ggplot(data = adj_data, aes(x=x,y=y))+
    geom_line(colour="#F8766D")+
    xlab("Iteration of the algorithm")+
    ylab("Average energy of the cells")+
    ggtitle("Average energy relaxation of the system")+
  ylim(0, max(adj_data$y) + 0.1 * max(adj_data$y))
  show(ploten)
  return(print(paste0("Regression function for the energy: a+b*(1-exp(-x/c) with a=",a,"b=",b,"c=",c)))
}

scutoids_analysis_oneiter<-function(pointsAx, pointsAy, pointsBx, pointsBy,ratio,ap_rad, n = 100,cylen,it){
  
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
    summarize(count=n(), .groups="drop")
  
  total_count <- sum(countdf$count)
  
  # Convert counts to percentages
  countdf$percent <- (countdf$count / total_count) * 100

  # Filter to only keep percentages greater  
  #countdf <- dplyr::filter(countdf, percent >= 1)
  
  scutoidsplot<-ggplot(countdf, aes(x = edgesA, y = edgesB, label=count))+
 #   geom_count(shape = "square", aes(color= percent,size=10))+
 # xlab("Edges on apical surface")+ylab("Edges on basal surface")+
 # ggtitle("Polygon class of apical and basal surfaces")+
 # guides(colour = "colorbar", size = "none")+
 # labs(color = "Percentages")+
 # scale_size_area(max_size = 70)+
 #  #geom_label()+
 # #scale_color_gradient(low = "yellow", high = "blue") +
 # scale_color_viridis(option = "D",direction = -1, limits = c(0,30)) +
 # xlim(4,8) +
 # ylim(4,8)
 #  
  geom_point(aes(color = percent), shape = 15, size=25) + # Square shape (15), size mapped to avg_count
    scale_color_viridis_c(
      option = "D",
      direction = -1,
      limits = c(0, 65),
      name = "Percentages"
    ) +
    scale_size_continuous(range = c(3, 9), guide = "none") + # Controlled size range
    scale_x_continuous(
      breaks = c(3, 4, 5, 6, 7, 8, 9),  # Custom x-axis ticks
      limits = c(3, 9.5),                 # Explicit x-axis limits
      minor_breaks = seq(3, 9, by = 1)   # Minor breaks for grid lines at each unit
    ) +
    scale_y_continuous(
      breaks = c(3, 4, 5, 6, 7, 8, 9),  # Custom y-axis ticks
      limits = c(3, 9),                 # Explicit y-axis limits
      minor_breaks = seq(3, 9, by = 1)   # Minor breaks for grid lines at each unit
    ) +
    labs(
      title = paste("Average number of edges at iteration",it),
      x = "Edges number on Apical Surface",
      y = "Edges number on Basal Surface"
    ) +
    #xlim(2.5, 9.5) +
    #ylim(2.5, 9.5) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "right"
    )
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

  perc<-data.frame(Iteration=integer(it+1), perc_sc=double(it+1))
  for (i in 1:(it+1)) {
    pointsAx<-filter(histpts, Iteration == i-1)$x
    pointsAy<-filter(histpts, Iteration == i-1)$y
    pointsBx <- pointsAx*ratio
    pointsBy <- pointsAy*ratio
    tslA<-deldir(pointsAx,pointsAy,rw=rect1)
    tilA<-tile.list(tslA)[(n+1):(2*n)]
    tslB<-deldir(pointsBx,pointsBy,rw=rect2)
    tilB<-tile.list(tslB)[(n+1):(2*n)]
    cellsdf<-data.frame(edgesA=integer(),edgesB=integer())
    for (j in 1:length(tilA)) {
      cellsdf[j,c(1,2)]<-c(length(tilA[[j]]$x),length(tilB[[j]]$x))
    }
    #percen<-length(filter(cellsdf, edgesA==edgesB)[[1]]) #numbers
    #perc[i,c(1,2)]<-c(i,percen)
    
    # Count the number of matching edges (where edgesA == edgesB)
    matching_tiles <- length(filter(cellsdf, edgesA == edgesB)[[1]])
    
    # Calculate the total number of tiles in this iteration
    total_tiles <- nrow(cellsdf)
    
    # Calculate the percentage of matching tiles
    perc_sc <- 100-(matching_tiles / total_tiles) * 100
    
    # Store the percentage for this iteration
    perc[i, c(1, 2)] <- c(i-1, perc_sc)
  }
  
  if(plotShow==TRUE){
    ploten<-ggplot(perc,aes(x=Iteration,y=perc_sc))+
      geom_line(colour="#F8766D")+
      xlab("Iteration of the algorithm")+
      ylab("Percentage of scutoids")+
      ylim(0, 100)+
      ggtitle("Evolution of the percentage of scutoids of the system")
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
  
  perc <- data.frame(Iteration = 0:it, percen = rep(0, (it+1)))
  
  # Use foreach for parallel execution
  foreach(j = 1:sim, .combine = 'rbind', .packages = c('deldir', 'dplyr'), .export = c('scutoids_percr')) %dopar% {
    
    # Fetch points for the current simulation
    percen <- scutoids_percr(results[[j, 1]]$points_evolution, n, it, ratio, ap_rad, cylen, plotShow = FALSE)$perc_sc
    
    return(percen)  # Return the perc_sc values for each simulation
  } -> perc_results
  
  # Stop the parallel cluster
  stopCluster(cl)
    # Summing the perc_results and averaging over simulations
  perc$percen <- colMeans(perc_results)
  
  # Plotting the results
  ploten <- ggplot(perc, aes(x = Iteration, y = percen)) +
    geom_line(colour = "#F8766D") +
    xlab("Iteration of the algorithm") +
    ylab("Percentage of scutoids") +
    ylim(0,100)+
    ggtitle("Evolution of the average percentage of scutoids")
  show(ploten)
  
  return(paste("Average percentage of scutoids: ",mean(perc_results)))
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

scutoids_analysis_simulations <- function(results, Ratio = 2.5, n = 100, sim = 100, it = 150, cylen,ap_rad,nsim=1){
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
    ptsx <- dplyr::filter(results[[i]]$points_evolution, Iteration==it)$x
    ptsy <- dplyr::filter(results[[i]]$points_evolution, Iteration==it)$y
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
  

  histdf_avcount <- dplyr::filter(histdf_avcount, avg_count > 1)

# 
#   scutoidsplot<-ggplot(histdf_avcount, aes(x = edgesA, y = edgesB, label =avg_count))+
#     geom_count(shape = "square", aes(color= avg_count,size=10))+
#     xlab("Average of edges on apical surface")+ylab("Average of edges on basal surface")+
#     ggtitle("Average polygon class of apical and basal surfaces")+
#     guides(colour = "colorbar", size = "none")+
#     labs(color = "Percentages")+
#     scale_size_area(max_size = 70)+
#     #geom_label()+
#     #scale_fill_gradient(low = "light blue", high = "deepskyblue")+
#     scale_color_viridis(option = "D",direction = -1, limits = c(0,35)) +
#     xlim(3,9)+
#     ylim(3,9)
#   # xlim(min(histdf_avcount$edgesA)-0.5, max(histdf_avcount$edgesA)-0.5)+
#   # ylim(min(histdf_avcount$edgesB)-0.5, max(histdf_avcount$edgesB)-0.5)
  
  scutoidsplot <- ggplot(histdf_avcount, aes(x = edgesA, y = edgesB, label = avg_count)) +
    geom_point(aes(color = avg_count), shape = 15, size=25) + # Square shape (15), size mapped to avg_count
    scale_color_viridis_c(
      option = "D",
      direction = -1,
      limits = c(0, 65),
      name = "Percentages"
    ) +
    scale_size_continuous(range = c(3, 9), guide = "none") + # Controlled size range
    scale_x_continuous(
      breaks = c(3, 4, 5, 6, 7, 8, 9),  # Custom x-axis ticks
      limits = c(3, 9.5),                 # Explicit x-axis limits
      minor_breaks = seq(3, 9, by = 1)   # Minor breaks for grid lines at each unit
    ) +
    scale_y_continuous(
      breaks = c(3, 4, 5, 6, 7, 8, 9),  # Custom y-axis ticks
      limits = c(3, 9),                 # Explicit y-axis limits
      minor_breaks = seq(3, 9, by = 1)   # Minor breaks for grid lines at each unit
    ) +
    labs(
      title = paste("Average number of edges at iteration",it),
      x = "Edges number on Apical Surface",
      y = "Edges number on Basal Surface"
    ) +
    #xlim(2.5, 9.5) +
    #ylim(2.5, 9.5) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "right"
    )
  show(scutoidsplot)
  
  output_file <- paste("results/figures/Edges at iteration",it, "_Average",nsim,"sim.")

  # Save the plot
  ggsave(
    filename = output_file,    # File name
    plot = scutoidsplot,                  # ggplot object
    device = "svg",            # File format
    width = 8,                 # Width in inches
    height = 6,                # Height in inches
    dpi = 300                  # Resolution (optional for SVG)
  )
  
  # Confirmation
  print(paste0("Plot saved:", output_file))
}









# ADDED 


energy_iteration <- function(pointsx, pointsy, n=100, Lay = 10,rad_coef,Radius,cylen,lambda,gamma,s0){
  
  
  
  tesener <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                         Elastic_energy = double(Lay),
                         Tension_energy = double(Lay),
                         Contractile_energy = double(Lay),
                         Bending_energy = double(Lay),
                         Total_energy=double(Lay))
  elasticener<-0
  tensionener<-0
  contractilener<-0
  totalener<-0
  
  A0 <- (2*pi*Radius*(cylen))/n
  gamad <- gamma
  lamad <- lambda
  
  cyl_thickness <- rad_coef*Radius - Radius
  rad <- list()
  rec <- list()
  for(k in 1:Lay){
    rad[[k]]<- Radius+(k-1)*(cyl_thickness/(Lay-1)) #the radius of the layer k
    rec[[k]]<-c(0,0+3*(2*pi*rad[[k]]),0,cylen)
  }
  
  for (i in 1:Lay) {
    tesel <- deldir(pointsx*(rad[[i]]/rad[[1]]),pointsy,rw=rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    
    perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
    areas <- sapply(tilest,function(x){x$area/A0})
    
    elener <- (sum((areas-1)^2))
    tenener <- sum(lamad*perims)
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/s0 )   
    contener <- sum((gam/2)*(perims^2))
    elasticener <- elasticener+elener
    tensionener <- tensionener+tenener
    contractilener <- contractilener+contener
    totalener<- totalener+elener+tenener+contener
  }
  elasticener <- elasticener/Lay
  tensionener <- tensionener/Lay
  contractilener <- contractilener/Lay
  totalener <- totalener/Lay
  
  return(c(elasticener,tensionener,contractilener,totalener))
}

energy_analisis_1sim <- function(histpts, it = 150, lay = 10, n =100,rad_coef,Radius,cylen,lambda,gamma,s0){
  
  histener<-data.frame(it = 1:(it+1), elen = double(it+1), tenen = double(it+1),
                       conten = double(it+1), toten = double(it+1))
  #print(head(histpts))
  for (i in 1:(it+1)) {
    pts <- filter(histpts, Iteration==(i-1))
    histener[i,c(2,3,4,5)] <- energy_iteration(pts$x,pts$y, n = n, Lay = lay,rad_coef = rad_coef,Radius = Radius,cylen = cylen,lambda=lambda,gamma=gamma,s0=s0)
  }
  
  return(histener)
}


energy_analisis_averages_par <- function(results, it = 150, lay = 10, n = 100,rad_coef,cylen,Radius,lambda,gamma,s0,nsim=1) {
  # Number of simulations
  N_SIM <- nrow(results)
  
  # Initialize a data frame to store averaged energy values
  average_energy <- data.frame(Iteration = 1:(it + 1), elen = double(it + 1), 
                               tenen = double(it + 1), conten = double(it + 1), 
                               toten = double(it + 1))
  
  # Set up parallel backend (adjust the number of cores as necessary)
  num_cores <- detectCores() - 1  # Leave one core free for other processes
  cl <- makeCluster(num_cores)
  registerDoParallel(cl) 
  
  # Use foreach for parallel execution
  energy_results <- 
    foreach(i = 1:N_SIM, .combine = 'rbind', .packages = c('dplyr','deldir'), .export = c('energy_analisis_1sim','energy_iteration')) %dopar% {
    histpts <- results[[i, 1]]$points_evolution
    energy_analisis_1sim(histpts, it = it, lay = lay, n = n, rad_coef = rad_coef,cylen = cylen,Radius=Radius,lambda=lambda,gamma = gamma,s0=s0)
  } 
  
  # Stop the parallel cluster
  stopCluster(cl)
  
  average_energy <- energy_results %>%
    group_by(it) %>%
    summarise(
      elen = mean(elen, na.rm = TRUE),
      tenen = mean(tenen, na.rm = TRUE),
      conten = mean(conten, na.rm = TRUE),
      toten = mean(toten, na.rm = TRUE)
   )
  
  # Summarize results
  average_energy <- energy_results %>%
    group_by(it) %>%
    summarize(across(c(elen, tenen, conten, toten), mean))
  
  average_energy$toten <- average_energy$toten 
  average_energy$tenen <- average_energy$tenen 
  average_energy$elen <- average_energy$elen 
  average_energy$conten <- average_energy$conten 
  
  # p<-ggplot(average_energy, aes(x = it))+
  #   geom_line(aes(y = toten, colour = "Total energy"))+
  #   geom_line(aes(y= tenen, colour = "Adhesion"))+
  #   geom_line(aes(y = elen, colour = "Elastic"))+
  #   geom_line(aes(y = conten, colour = "Contractility"))+
  #   # geom_line(aes(x = Layer, y = Bending_energy, colour = "Bending energy"))+
  #   labs(title = "Decomposition of system energies",
  #        x = "Iteration", y = "Average energy per cell",
  #        color = "Energy type") +
  #   scale_colour_manual("",
  #                       breaks = c("Total energy",
  #                                  "Adhesion",
  #                                  "Elastic",
  #                                  "Contractility"),
  #                       # "Bending energy"),
  #                       values = c("red",
  #                                  "blue",
  #                                  "darkgreen",
  #                                  "purple"))
  # 
  p <- ggplot(average_energy, aes(x = it)) +
    geom_line(aes(y = toten, colour = "Total"), linewidth = 1.2) +
    geom_line(aes(y = tenen, colour = "Adhesion"), linewidth = 1.2, linetype = "dashed") +
    geom_line(aes(y = elen, colour = "Elastic"), linewidth = 1.2, linetype = "dashed") +
    geom_line(aes(y = conten, colour = "Contractility"), linewidth = 1.2, linetype = "dashed") +
    labs(
      title = "Average tissue energy",
      #subtitle = "Evolution of energy components over iterations",
      x = "Iteration",
      y = "Average Total Energy",
      color = "Energy Type"
    ) +
    scale_colour_manual(
      "",
      breaks = c("Total", "Adhesion", "Elastic", "Contractility"),
      values = c("red", "blue", "darkgreen", "purple")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "top",
      legend.key.size = unit(1, "lines"))+
      ylim(0, max(average_energy$toten) + 0.1 * max(average_energy$toten))
    
  show(p)
  
  output_file <- paste0("results/figures/Energies_Average",N_SIM,"sim..svg")
  
  # Save the plot
  ggsave(
    filename = output_file,    # File name
    plot = p,                  # ggplot object
    device = "svg",            # File format
    width = 8,                 # Width in inches
    height = 6,                # Height in inches
    dpi = 300                  # Resolution (optional for SVG)
  )
  
  # Confirmation
  print(paste0("Plot saved:", output_file))
  
  return(average_energy)
}


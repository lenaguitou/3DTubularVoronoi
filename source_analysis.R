library(deldir)
library(ggplot2)
#library(ggvoronoi)
#library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)
library(tidyr)


rad_coef = parameters$ratio_rad
steps = parameters$n_steps
Radius = parameters$apical_rad #5/(2*pi) #params$radiusA
Radius2 <- rad_coef*Radius
cyl_thickness <- Radius2-Radius
cyl_length = parameters$cyl_length
cyl_width = Radius*(2*pi)
rec <- list()
rad <- list()
L = parameters$n_layers
n = parameters$n_cells
A0 <- ((Radius2+Radius)*pi*cyl_length)/n
xmin <- 0
ymin <- 0
xmax <- cyl_width
ymax <- cyl_length
for(k in 1:L){
  rad[[k]]<- Radius+(k-1)*(cyl_thickness/(L-1)) #the radius of the layer k
  rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
}
lamad = parameters$lambda
gamad = parameters$gamma


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

energy_analysis_nobend <- function(histpts, it=150, n=100, Lay=10){
  
  tesener <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                         Elastic_energy = double(Lay),
                         Tension_energy = double(Lay),
                         Contractile_energy = double(Lay),
                         Total_energy=double(Lay))
  tesener2 <- data.frame( Layer = seq(1, rad_coef, by = ((rad_coef-1)/(Lay-1))),
                          Elastic_energy = double(Lay),
                          Tension_energy = double(Lay),
                          Contractile_energy = double(Lay),
                          Total_energy=double(Lay))
  
  points <- filter(histpts, Frame == 1)
  pointslast <- filter(histpts, Frame == it)
  
  for (i in 1:Lay) {
    tesel <- deldir(points$x*(rad[[i]]/rad[[1]]),points$y,rw=rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    
    perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
    areas <- sapply(tilest,function(x){x$area/A0})
    
    elener <- (sum((areas-1)^2)/n)
    tenener <- sum(lamad*perims)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contener <- sum((gam/2)*(perims^2))/n
    tesener[i,c(2,3,4,5)] <- c(elener/(Lay/n),
                               tenener/(Lay/n),
                               contener/(Lay/n),
                               (elener+tenener+contener)/(Lay/n))
    
    tesellast <- deldir(pointslast$x*(rad[[i]]/rad[[1]]),pointslast$y,rw=rec[[i]])
    tilestlast <- tile.list(tesellast)[(n+1):(2*n)]
    
    perimsl <- (tilePerim(tilestlast)$perimeters)/sqrt(A0)
    areasl <- sapply(tilestlast,function(x){x$area/A0})
    
    elenerl <- (sum((areasl-1)^2)/n)
    tenenerl <- sum(lamad*perimsl)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contenerl <- sum((gam/2)*(perimsl^2))/n
    tesener2[i,c(2,3,4,5)] <- c(elenerl/(Lay/n),
                                tenenerl/(Lay/n),
                                contenerl/(Lay/n),
                                (elenerl+tenenerl+contenerl)/(Lay/n))
  }
  
  #tesener_long <- tesener %>%
  #  pivot_longer(cols = c(Tension_energy, Elastic_energy, Contractile_energy, Total_energy),
  #               names_to = "Energy_Type", values_to = "Energy_Value")
  
  tesener_long <- tesener %>%
    pivot_longer(
      cols = c(Total_energy,Tension_energy, Elastic_energy, Contractile_energy),
      names_to = "Energy_Type",
      values_to = "Energy_Value"
    ) %>%
    mutate(
      Energy_Type = case_when(
        Energy_Type == "Total_energy" ~ "Total Energy",
        Energy_Type == "Tension_energy" ~ "Adhesion",
        Energy_Type == "Elastic_energy" ~ "Elastic",
        Energy_Type == "Contractile_energy" ~ "Contractility",
        TRUE ~ Energy_Type # Default case if no match
      )
    )
  
  tesener_long <- tesener_long %>%
    mutate(
      Energy_Type = factor(
        Energy_Type,
        levels = c("Total Energy", "Adhesion", "Elastic", "Contractility") # Desired order
      )
    )

  p <- ggplot(tesener_long, aes(x = Layer, y = Energy_Value, colour = Energy_Type)) +
    geom_line()  +
    scale_color_manual(
      values = c(
        "Total Energy" = "red",       # Red for Total Energy
        "Adhesion" = "blue",          # Blue for Adhesion
        "Elastic" = "darkgreen",          # Green for Elastic
        "Contractility" = "purple"    # Purple for Contractility
      )
    ) +
    labs(title = paste0("Energy before ", it, " iterations of the algorithm"),
         x = "Layer ratio (R/Ra)", y = "Average energy of the cell",
         color = "Energy type")
  show(p)
  
    #tesener2_long <- tesener2 %>%
    #pivot_longer(cols = c(Total_energy,Tension_energy, Elastic_energy, Contractile_energy ),
    #            names_to = "Energy_Type", values_to = "Energy_Value")
    
    tesener2_long <- tesener2 %>%
      pivot_longer(
        cols = c(Total_energy,Tension_energy, Elastic_energy, Contractile_energy),
        names_to = "Energy_Type",
        values_to = "Energy_Value"
      ) %>%
      mutate(
        Energy_Type = case_when(
          Energy_Type == "Total_energy" ~ "Total Energy",
          Energy_Type == "Tension_energy" ~ "Adhesion",
          Energy_Type == "Elastic_energy" ~ "Elastic",
          Energy_Type == "Contractile_energy" ~ "Contractility",
          TRUE ~ Energy_Type # Default case if no match
        )
      )
    
    tesener2_long <- tesener2_long %>%
      mutate(
        Energy_Type = factor(
          Energy_Type,
          levels = c("Total Energy", "Adhesion", "Elastic", "Contractility") # Desired order
        )
      )
  
  p2 <- ggplot(tesener2_long, aes(x = Layer, y = Energy_Value, colour = Energy_Type)) +
    geom_line() +
    scale_color_manual(
      values = c(
        "Total Energy" = "red",       # Red for Total Energy
        "Adhesion" = "blue",          # Blue for Adhesion
        "Elastic" = "darkgreen",          # Green for Elastic
        "Contractility" = "purple"    # Purple for Contractility
      )
    ) +
    labs(title = paste0("Energy after ", it, " iterations of the algorithm"),
         x = "Layer ratio (R/Ra)", y = "Average energy of the cell",
         color = "Energy type")
  show(p2)
  
  
  
  
  return(list(tesener,tesener2))
}

energy_iteration <- function(pointsx, pointsy, n=100, Lay = 10){
  
  
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
  
  for (i in 1:Lay) {
    tesel <- deldir(pointsx*(rad[[i]]/rad[[1]]),pointsy,rw=rec[[i]])
    tilest <- tile.list(tesel)[(n+1):(2*n)]
    
    perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
    areas <- sapply(tilest,function(x){x$area/A0})
    
    elener <- (sum((areas-1)^2)/n)
    tenener <- sum(lamad*perims)/n
    gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
    contener <- sum((gam/2)*(perims^2))/n
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

plot_energyss<-function(enerhist){
  
  p<-ggplot(enerhist, aes(x = it))+
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
}

energy_layers_sim1 <- function(histpts, it = 150, Lay = 10, n = 100, rect = rec) {
  
  tesener <- data.frame(iter = 1:(it + 1))
  
  for (i in 1:(it + 1)) {
    pts <- filter(histpts, Frame == i)
    ener <- rep(0, Lay)
    
    layers_energy <- seq(1, Lay, by = 1)
    
    for (j in layers_energy) {
      ener[j] <- energy_1_layer(pts$x * (rad[[j]] / rad[[1]]), pts$y, rect[[j]],
                                i = j, Lay = Lay, n = n)
    }
    
    tesener[i, 2:(Lay + 1)] <- ener
  }
  
  colnames(tesener)[2:(Lay + 1)] <- paste0("Layer ", 1:Lay)
  
  return(tesener)
}


energy_1_layer <- function(pointsx, pointsy, rec, i, Lay =10, n=100){
  
  tesel <- deldir(pointsx,pointsy,rw=rec)
  tilest <- tile.list(tesel)[(n+1):(2*n)]
  
  perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
  areas <- sapply(tilest,function(x){(x$area/A0)-1})
  
  gam<-gamad*exp((1-(rad[[i]]/rad[[1]]))/1 )
  
  tesener<-sum((areas)^2+(gam/2)*(perims^2)+
                 lamad*perims)/(Lay*n)
  return(tesener)
}

plot_energy_decomp <- function(layers_hist) {
  df_largo <- pivot_longer(layers_hist, cols = starts_with("Layer "), names_to = "variable", values_to = "valor")
  
  ggplot(df_largo, aes(x = iter, y = valor, color = variable)) +
    geom_line() +
    labs(title="Energy relaxation by layer",x = "Iteration", y = "Total Energy", color = "Layer")
}



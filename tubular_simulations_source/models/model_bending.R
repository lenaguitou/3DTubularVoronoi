library(deldir)
#library(ggplot2)
#library(gganimate)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)

# VERSION WITH BENDING  
#
# The main difference of this algorithm respect the original one is in the
# MOVEMENTS of the cells, since each layer moves separately.
# The other difference is in the ENERGY COMPUTATION, since we add a bending
# energy component, so the formula changes (and the function to compute it).
#
# All of this causes a change in the storage of the points.
# In the original algorithm, we just store the position of the apical layer, 
# since we know the position of each cell center in a layer by projections. In 
# this case, we have to store the position of the cell center in EVERY layer,
# because we don't project anymore. 
# That is why insted of a dataframe for the points, here we use a list of 
# dataframes.
# The rest of the code is similar to the original algorithm and works the same.
# Here we document the differences, to see the rest go to N_cylinder_algorithm.R
# in the N-cylinder folder.


  bending_move_points<-function(pt, wid, len, rc, n_cells = 100, n_layers = 3, rad){
    
    #We choose a cell
    ind <-sample(1:n_cells,1)
    
    #Here we apply a movement to each layer. To do it efficiently we use lapply
    pt <- lapply(1:n_layers, function(i){
      
      ptinx <- pt[[i]]$x[[ind]]+rnorm(1,mean=0,sd=rc)
      ptiny <- pt[[i]]$y[[ind]]+rnorm(1,mean=0,sd=rc)
      s <- (rad[[i]]/rad[[1]])
      
      # Check if the movement is inside the cylinder limits
      while((ptinx < 0 || ptinx > (s*wid)) ||
            (ptiny < 0 || ptiny > len)){
        ptinx <- pt[[i]]$x[[ind]] + rnorm(1, mean=0, sd=rc)
        ptiny <- pt[[i]]$y[[ind]] + rnorm(1, mean=0, sd=rc)
      }
      
      #Replicate the movement in the 2 copies of the cylinder rectangle
      pt[[i]]$x[c(ind,ind+n_cells,ind+2*n_cells)] <- c(ptinx, ptinx+s*wid, ptinx+2*s*wid)
      pt[[i]]$y[c(ind,ind+n_cells,ind+2*n_cells)] <- ptiny
      pt[[i]]
    })
    return(pt)
  }
  
  bending_tesellation_energy_N <- function(points, A0, rec, rad, gamma, lambda,
                                           omega = 1, n_cells, n_layers = 3, s0_ratio=1){
    
    tesener <- numeric(n_layers)
    
    #We compute energy as in the original algorithm
    tesener <- sapply(1:n_layers, function(i){
      
      tesel <- deldir(points[[i]]$x,points[[i]]$y,rw=rec[[i]])
      tilest <- tile.list(tesel)[(n_cells+1):(2*n_cells)]
      
      perims <- (tilePerim(tilest)$perimeters)/sqrt(A0)
      areas <- sapply(tilest,function(x){x$area/A0})
      
      gam<-gamma*exp((1-(rad[[i]]/rad[[1]]))/s0_ratio)
      
      sum((areas-1)^2+(gam/2)*(perims^2)+lambda*perims)
    })
    
    #We add a bending energy component, depending on the angle that the center 
    # does with the center above and below itself. That is why the apical and 
    # basal centers don't have bending energy
    bendener <- sapply(2:(n_layers-1), function(i){
      
      #First we compute the angles of every cell with the scalar product
      angles <- sapply(1:n_cells, function(j){
        
        #we compute the position of a cell center in a layer of the 3d cylinder
        ptcentral <- c(points[[i]]$y[j],
                       rad[[i]]*cos((1/rad[[i]])*points[[i]]$x[j]),
                       rad[[i]]*sin((1/rad[[i]])*points[[i]]$x[j]))
        
        #We compute the position of the cell center bellow
        ptinf <- c(points[[i-1]]$y[j],
                   rad[[i-1]]*cos((1/rad[[i-1]])*points[[i-1]]$x[j]),
                   rad[[i-1]]*sin((1/rad[[i-1]])*points[[i-1]]$x[j]))
        
        #we compute the position of the cell center above
        ptsup <- c(points[[i+1]]$y[j],
                   rad[[i+1]]*cos((1/rad[[i+1]])*points[[i+1]]$x[j]),
                   rad[[i+1]]*sin((1/rad[[i+1]])*points[[i+1]]$x[j]))
        
        # we compute the scalar product and norm, with some adjustments so that 
        #the cos is between -1 and 1 (without pmin and pmax sometime we get a
        #slightly higher or lower result)
        vec1 <- ptsup - ptcentral
        vec2 <- ptinf - ptcentral
        
        #We use pmin and pmax to avoid errors in the arc-cosine computation
        v <- pmin(pmax(((vec1%*%vec2)[1,1])/
                       (norm(vec1,type = "2")*norm(vec2,type = "2")),-1.0),1.0)
        ang <- acos(v)        
        return(ang)
      })
      #we use the angles to apply the bending energy formula
      return(sum(omega*((angles-pi)^2)))
    })
    #we sum energies and return it, in this case the energy is normalized by layer.
    return((sum(tesener)+sum(bendener))/n_layers)
  }
  
  choice_metropolis<-function(delta,beta){
    
    p<-numeric(1)
    
    if(delta<=0){p<-1}else if(exp(-delta*beta)==Inf){p<-0}else{p<-exp(-delta*beta)}
    
    a=sample(c(0,1), size = 1, replace=TRUE, prob = c(1-p,p))
    
    return(a)
  }

  ggplotvor<-function(plotpoints,tit,wid = 5){
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
      labs(title = tit) # Title text
    show(pl)
  }
  
  areasideplots<-function(px,py,rect,tit1,tit2,A0){
    tsl<-deldir(px,py,rw=rect)
    til<-tile.list(tsl)[(n+1):(2*n)]
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
  
  nu_sq <- function(points, rec, n_cells=100){
    Lay <- length(points)
    teselap <- deldir(points[[1]]$x, points[[1]]$y, rw = rec[[1]])
    teselba <- deldir(points[[n_layers]]$x, points[[n_layers]]$y, rw = rec[[n_layers]])
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
  
  #comienza el programa
  
  
  metropolisad_ben<-function(seed = 666, n_steps = 250, n_cells = 100, n_layers = 5,
                         apical_rad = 5/(2*pi), ratio_rad = 2.5, cyl_length = 20,
                         gamma = 0.15, lambda = 0.04, s0_ratio = 1,
                         omega = 1, beta = 100){
    
    
    #We define our variables
    
    
    RadiusB <-ratio_rad*apical_rad
    cyl_width_A <- 2*pi*apical_rad
    cyl_width_B <- 2*pi*RadiusB
    
    cyl_thickness <- RadiusB-apical_rad
    
    #We define the vertices of the apical plane
    
    xmin <- 0
    xmax <- cyl_width_A
    ymin <- 0
    ymax <- cyl_length
    
    r <- cyl_width_A/n_cells #radius to make the moves
    Am <- ((apical_rad+RadiusB)*pi*cyl_length)/n_cells
    
    rec <- list()
    rad <- list()
    
    for(k in 1:n_layers){
      rad[[k]]<- apical_rad+(k-1)*(cyl_thickness/(n_layers-1)) #the radius of the layer k
      rec[[k]]<-c(xmin,xmin+3*(2*pi*rad[[k]]),ymin,ymax)
    }
    
    #Start, first iteration
    
    set.seed(seed)
    
    x1 <- runif(n_cells,xmin,xmax)
    y1 <- runif(n_cells,ymin,ymax)
    x <-c(x1,x1+cyl_width_A,x1+2*cyl_width_A)
    y <-c(y1,y1,y1)
    
    points <- vector(mode = "list", length = n_layers)
    points <- lapply(1:n_layers, function(i){data.frame(x = (rad[[i]]/rad[[1]])*x, y = y )})
    
    pointsinit <- points
    energytesel <- bending_tesellation_energy_N(points, Am, rec , rad,
                                        gamma, lambda, omega, n_cells, n_layers, s0_ratio)
    
    #We create the variables to store the results
    
    energhist <- data.frame(Frame=numeric(n_steps), energy=numeric(n_steps))
    energhist[1,c(1,2)] <- c(0,energytesel)
    
    histpts <- vector(mode="list", length = n_layers)
    
    for (i in 1:n_layers) {
      histpts[[i]] <- data.frame(x = double(3*n_cells*n_steps), y = double(3*n_cells*n_steps), Frame = double(3*n_cells*n_steps))
      histpts[[i]][1:(3*n_cells),c(1,2)] <- pointsinit[[i]]
      histpts[[i]][1:(3*n_cells),3] <- 1
    }

    #Start of the loop
    
    for (j in 1:n_steps) {
      for(l in 1:n_cells) {
        points2 <- bending_move_points(points, wid = cyl_width_A, len = cyl_length, rc = r, n_cells=n_cells, n_layers=n_layers, rad= rad)
        energytesel2 <- bending_tesellation_energy_N(points2, Am, rec, rad, gamma, lambda,omega,n_cells,n_layers,s0_ratio)
        c <- choice_metropolis(energytesel2-energytesel, beta)
        cond <- c==1
        if(cond){
          points <- points2
          energytesel <- energytesel2
        }
        gc()
      }
      for (i in 1:n_layers) {
        histpts[[i]][(j*3*n_cells+1):(j*3*n_cells+3*n_cells),c(1,2)] <- points[[i]]
        histpts[[i]][(j*3*n_cells+1):(j*3*n_cells+3*n_cells),3] <- j+1
      }
      energhist[j+1,c(1,2)] <- c(j,energytesel)
      gc()
    }
    #save(histpts, file = paste0("results_", i, ".Rds"))
    #nu2 <- nu_sq(points = points, rec = rec, n_cells = 100)
    return(list(points_evolution=histpts,
                energy_evolution=energhist))
  }
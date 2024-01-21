parameters_static <- list(
  seed = 42, # 
  steps = 5, # simulation steps
  n = 50, # number of cells
  L=4, # numbers of layers
  RadiusA = 5/(2*pi), # apical radius
  Ratio = 2.5, # ratio between apical and basal radius 
  cyl_length = 20, # cylinder length
  gamma_ad = 0.15, #  contractility
  lambda_ad = 0.04, # adhesion 
  # elastic K = 1 cause adimentionality
  beta = 100, #constant involved in jumps because of perturbations
  s0=1 # 
)

parameters_bending <- list(
  seed = 42, # 
  steps = 5, # simulation steps
  n = 50, # number of cells
  L=4, # numbers of layers
  RadiusA = 5/(2*pi), # apical radius
  Ratio = 2.5, # ratio between apical and basal radius 
  cyl_length = 20, # cylinder length
  gamma_ad = 0.15, #  contractility
  lambda_ad = 0.04, # adhesion 
  # elastic K = 1 cause adimentionality
  omega = 1, # Bending Constant (Omega)
  beta = 100, #constant involved in jumps because of perturbations
  s0=1 # 
)
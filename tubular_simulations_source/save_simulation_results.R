save_results <- function(simulation,save_parallel=FALSE) {
  
  # FUNCTION TO SAVE RESULTS
  # It receives a list simulation, which is the direct return of the function
  # perform_simulations with "return_params=TRUE" as input. This function is found
  # in tubular_simulations_source/model_execution.R. Otherwise, simulation is a
  # list containing:
  #   - simulation$results. The direct result of perform_simulation without
  #     "return_params=TRUE".
  #   - simulation$parameters. The parameters used by this function
  #   - simulation$algorithm. The algorithm used in the simulation (could be
  #     "static" or "bending").
  
  if (save_parallel == TRUE) {
    date <- format(Sys.Date(), "%y_%m_%d")
    if (!dir.exists(paste0("results/parallel/"))) {
      dir.create(paste0("results/parallel/"), recursive = TRUE)
    }
    filename <- paste0("results/parallel/","/results_simulation_", date, ".RData")
    
    # verify if file already exists
    counter <- 1
    while (file.exists(filename)) {
      filename <- paste0("results/parallel/", "/results_simulation_", date, "_", counter, ".RData")
      counter <- counter + 1
    }
    
    # Save
    saveRDS(simulation, file = filename)
    print(paste0("Results saved successfully in path: ", filename))
    return(NULL)
  }
  
  algorithm <- as.character(simulation$algorithm)
  date <- format(Sys.Date(), "%y_%m_%d")
  
  # create directories
  if (!dir.exists(paste0("results/", algorithm))) {
    dir.create(paste0("results/", algorithm), recursive = TRUE)
  }
  
  # set filename
  filename <- paste0("results/", algorithm, "/results_simulation_", date, ".RData")
  
  # verify if file already exists
  counter <- 1
  while (file.exists(filename)) {
    filename <- paste0("results/", algorithm, "/results_simulation_", date, "_", counter, ".RData")
    counter <- counter + 1
  }
  
  # Save
  saveRDS(simulation, file = filename)
}

load_results <- function(filename) {
  # FUNCTION TO RECOVER THE SIMULATION RESULTS
  # It recovers the simulation saved with the previous function save_results
  recovered_simulation <- readRDS(filename)
  return(recovered_simulation)
}


library(foreach)
library(doParallel)
source("tubular_simulations_source/models/model_static.R")
#source("tubular_simulations_source/models/model_bending.R")
source("tubular_simulations_source/save_simulation_results.R")

perform_simulations<-function(
    parameters,
    algorithm=TRUE,
    return_params=TRUE,
    save=TRUE,
    save_parallel=FALSE
    ){
  
  # Function that perform the simulation. Parameters:
  #   - parameters. List with the necessary parameters, if it is not given, it uses
  #     the parameters from the file 
  #     "tubular_simulations_source/models/model_parameters.R"
  #   - algorithm. String with the algorithm that we want to use. It can be 
  #     whether "static" or "bending", or else it will fail.
  #   - return_params. If it is set to TRUE, it returns a list with the 
  #     results, the parameters, and the algorithm used. If FALSE (by default),
  #     it will return just the results. It is convenient to set it to TRUE if 
  #     we want to store the results.
  #   - save. If save is set to TRUE, the results produced by the simulation
  #     will be automatically saved with the save_result functions located in
  #     file save_simulation_results
  
  if (algorithm=="static") {
    result_alg <- do.call(metropolisad, parameters)
  }
  #else if (algorithm=="bending") {
  #  result_alg <- do.call(metropolisad_ben, parameters)
  #}
  else{
    stop("Introduce a valid value for parameter algorithm: 
         'static' or 'bending'")
  }
  results_with_params = list(results=result_alg,
                             parameters=parameters,
                             algorithm=algorithm)
  #if (save == TRUE){
  save_results(results_with_params,save_parallel)
  #}
  if(return_params){
    return(results_with_params)
  }
  else{
    return(result_alg)
  }
}

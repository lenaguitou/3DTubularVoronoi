library(jsonlite)
library(foreach)
library(doParallel)
source("tubular_simulations_source/models/model_static.R")
source("tubular_simulations_source/models/model_bending.R")
source("tubular_simulations_source/models/model_parameters.R")

perform_simulations<-function(
    parameters=NULL,
    algorithm="static",
    return_params=FALSE
){
  
  if (algorithm=="static") {
    if (is.null(parameters)) {
      parameters = parameters_static
    }
    result_alg <- do.call(metropolisad, parameters)
  }
  else if (algorithm=="bending") {
    if (is.null(parameters)){
      parameters = parameters_bending
    }
    result_alg <- do.call(metropolisad_ben, parameters)
  }
  else{
    stop("Introduce a valid value for parameter algorithm: 
         'static' or 'bending'")
  }
  if(return_params){
    return(list(results=result_alg, parameters=parameters, algorithm=algorithm))
  }else{
    return(result_alg)
  }
}

#'@title  Basic Processing Tools for Enriching get_weather outputs
#'
#'
#' @description A wraper of Param_Radiation(), Param_Atmospheric(), and Param_Temperature() for evaluating get_weather() datasets.
#' @author Germano Costa Neto
#' @param weather.data data.frame. A get_weather output
#' @export

processWTH <- function(weather.data){

  # computing radiation paramters
  weather.data <-Param_Radiation(weather.data=weather.data, merge = T)
  # computing atmospheric paramters
  weather.data <-Param_Atmospheric(weather.data=weather.data, merge = T)
  weather.data <- Param_Temperature(weather.data=weather.data,merge=T)

  return(weather.data)
}



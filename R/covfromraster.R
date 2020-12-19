#'@title  Easily Extraction of Worldwide Raster-based Data
#'
#' @description A wrapper of \code{raster::getData()} for extracting variables from rasters based on latitude and loongitude of experiments.
#'
#' @author Germano Costa Neto
#'
#' @param covraster RasterLayer. A raster from which data is going to be extracted.
#' @param env.data data.frame. A \code{get_weather} output.
#' @param Latitude character, the name of latitude column in env.data
#' @param Longitude character, the name of longitude column in env.data
#' @param env.id character, the name of environmental identification column in env.data
#' @param name.out character, output name
#'
#' @return
#' The original dataframe with additional geographic/weather information.
#'
#' @details
#' TODO
#'
#' @examples
#' # TODO
#' # Extraction of clay soil content for Nairobi, Kenya
#'require(EnvRtype)
#'data("clay_5_15")
#'env.data = data.frame(LAT = -1.367, LON = 36.834, env = 'NAIROBI')
#'env.data = extract_GIS(covraster = clay_5_15,name.out = 'clay_5_15',env.data = env.data)
#'head(env.data)
#'
#' @importFrom raster extract merge
#' @importFrom sp proj4string CRS coordinates<- proj4string<- coordinates
#'
#' @export

extract_GIS <- function(covraster=NULL,Latitude=NULL, Longitude=NULL,env.data=NULL,env.id=NULL,name.out = NULL){

  if(is.null(name.out)) name.out = 'ALT'
  if(is.null(Latitude)) Latitude <- 'LAT'
  if(is.null(Longitude)) Longitude <-'LON'
  if(is.null(env.id)) env.id <- 'env'
  loc  = data.frame(x=env.data[,Longitude],y=env.data[,Latitude])
  env = env.data[,env.id]
  sp::coordinates(loc)= ~x+y
  sp::proj4string(loc) = sp::CRS("+proj=longlat +datum=WGS84") # acho interessante colocar essa informacao no help da funcao.

  for(i in 1:length(names(covraster))) env = cbind(env,data.frame(raster::extract(covraster[[i]], loc)))
  names(env)[-1] = names(covraster)
  env<-data.frame(unique(env))
  names(env)[!names(env) %in% env.id] <- name.out
  return(raster::merge(env,env.data,by=env.id))
}

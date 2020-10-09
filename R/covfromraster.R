#'@title  Easily Extraction of Worldwide Raster-based Data
#'
#' @description A wrapper of \code{raster::getData()} for extracting variables from rasters based on latitude and loongitude of experiments.
#'
#' @author Germano Costa Neto
#'
#' @param covraster RasterLayer. A raster from which data is going to be extracted.
#' @param env.data data.frame. A \code{get_weather} output.
#'
#' @return
#' The original dataframe with additional geographic/weather information.
#'
#' @details
#' TODO
#'
#' @examples
#' # TODO
# env.data = get_weather(lat = -13.05, lon = -56.05, country = 'BRA')
# srtm = raster::getData('worldclim', var='tmin', lat = -13.05, lon = -56.05, res = 0.5)
# env.data = Extract_GIS(covraster = srtm, env.data = env.data)
#'
#' @importFrom raster extract merge
#' @importFrom sp proj4string CRS coordinates<- proj4string<- coordinates
#'
#' @export

Extract_GIS <- function(covraster=NULL,Latitude=NULL, Longitude=NULL,env.data=NULL,env.id=NULL,name.out = NULL){
  
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

#RasterToCov <-function(raster=NULL,env.data=NULL, lon=NULL,
                     #    lat=NULL, env.id=NULL, .crs=NULL,.path=NUL,covname=NULL){

 # if(is.numeric(lon) | is.numeric(lat)) loc <- data.frame(x=lon,y=lat)
  #if(!is.numeric(lon) | !is.numeric(lat)) loc  <- data.frame(x=env.data[,lon],y=env.data[,lat])
  #if(is.null(env.id)) env <- paste0('env_',1:length(lat))
  #if(!is.null(env.data) | !is.null(env.id)) env <- env.data[,env.id]

  #coordinates(loc)= ~x+y
  #proj4string(loc) = CRS("+proj=longlat +datum=WGS84")

  #for(i in 1:length(names(raster))) env = cbind(env,data.frame(extract(raster[[i]], loc)))
  #names(env)[-1] = names(raster)
  #env<-data.frame(unique(env))
  #if(!is.null(covname)) names(env)[!names(env) %in% 'env'] <- covname
  #return(merge(env,env.data,by='env'))
#}

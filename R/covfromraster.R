#'@title  Easily Extraction of Worldwide Raster-based Data
#' @description Extracts variables from rasters.
#' @author Germano Costa Neto
#' @param covraster RasterLayer. A raster from which data is going to be extracted.
#' @param weather.data data.frame. A get_weather output
#' @param K_E list of envirotype-related kernels (n x n genotypes-environment).
#' If NULL, benchmarck genomic kernels are built.
#' @param K_G list of genomic enabled kernels (p x p genotypes)
#' @param Y data.frame contaning the following colunms: environemnt, genotype, trait value
#' @param model model structure for genomic predicion. It can be c('MM','MDs','E-MM','E-MDs'),
#' which MM (main effect model or Y=fixed + G) amd MDs (Y=fixed+G+GxE)
#' @param reaction boolean, inclusion of a reaction norm based GxE kernel (default = FALSE)
#' @param intercept.random boolean, inclusion of a genomic random intercepet (default = FALSE)
#' @importFrom raster extract merge
#' @importFrom sp proj4string CRS coordinates<- proj4string<- coordinates
#' @export

Extract_GIS <- function(covraster=NULL,weather.data=NULL){

  loc  = data.frame(x=weather.data[,'LON'],y=weather.data[,'LAT'])
  env = weather.data[,'env']
  sp::coordinates(loc)= ~x+y
  sp::proj4string(loc) = sp::CRS("+proj=longlat +datum=WGS84") # acho interessante colocar essa informacao no help da funcao.

  for(i in 1:length(names(covraster))) env = cbind(env,data.frame(raster::extract(covraster[[i]], loc)))
  names(env)[-1] = names(covraster)
  env<-data.frame(unique(env))
  names(env)[!names(env) %in% 'env'] <- 'ALT'
  return(raster::merge(env,weather.data,by='env'))
}
#RasterToCov <-function(raster=NULL,weather.data=NULL, lon=NULL,
                     #    lat=NULL, env.id=NULL, .crs=NULL,.path=NUL,covname=NULL){

 # if(is.numeric(lon) | is.numeric(lat)) loc <- data.frame(x=lon,y=lat)
  #if(!is.numeric(lon) | !is.numeric(lat)) loc  <- data.frame(x=weather.data[,lon],y=weather.data[,lat])
  #if(is.null(env.id)) env <- paste0('env_',1:length(lat))
  #if(!is.null(weather.data) | !is.null(env.id)) env <- weather.data[,env.id]

  #coordinates(loc)= ~x+y
  #proj4string(loc) = CRS("+proj=longlat +datum=WGS84")

  #for(i in 1:length(names(raster))) env = cbind(env,data.frame(extract(raster[[i]], loc)))
  #names(env)[-1] = names(raster)
  #env<-data.frame(unique(env))
  #if(!is.null(covname)) names(env)[!names(env) %in% 'env'] <- covname
  #return(merge(env,weather.data,by='env'))
#}

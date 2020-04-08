#'@title  Easily Extraction of Worldwide Raster-based Data
#'
#' @description Get multuple genomic and envirotype-informed kernels for bayesian genomic prediciton
#' @author Germano Costa Neto
#' @param K_E list of envirotype-related kernels (n x n genotypes-environment).
#' If NULL, benchmarck genomic kernels are built.
#' @param K_G list of genomic enabled kernels (p x p genotypes)
#' @param Y data.frame contaning the following colunms: environemnt, genotype, trait value
#' @param model model structure for genomic predicion. It can be c('MM','MDs','E-MM','E-MDs'),
#' which MM (main effect model or Y=fixed + G) amd MDs (Y=fixed+G+GxE)
#' @param reaction boolean, inclusion of a reaction norm based GxE kernel (default = FALSE)
#' @param intercept.random boolean, inclusion of a genomic random intercepet (default = FALSE)
#' @importFrom raster extract
#' @importFrom reshape2 merge


Extract_GIS <- function(covraster=NULL,reference=NULL, lon=NULL,
                        lat=NULL, env.id=NULL, .crs=NULL,.path=NUL,covname=NULL){
  require(raster)
  require(reshape2)
  if(is.null(lat)) lat <-'LAT'
  if(is.null(lon)) lon <- 'LON'
  if(is.null(env.id)) env.id <- 'env'

  loc  = data.frame(x=reference[,lon],y=reference[,lat])
  env = reference[,env.id]
  coordinates(loc)= ~x+y
  proj4string(loc) = CRS("+proj=longlat +datum=WGS84")

  for(i in 1:length(names(covraster))) env = cbind(env,data.frame(extract(covraster[[i]], loc)))
  names(env)[-1] = names(covraster)
  env<-data.frame(unique(env))
  if(!is.null(covname)) names(env)[!names(env) %in% 'env'] <- covname
  return(merge(env,reference,by='env'))
}
#RasterToCov <-function(raster=NULL,reference=NULL, lon=NULL,
                     #    lat=NULL, env.id=NULL, .crs=NULL,.path=NUL,covname=NULL){

 # if(is.numeric(lon) | is.numeric(lat)) loc <- data.frame(x=lon,y=lat)
  #if(!is.numeric(lon) | !is.numeric(lat)) loc  <- data.frame(x=reference[,lon],y=reference[,lat])
  #if(is.null(env.id)) env <- paste0('env_',1:length(lat))
  #if(!is.null(reference) | !is.null(env.id)) env <- reference[,env.id]

  #coordinates(loc)= ~x+y
  #proj4string(loc) = CRS("+proj=longlat +datum=WGS84")

  #for(i in 1:length(names(raster))) env = cbind(env,data.frame(extract(raster[[i]], loc)))
  #names(env)[-1] = names(raster)
  #env<-data.frame(unique(env))
  #if(!is.null(covname)) names(env)[!names(env) %in% 'env'] <- covname
  #return(merge(env,reference,by='env'))
#}

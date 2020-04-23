#'@title  Basic Processing Tools for Enriching get_weather outputs
#'
#'
#' @description Summarize get_weather() outputs based on environments and defined time intervals (e.g.,phenology)
#' @author Germano Costa Neto
#' @param x data.frame. A get_weather output
#' @param lon numeric. Longitude values in WGS84
#' @param lat numeric. Latitude values in WGS84
#' @param env.id character. Identification of the site or environment.
#' @param DOY numeric. Julian day. DOY column as benchmark (from get_weather)
#' @param download.ALT boolean. Default as TRUE, to download Altitude data from SRTM database of CGIAR.
#' @param country character. ID for country. Default is 'BRA'. For more detais see getData() from raster package.

#' @importFrom raster getData



processWTH <- function(x,lon=NULL,lat=NULL,env.id=NULL,DOY=NULL,download.ALT=TRUE,country=NULL){

  if(is.null(country)) country <- 'BRA'
  ALT <- 600
  if(isTRUE(download.ALT)){
    if (!require(raster)) install.packages("raster");require(raster)
    cat('------------------------------------------------ \n')
    cat('ATTENTION: This function requires internet access \n')
    cat('for more detailes see raster package \n')
    cat('------------------------------------------------  \n')

    # collecting altitude data from raster SRTM database
    srtm <- getData('alt', country="BRA",mask=TRUE)
    df<-Extract_GIS(covraster = srtm,reference = x,lon = lon,lat = lat ,env.id = env.id,covname = 'ALT')
    ALT <- 'ALT'
  }

  # computing radiation paramters
  df<-Param_Radiation(df = df,DOY=DOY,LAT=lat,merge = T)
  # computing atmospheric paramters
  df <-Param_Atmospheric(df = df,Alt=ALT,merge = T)

  df <- Param_Temperature(df=df,merge=T)

  return(df)
}



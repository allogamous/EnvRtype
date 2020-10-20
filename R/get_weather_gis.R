#'@title  Easily Collection of Worldwide Daily Weather Data
#'
#'
#' @description Imports daily-scale weather data from the NASA-POWER GIS and geographic data from SRTM database.
#'
#' @author Germano Costa Neto
#'
#' @param env.id vector (character or level). Identification of the site/environment (e.g. Piracicaba01).
#' @param lat vector (numeric). Latitude values of the site/environment (e.g. -13.05) in WGS84.
#' @param lon vector (numeric). Longitude values site/environment (e.g. -56.05) in WGS84.
#' @param variables.names vector (character). Name of the variables. Should be "T2M","T2M_MAX","T2M_MIN","PRECTOT", "WS2M","RH2M","T2MDEW", "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN", and/or "ALLSKY_TOA_SW_DWN". See Details for more information.
#' @param start.day vector (character). First date in which weather/geographic data should be collected (e.g. "2015-02-15").
#' @param end.day vector (character). Last date in which weather/geographic data should be collected (e.g. "2015-06-15").
#' @param dir.path character. Directory for the output. If not informed, the output will be saved in the current workind directory.
#' @param save bollean. If TRUE, save each environmental data.frame as .csv in dir.path.
#' @param temporal.scale character. Default = 'DAILY'. See \code{get_power()} function in nasapower package for more details.
#' @param country character. Country in which the lat and lon values are positioned (e.g. 'BRA'). A single country should be informed. Wrapper of raster::getData.
#'
#' @return A data.frame with selected \code{variable.names} collected from a \code{start.day} to a \code{end.day} at the informed \code{lat} and \code{lon}.
#'
#' @details
#' The available variables are:
#' \itemize{
#'  \item T2M: Temperature at 2 Meters
#'  \item T2M_MAX: Maximum Temperature at 2 Meters
#'  \item T2M_MIN: Minimum Temperature at 2 Meters
#'  \item PRECTOT: Precipitation
#'  \item WS2M: Wind Speed at 2 Meters
#'  \item RH2M: Relative Humidity at 2 Meters
#'  \item T2MDEW: Dew/Frost Point at 2 Meters
#'  \item ALLSKY_SFC_LW_DWN: Downward Thermal Infrared (Longwave) Radiative Flux
#'  \item ALLSKY_SFC_SW_DWN: All Sky Insolation Incident on a Horizontal Surface
#'  \item ALLSKY_TOA_SW_DWN: Top-of-atmosphere Insolation
#' }
#'
#' @examples
#' ## Temperature for a single location:
#' get_weather(env.id = "NM", lat = -13.05, lon = -56.05,
#'             start.day = "2015-02-15", end.day = "2015-06-15",
#'             variables.names = c("T2M"))
#'
#' ## All variables for two locations:
#' env = c("NM","SO")
#' lat = c(-13.05,-12.32); lon = c(-56.05,-55.42)
#' plant.date = c("2015-02-15",'2015-02-13')
#' harv.date = rep("2015-06-15", 2)
#' get_weather(env.id = env, lat = lat, lon = lon,
#'             start.day = plant.date, end.day = harv.date)
#'
#' @references
#' Sparks A (2019). _nasapower: NASA-POWER Data from R_. R package version 1.1.3, <URL:https://CRAN.R-project.org/package=nasapower>.
#'
#' @importFrom utils write.csv
#' @importFrom nasapower get_power
#' @importFrom plyr ldply
#' @importFrom utils install.packages
#'
#' @export

#----------------------------------------------------------------------------------------
# getting weather data from NASAPOWER GIS database
# adaptation from nansapower package's get_power function
#----------------------------------------------------------------------------------------

get_weather=function (env.id = NULL, lat = NULL, lon = NULL, start.day = NULL, 
                     end.day = NULL, variables.names = NULL, dir.path = NULL, 
                     save = FALSE, temporal.scale = "DAILY", country = NULL) 
{
  cat("------------------------------------------------ \n")
  cat("ATTENTION: This function requires internet access \n")
  cat("------------------------------------------------  \n")
  if (is.null(env.id)) {
    env.id <- paste0("env", seq_along(lat))
  }
  if (!(is.character(env.id) || is.factor(env.id))) {
    stop("env.id should be a vector of characters (e.g. 'env1') or factors")
  }
  if (!requireNamespace("nasapower", quietly = TRUE)) {
    utils::install.packages("nasapower")
  }
  if (!requireNamespace("plyr", quietly = TRUE)) {
    utils::install.packages("plyr")
  }
  if (is.null(dir.path)) {
    dir.path = getwd()
  }
  if (is.null(start.day)) {
    start.day <- Sys.Date() - 1000
    cat(paste0("start.day is NULL", "\n"))
    cat(paste0("matched as ", start.day, "\n"))
  }
  if (is.null(end.day)) {
    end.day <- start.day + 30
    cat(paste0("end.day is NULL", "\n"))
    cat(paste0("matched as ", end.day, "\n"))
  }
  if (is.null(variables.names)) {
    variables.names = c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOT", 
                        "WS2M", "RH2M", "T2MDEW", "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN", 
                        "ALLSKY_TOA_SW_DWN")
  }
  env.id = as.factor(env.id)
  .Ne = length(env.id)
  .C = vector(length = .Ne, "list")
  for (.E in 1:.Ne) {
    CL = data.frame(nasapower::get_power(community = "AG", 
                                         lonlat = c(lon[.E], lat[.E]), pars = variables.names, 
                                         dates = c(start.day[.E], end.day[.E]), temporal_average = temporal.scale))
    CL$daysFromStart = 1:nrow(CL)
    .C[[.E]] = CL
    names(.C)[[.E]] = env.id[.E]
    if (isTRUE(save)) {
      utils::write.csv(file = paste(env.id[.E], ".csv", 
                                    sep = ""), x = CL)
    }
    cat(paste0("Environment ", env.id[.E], " downloaded \n"))
  }
  names(.C) = env.id
  .C = plyr::ldply(.C)
  names(.C)[names(.C) %in% ".id"] = "env"
  df  = .C
  try(if (!is.null(country)) {
    if (!requireNamespace("raster")) 
      utils::install.packages("raster")
    srtm <- raster::getData("alt", country = country, mask = TRUE)
    df <- extract_GIS(covraster = srtm, env.data = df)
  })
  return(df)
}

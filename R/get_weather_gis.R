#==================================================================================================
# Title.    : Collecting daily weather from NASA POWER
# Author.   : G Costa-Neto
# Created at: 2020-12-30
# Updated at: 2023-11-21 
# Previous versions: EnvRtype::extract_GIS()
# 
#
# get_weather(), based on nasapower::get_power() from Sparks et al 2018
#==================================================================================================
#'
#'@title  Easily Collection of Worldwide Daily Weather Data.
#'
#'
#' @description Imports daily-scale weather data from the NASA-POWER GIS and geographic data from SRTM database.
#'
#' @author Germano Costa Neto and Giovanni Galli, modified by Tiago Olivoto
#'
#' @param env.id vector (character or level). Identifies  the site/environment (e.g. Piracicaba01).
#' @param lat vector (numeric). Latitude values of the site/environment (e.g. -13.05) in WGS84.
#' @param lon vector (numeric). Longitude values site/environment (e.g. -56.05) in WGS84.
#' @param variables.names vector (character). Name of the variables. Should be "T2M","T2M_MAX","T2M_MIN","PRECTOT", etc. See Details for more information.
#' @param start.day vector (character). First date in which weather/geographic data should be collected (e.g. "2015-02-15").
#' @param end.day vector (character). Last date in which weather/geographic data should be collected (e.g. "2015-06-15").
#' @param dir.path character. Directory for the output. If not informed, the output will be saved in the current working directory.
#' @param save bollean. If TRUE, save each environmental data.frame as .csv in dir.path.
#' @param parallel bollean. If TRUE, a parallel strategy is implemented. The
#'   vectors are split into chunks with `chunk_size` elements (30 by default)
#'   where the data is downloaded. The function then sleeps for `sleep` seconds
#'   to ensure the available number of requests per minute.
#'   (https://github.com/ropensci/nasapower/issues/57)
#' @param workers The number of processes in parallel. Defaults to 90% of
#'   available cores.
#' @param chunk_size The size of the chunks where the parallel strategy is
#'   implemented. Defaults to 29. Increasing this number may exceed the limit of
#'   queries per minute of the API.
#' @param sleep The time (in seconds) to sleep the function after each chunk has
#'   been downloaded. Defaults to 60. Decreasing this number may exceed the
#'   limit of queries per minute of the API.
#' @return A data.frame with selected \code{variable.names} collected from a \code{start.day} to a \code{end.day} at the informed \code{lat} and \code{lon}.
#'
#' @details
#' The available variables are (from NASAPOWER & Computed):
#' \itemize{
#'  \item T2M: Temperature at 2 Meters. C
#'  \item T2M_MAX: Maximum Temperature at 2 Meters, C
#'  \item T2M_MIN: Minimum Temperature at 2 Meters, C
#'  \item PRECTOT: Precipitation (mm)
#'  \item WS2M: Wind Speed at 2 Meters, m/s
#'  \item RH2M: Relative Humidity at 2 Meters, percentage
#'  \item QV2M: Specific Humidity, the ratio of the mass of water vapor to the total mass of air at 2 meters (kg water/kg total air).
#'  \item T2MDEW: Dew/Frost Point at 2 Meters, C
#'  \item ALLSKY_SFC_LW_DWN: Downward Thermal Infrared (Longwave) Radiative Flux
#'  \item ALLSKY_SFC_SW_DWN: All Sky Insolation Incident on a Horizontal Surface
#'  \item ALLSKY_SFC_SW_DNI All Sky Surface Shortwave Downward Direct Normal Irradiance
#'  \item ALLSKY_SFC_UVA: All Sky Surface Ultraviolet A (315nm-400nm)  Irradiance
#'  \item ALLSKY_SFC_UVB: All Sky Surface Ultraviolet B (280nm-315nm)  Irradiance
#'  \item ALLSKY_SFC_PAR_TOT: All Sky Surface Photosynthetically Active Radiation (PAR) Total
#'  \item FROST_DAYS: If it was a frost day (temperature less than 0C or 32F)
#'  \item GWETROOT: Root Zone Soil Wetness (layer from 0 to 100cm),ranging from 0 (water-free soil) to 1 (completely saturated soil).
#'  \item GWETTOP: Surface Soil Wetness (layer from 0 to 5cm),ranging from 0 (water-free soil) to 1 (completely saturated soil).
#'  \item EVPTRNS: Evapotranspiration energy flux at the surface of the earth.
#'  \item P-ETP: Computed difference between preciptation and evapotranspiration.
#'  \item VPD: Vapour pressure deficit (kPa), computed from T2MDEW,T2M_MAX and T2M_MIN
#'  \item N: Photoperiod (in h), computed from latitude and julian day (day of the year)
#'  \item n: number of sunny hours in the day for a clear sky (no clouds)
#'  \item RTA: Radiation on the top of the atmosphere, computed from latitude and julian day (day of the year)
#'  \item TH1: Temperature-Humidity Index by NRC 1971
#'  \item TH2: Temperature-Humidity Index by Yousef 1985
#'  \item PAR_TEMP: computed ratio between ALLSKY_SFC_PAR_TOT and  T2M
#' }
#'
#' @examples
#'\dontrun{
#' ## Temperature for a single location
#' get_weather(env.id = "NM", lat = -13.05, lon = -56.05,
#'             start.day = "2015-02-15", end.day = "2015-06-15",
#'             variables.names = c("T2M"))
#'
#' ## GWETROOT (Root Zone Soil Wetness) for a single location
#' get_weather(env.id = "NM", lat = -13.05, lon = -56.05,
#'             start.day = "2015-02-15", end.day = "2015-06-15",
#'             variables.names = c("GWETROOT"))
#'
#' ## All variables for two locations
#' env = c("NM","SO")
#' lat = c(-13.05,-12.32); lon = c(-56.05,-55.42)
#' plant.date = c("2015-02-15",'2015-02-13')
#' harv.date = rep("2015-06-15", 2)
#' get_weather(env.id = env, lat = lat, lon = lon,
#'             start.day = plant.date, end.day = harv.date)
#'}
#'
#' @references
#' Sparks A (2018). _nasapower: NASA-POWER Data from R_. R package version 1.1.3, <URL:https://CRAN.R-project.org/package=nasapower>.
#' Yousef M.K.(1985) _Stress Physiology in Livestock_, CRC Press, Boca Raton, FL
#' NRC (1971) _A guide to environmental research on animals_, Natl. Acad. Sci., Washington, DC (1971)
#' @importFrom utils write.csv
#' @importFrom nasapower get_power
#' @importFrom plyr ldply
#' @importFrom utils install.packages
#' @importFrom raster getData
#'
#' @export

#----------------------------------------------------------------------------------------
# getting weather data from NASAPOWER GIS database
# adaptation from nansapower package's get_power function
#----------------------------------------------------------------------------------------
get_weather <- function(env.id = NULL,
                        lat = NULL,
                        lon = NULL,
                        start.day = NULL,
                        end.day = NULL,
                        variables.names = NULL,
                        dir.path = NULL,
                        save = FALSE,
                        parallel = TRUE,
                        workers = NULL,
                        chunk_size = 29,
                        sleep = 60)
{
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    utils::install.packages("doParallel")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    utils::install.packages("doParallel")
  }

  if (!requireNamespace("foreach", quietly = TRUE)) {
    utils::install.packages("foreach")
  }

  ########### helpers
  # helper function to download the data
  get_helper <- function(lon, lat, variables.names, start.day, end.day, env.id, save){
    CL = data.frame(nasapower::get_power(community = "ag",
                                         lonlat = c(lon, lat),
                                         pars = variables.names,
                                         dates = c(start.day, end.day),
                                         temporal_api = "daily"))
    cor_rain_name = which(names(CL) %in% "PRECTOTCORR")
    names(CL)[cor_rain_name] = "PRECTOT"
    CL$daysFromStart = 1:nrow(CL)
    CL$env <- env.id
    CL <- CL[, c(which(colnames(CL) == "env"), which(colnames(CL) != "env"))]
    if (isTRUE(save)) {
      utils::write.csv(file = paste(env.id, ".csv",
                                    sep = ""), row.names = F, x = CL)
    }
    return(CL)
  }

  # got from https://github.com/TiagoOlivoto/pliman/blob/master/R/utilities.R
  sec_to_hms <- function(t){
    paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
          formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
          formatC(t %% 60, width = 2, format = "d", flag = "0"),
          sep = ":"
    )
  }

  progress <- function(min = 0,
                       max = 100,
                       leftd = "|",
                       rightd = "|",
                       char = "=",
                       style = 2,
                       width = getOption("width"),
                       time = Sys.time()){
    # Adapted from https://stackoverflow.com/a/26920123/15245107
    return(list(min = min,
                max = max,
                leftd = leftd,
                rightd = rightd,
                char = char,
                style = style,
                width = width,
                time = time))
  }

  run_progress <- function(pb,
                           actual,
                           text = "",
                           digits = 0,
                           sleep = 0){
    Sys.sleep(sleep)
    elapsed <- sec_to_hms(as.numeric(difftime(Sys.time(), pb$time, units = "secs")))
    temp <- switch(
      pb$style,
      list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd),
           text = paste(text, paste(pb$leftd, '%s%s', pb$right, sep = ""))),
      list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 6,
           text =  paste(text, paste(pb$leftd, '%s%s', pb$right, sep = ""), '% s%%')),
      list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 9,
           text = paste(text, paste(pb$leftd, '%s%s', pb$rightd, sep = ""), elapsed)),
      list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 15,
           text = paste(text, paste(pb$leftd, '%s%s', pb$rightd, sep = ""), '% s%%', elapsed))
    )
    step <- round(actual / pb$max * (pb$width - temp$extra))
    message(sprintf(temp$text,
                strrep(pb$char, step),
                strrep(' ', pb$width - step - temp$extra),
                round(actual / pb$max * 100, digits = digits)), "\r")
    if(actual == pb$max){
      message("\n")
    }
  }
  # split a vector into chunks with a given length
  split_chunk <- function(vec, length){
    split(vec, ceiling(seq_along(vec) / length))
  }




  message("---------------------------------------------------------------")
  message('get_weater() - Pulling Daily Weather Features from NASA POWER')
  message('Connecting to the API client')
  message('https://docs.ropensci.org/nasapower')
  message("---------------------------------------------------------------")

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
    message(paste0("start.day is NULL", "\n"))
    message(paste0("matched as ", start.day, "\n"))
  }
  if (is.null(end.day)) {
    end.day <- start.day + 30
    message(paste0("end.day is NULL", "\n"))
    message(paste0("matched as ", end.day, "\n"))
  }
  if (is.null(variables.names)) {
  # variables.names = c("T2M", "T2M_MAX", "T2M_MIN",
  #                      "PRECTOTCORR", "WS2M", "RH2M",
  #                      "T2MDEW", "ALLSKY_SFC_LW_DWN",
   #                     "ALLSKY_SFC_SW_DWN",
  #                      'GWETROOT','GWETPROF','GWETTOP',
   #                     'EVLAND','FROST_DAYS','FRSEAICE',
  #                      'FRSNO','SNODP','QV2M','EVPTRNS','CDD10','TO3')

    variables.names = c("T2M", "T2M_MAX", "T2M_MIN","T2MDEW",#'CDD10',#'TS',
                       "ALLSKY_SFC_LW_DWN",
                        "ALLSKY_SFC_SW_DWN",'ALLSKY_SFC_SW_DNI',
                        'ALLSKY_SFC_PAR_TOT',
                       'ALLSKY_SFC_UVA','ALLSKY_SFC_UVB',
                        "PRECTOTCORR",'EVPTRNS','QV2M', "RH2M",
                        'GWETROOT',#'GWETPROF',
                        'GWETTOP',
                    #    'EVLAND',
                        'FROST_DAYS',#'FRSEAICE',
                        #'FRSNO','SNODP',
                         "WS2M")
  }
  variables.names[grepl(variables.names, pattern = "PRECTOT")] = "PRECTOTCORR"
  envs_to_pull = length(unique(env.id))
  env.id = as.factor(env.id)

  envs_to_pull <- unique(env.id)
  startTime <- Sys.time()
  message(paste0('Start at...........................', format(startTime, "%a %b %d %X %Y"),'\n'))
  message(paste0('Number of Environmental Units ........',length( envs_to_pull)))
  message(paste0('Parallelization.....................[ ',ifelse(isTRUE(parallel),'x',''),' ]'))


  # sequential strategy
  # sleeps for 'sleep' s after 30 or more queries in 60 s
  if(parallel == FALSE){
    results <- list()
    pb <- progress(max = length(env.id), style = 4)
    init_time <- Sys.time()
    iter <- 0
    for (i in 1:length(env.id)) {
      iter <- iter + 1
      query_issue <- as.numeric(difftime(Sys.time(), init_time, units = "secs")) > 60
      if(iter >= 30 & query_issue){
        message("Waiting ", sleep, "s for a new query to the API.")
        Sys.sleep(sleep)
        iter <- 0
        init_time <- Sys.time()
      }
      results[[i]] <-
        get_helper(lon = lon[i],
                   lat = lat[i],
                   variables.names = variables.names,
                   start.day = start.day[i],
                   end.day = end.day[i],
                   env.id = env.id[i],
                   save = save)
      msg <- paste0("Env ", env.id[i], " (", i, "/", length(env.id), ") ", "downloaded")
      run_progress(pb, actual = i, text = msg)

    }
  }
  # parallel strategy
  # split the vectors into chunks with 'chunk_size' points
  # apply the function in parallel using 'workers' cores
  if(parallel == TRUE){
    env.id_par = split_chunk(env.id, length = chunk_size)
    lat_par =  split_chunk(lat, length = chunk_size)
    lon_par =  split_chunk(lon, length = chunk_size)
    start.day_par =  split_chunk(start.day, length = chunk_size)
    end.day_par =  split_chunk(end.day, length = chunk_size)
    nworkers <- ifelse(is.null(workers),
                       trunc(parallel::detectCores() * 0.9), workers)
    clust <- parallel::makeCluster(nworkers)
    message(paste0('Number of threads.....................',nworkers))

    on.exit(parallel::stopCluster(clust))
    results <- list()
    pb <- progress(max = length(env.id_par), style = 4)
    for (i in 1:length(env.id_par)) {
      env.id_par_tmp <- env.id_par[[i]]
      lat_par_tmp <- lat_par[[i]]
      lon_par_tmp <- lon_par[[i]]
      start.day_par_tmp <- start.day_par[[i]]
      end.day_par_tmp <- end.day_par[[i]]
      parallel::clusterExport(clust,
                              varlist = c("get_helper", "lat_par_tmp", "lon_par_tmp", "variables.names",
                                          "start.day_par_tmp", "end.day_par_tmp", "env.id_par_tmp"),
                              envir = environment())
      length_chunk <- length(env.id_par[[i]])
      temp <- parallel::parLapply(clust, 1:length_chunk, function(j) {
        get_helper(lon = lon_par_tmp[j],
                   lat = lat_par_tmp[j],
                   variables.names = variables.names,
                   start.day = start.day_par_tmp[j],
                   end.day = end.day_par_tmp[j],
                   env.id = env.id_par_tmp[j],
                   save = save)
      })
      results[[i]] <- plyr::ldply(temp)
      if(i < length(env.id_par)){
        message("Waiting ", sleep, "s for a new query to the API.")
      #  msg <- paste0("Chunk ", i, "/", length(env.id_par), " (",length_chunk,  " points) downloaded")
       # run_progress(pb, actual = i, text = msg)
        Sys.sleep(sleep)
      }
    }
 #   message("NASA POWER: Done!")
  }
  .output_weather <-  plyr::ldply(results)

  message(paste0('NASA POWER done! Thank you, Sparks et al 2018! \n'))


  variables.names[grepl(variables.names, pattern = "PRECTOTCORR")] = "PRECTOT"

  comp1 <-''

  if(isTRUE("PRECTOT" %in% variables.names) & isTRUE('EVPTRNS' %in% variables.names))
  {
   # message(paste0('Computing P-ETP......................'))

    .output_weather$P_ETP <- .output_weather$PRECTOT - .output_weather$EVPTRNS
    variables.names <- c( variables.names,'P_ETP')
    comp1 <- 'x'
  }

  message(paste0('Computing P-ETP.....................[ ',comp1,' ]'))
  comp1 <-''

  if(isTRUE("T2MDEW" %in% variables.names) & isTRUE("T2M_MIN" %in% variables.names) & isTRUE("T2M_MAX" %in% variables.names))
  {
    teten <- function(Temp) return(.61078*exp((17.27*Temp)/(Temp+237.3)))

    .output_weather$VPD <- ((teten(.output_weather$T2M_MIN) +teten(.output_weather$T2M_MAX))/2)-teten(.output_weather$T2MDEW)
    comp1 <-'x'
    variables.names <- c( variables.names,'VPD')
  }
  message(paste0('Computing VPD.......................[ ',comp1,' ]'))
  comp1 <-''
  ## computing global radiation and day lenght
  Ra_fun <- function(J,lat){

    deg2rad <- function(deg) (deg * pi) / (180)


    rlat = deg2rad(lat)
    fi = .409*sin((2*pi/365)*J-1.39)
    dr = 1+ .033*cos((2*pi/365)*J)
    ws = acos(-tan(rlat)*tan(fi))
    Ra = (1440/pi)*.0820*dr*(ws*sin(rlat)*sin(fi)+cos(rlat)*cos(fi)*sin(ws))
  #  N = (24/pi)*ws
    # Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahl, Hsin-i Wu and Robert M. Schoolfield, 1995.
    #A model comparison for daylength as a function of latitude and day of the year. Ecological Modeling 80:87-95
    P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860*(J-186)))))
    P <-  (sin(0.8333 * pi/180) + sin(lat * pi/180) * sin(P)) / (cos(lat * pi/180) * cos(P))
    P <- pmin(pmax(P, -1), 1)
    DL <- 24 - (24/pi) * acos(P)

    return(data.frame(Ra=Ra,N=DL))
  }

  comp1 <-''

  DL <-    Ra_fun(J =.output_weather$DOY,lat=.output_weather$LAT)
  .output_weather$N <- round(DL$N,3)
  .output_weather$RTA <- round(DL$Ra,3)
  variables.names <- c(variables.names,'N','RTA')

  if(isTRUE("RTA" %in% variables.names) & isTRUE('ALLSKY_SFC_SW_DNI' %in% variables.names))
  {
    .output_weather$n <-  round(abs(.output_weather$N*(.output_weather$ALLSKY_SFC_SW_DNI /.output_weather$RTA)),3)
    variables.names <- c(variables.names,'n')
    comp1 <-'x'
  }
  message(paste0('Computing N,n,RTA...................[ ',comp1,' ]'))


  comp1 <-''

  if(isTRUE('T2M' %in% variables.names)& isTRUE('RH2M' %in% variables.names))
  {
    # THI1 = [1.8 × (AT) + 32] − [0.55 − 0.0055 × (RH)] × [1.8 × (AT) − 26];
    .output_weather$TH1 <- (1.8*.output_weather$T2M + 32) - (0.55 - 0.0055* .output_weather$RH2M)*(1.8*.output_weather$T2M-26)

    variables.names <- c(variables.names,'TH1')
    comp1 <-'x'
  }

  message(paste0('Computing TH1.......................[ ',comp1,' ]'))


  comp1 <-''

  if(isTRUE("T2MDEW" %in% variables.names) & isTRUE('T2M' %in% variables.names))
  {
    # THI2 = AT + (0.36 × DP) + 41.2.
    .output_weather$TH2 <- .output_weather$T2M + (0.36*.output_weather$T2MDEW)+41.2

    variables.names <- c(variables.names,'TH2')
    comp1 <-'x'
  }

  message(paste0('Computing TH2.......................[ ',comp1,' ]'))


  comp1 <-''

  if(isTRUE("ALLSKY_SFC_PAR_TOT" %in% variables.names) & isTRUE('T2M' %in% variables.names))
  {
    # THI2 = AT + (0.36 × DP) + 41.2.
    .output_weather$PAR_TEMP <- .output_weather$ALLSKY_SFC_PAR_TOT/.output_weather$T2M

    variables.names <- c(variables.names,'PAR_TEMP')
    comp1 <-'x'
  }
  message(paste0('Computing PAR_TEMP..................[ ',comp1,' ]'))
  comp1 <-''



  id.var = which(names( .output_weather) %in% variables.names)
  id.env =  which(names( .output_weather) %in% c("env",'LON','LAT','YYYYMMDD',"DOY",'daysFromStart'))
  .output_weather[, c(id.var )][ .output_weather[,  c(id.var )] == -999] = NA
  .output_weather <- .output_weather[,c(id.env,id.var )]

 # ids <- which(names(.output_weather) %in% c("env",'LON','LAT','YYYYMMDD','daysFromStart',variable.names))
 #.output_weather <- .output_weather[,ids]
  endTime <- Sys.time()
  message(paste0('\nEnvironmental Units downloaded.......',  length(unique(  .output_weather$env)),'/',length( envs_to_pull)))
  message(paste0('Daily Weather Features...............',  length(variables.names)))
  message(paste0('Environmental data...................',  length(variables.names)* length(unique(  .output_weather$env)),'\n'))

  message(paste0('Done!..............................', format(endTime, "%a %b %d %X %Y")))
  message(paste0('Total Time.........................',  round(difftime(endTime,startTime,units = 'secs'),2),'s'))

  return(  .output_weather)

}










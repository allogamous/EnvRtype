#'@title  Easily Collection of Worldwide Daily Weather Data
#'
#'
#' @description Imports daily-scale weather data from the NASA-POWER GIS and geographic data from SRTM database.
#'
#' @author Germano Costa Neto, modified by Tiago Olivoto
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
#' @param country vector (character). Country in which the lat and lon values are positioned (e.g. 'BRA'). For USA continental area, use USA1. Wrapper of raster::getData.
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
#'\dontrun{
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
#'}
#'
#' @references
#' Sparks A (2018). _nasapower: NASA-POWER Data from R_. R package version 1.1.3, <URL:https://CRAN.R-project.org/package=nasapower>.
#'
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
                         temporal.scale = "daily",
                         country = NULL,
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
    cat(sprintf(temp$text,
                strrep(pb$char, step),
                strrep(' ', pb$width - step - temp$extra),
                round(actual / pb$max * 100, digits = digits)), "\r")
    if(actual == pb$max){
      cat("\n")
    }
  }
  # split a vector into chunks with a given length
  split_chunk <- function(vec, length){
    split(vec, ceiling(seq_along(vec) / length))
  }

  cat("------------------------------------------------ \n")
  cat("ATTENTION: This function requires internet access \n")
  cat("------------------------------------------------  \n")
  cat('Connecting to the NASA POWER API Client, Sparks et al 2018 \n')
  cat('https://docs.ropensci.org/nasapower \n')
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
    variables.names = c("T2M", "T2M_MAX", "T2M_MIN",
                        "PRECTOTCORR", "WS2M", "RH2M",
                        "T2MDEW", "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN")
  }
  variables.names[grepl(variables.names, pattern = "PRECTOT")] = "PRECTOTCORR"
  env.id = as.factor(env.id)
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
        msg <- paste0("Chunk ", i, "/", length(env.id_par), " (",length_chunk,  " points) downloaded")
        run_progress(pb, actual = i, text = msg)
        Sys.sleep(sleep)
      }
    }
    message("\nNASA POWER: Done!")
  }
  df <-  plyr::ldply(results)

  if(is.null(country))
  {

    cat("\n")
    cat('country = NULL ; Not possible to download ALT data. Please provide a country ID (e.g., USA1, BRA, FRA)')

    variables.names[grepl(variables.names, pattern = "PRECTOTCORR")] = "PRECTOT"
    ids = which(names(df) %in% variables.names)
    df[, ids][df[, ids] == -999] = NA
    return(df)
  }
  if (!is.null(country)) {


    if (!requireNamespace("raster"))
      utils::install.packages("raster")
    suppressMessages('raster')
    suppressWarnings('raster')

    cat("\n")
    cat('Connecting to https://biogeo.ucdavis.edu/data/ using Hijmans 2021')
    cat("------------------------------------------------  \n")

    unique_country <- unique(country)
    raster_alt <- vector(mode = 'list',length = length(unique_country))
    names(raster_alt) = unique_country

    for(n in 1:length(unique_country))
    {
      if(isTRUE(grepl(x = unique_country[n],pattern = 'USA')))
      {
        id_region = as.numeric(gsub(x = unique_country[n],pattern = 'USA',replacement = ''))
        raster_alt[[n]] <- suppressMessages(raster::getData("alt", country = 'USA',mask = TRUE)[[id_region]])
      }
      else
      {
        raster_alt[[n]] <-  suppressMessages(raster::getData("alt", country = unique_country[n], mask = TRUE))
      }
    }


    df_no_alt <- df
    df <- c()
    for(i in 1:length(env.id))
    {
      country_raster <-raster_alt[[which(names(raster_alt)%in% country[i])]]
      id <- which(df_no_alt$env %in% env.id[i])
      df <- rbind(df,extract_GIS(covraster = country_raster, env.data = df_no_alt[id,],name.out = 'ALT'))
    }

    variables.names <- c(variables.names,'ALT')
    message("\nSRTM: Done!")

    variables.names[grepl(variables.names, pattern = "PRECTOTCORR")] = "PRECTOT"
    ids = which(names(df) %in% variables.names)
    df[, ids][df[, ids] == -999] = NA
    return(df)
  }

}









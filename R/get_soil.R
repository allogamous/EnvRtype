#==================================================================================================
# Title.    : Collecting soil data from SoilGrids
# Author.   : G Costa-Neto, based on B. Monier codes
# Created at: 2022-09-15
# Updated at 2025-01-02 (envirotypeR version 0.1.3)
# Previous versions: -
# Current Version: 0.1.3 (envirotypeR)
#
# get_soil()
#==================================================================================================
#
#'@title  Easily Collection of Soil Features (250m res)
#'
#'
#' @description Imports soil features from SoilGrids data set
#'
#' @author Brandon Monier, modified by Germano Costa Neto
#'
#' @param env.id vector (character or level). Identification of the site/environment (e.g. Piracicaba01).
#' @param lat vector (numeric). Latitude values of the site/environment (e.g. -13.05) in WGS84.
#' @param lon vector (numeric). Longitude values site/environment (e.g. -56.05) in WGS84.
#' @param variables.names vector (character).
#' @param properties A SoilGrid property. Defaults to \code{clay}. Possible
#' @param verbose Report messages to console? Defaults to \code{FALSE}.
#'
#' @return A `data.frame` object
#'
#'
#' @examples
#'\dontrun{
#' ## Downloading clay and nitrogen for a given location
#' get_soil(env.id = "NM", lat = -13.05, lon = -56.05,variables.names = c('clay','nitrogen'))
#'
#' ## download all variables (default) for a single location
#' get_soil(env.id = "NM", lat = -13.05, lon = -56.05)
#'
#' ## All variables for two locations
#' env = c("NM","SO")
#' lat = c(-13.05,-12.32); lon = c(-56.05,-55.42)
#' get_soil(env.id = env, lat = lat, lon = lon)
#'}
#'
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom jsonlite fromJSON
#'
#' @references
#' Poggio et al (2021). _SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty_. European Geosciences Union, <https://doi.org/10.5194/soil-7-217-2021>.
#'
#' @importFrom utils write.csv
#' @importFrom nasapower get_power
#' @importFrom plyr ldply
#' @importFrom utils install.packages
#' @importFrom raster getData
#'
#----------------------------------------------------------------------------------------
# getting soil data from SoiGrids
#----------------------------------------------------------------------------------------


get_soil = function(env.id,
                    lat,
                    lon,
                    home.path = NULL,  # home directory. If null, getwd()
                    variables.names='clay',
                    #     output_name = NULL, # character sdenoting the name of the model under test
                    sleep = 20,
                    logfile=NULL) # if is TRUE, and n.core is null, then use detectCores() - 1)
{


  if(sleep < 15) sleep <- 15


  if (!requireNamespace("plyr", quietly = TRUE)) {
    utils::install.packages("plyr")
  }


  if (!requireNamespace("foreach", quietly = TRUE)) {
    utils::install.packages("foreach")
  }


  if (!requireNamespace("httr", quietly = TRUE)) {
    utils::install.packages("httr")
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    utils::install.packages("jsonlite")
  }

  library(jsonlite, include.only = c("fromJSON"))
  library(httr, include.only = c("content", "GET"))

  if(is.null(logfile)) logfile <-'get_soil_log.txt'


  # === Support functions =======================================================



  .features_id_creation <- function(.data)
  {
    .data <- data.frame(.data)
    .data <- data.frame(.data)
    .data$depth <- gsub(.data$depth,
                        pattern = "-", replacement = "_")
    .data$feature <- paste0(.data$property,
                            "|", .data$depth)
    return(.data)
  }

  ## General parameters (that can be integrated into funtion)
  configParams <- list(
    urlTemplate = "https://rest.isric.org/soilgrids/v2.0/properties/query?lon=%f&lat=%f%s&depth=0-5cm&depth=0-30cm&depth=5-15cm&depth=15-30cm&depth=30-60cm&depth=60-100cm&depth=100-200cm&value=Q0.05&value=Q0.5&value=Q0.95&value=mean&value=uncertainty",
    acceptedProps = c(
      "bdod",
      "cec",
      "clay",
      "cfvo",
      "nitrogen",
      "ocd",
      "ocs",
      "phh2o",
      "sand",
      "silt",
      "soc",
      "wv0010", "wv0033", "wv1500"
    )
  )

  formatSoilGrid <- function(
    json = NULL
  ) {
    if (is.null(json)) stop("Missing JSON data.")
    if (!is.list(json)) stop("This does not appear to be JSON data.")
    if (is.null(json$properties$layers)) stop("This does not appear to be SoilGrid data.")

    tmpProc <- lapply(seq_len(length(json$properties$layers$depths)), function(i) {
      tmpDf      <- json$properties$layers$depths[[i]]

      tmpDf <- cbind(tmpDf$label, tmpDf$values)

      tmpDf$name <- json$properties$layers$name[i]
      tmpDf$unit <- json$properties$layers$unit_measure$mapped_units[i]
      tmpDf$lon  <- json$geometry$coordinates[1]
      tmpDf$lat  <- json$geometry$coordinates[2]

      colnames(tmpDf) <- c(
        "depth", "q_5", "q_50", "q_95",
        "mean", "uncertainty", "property",
        "unit", "lon", "lat"
      )
      return(tmpDf)
    })

    return(do.call("rbind", tmpProc))
  }


  getSoilData <- function(
    lon = 0,
    lat = 0,
    properties = "clay",
    verbose = F
  ) {
    if (is.null(url)) stop("Missing URL signature.")

    if (any(!properties %in% configParams$acceptedProps)) stop("Incorrect properties.")

    url <- configParams$urlTemplate
    propString <- paste("&property=", properties, collapse = "", sep = "")

    finalURL <- sprintf(url, lon, lat, propString)

    if (verbose) message("Getting query...")

    getReq  <- httr::GET(finalURL)
    getReq2 <- httr::content(getReq, "text", encoding = "ISO-8859-1")

    jsonReq <- jsonlite::fromJSON(getReq2)

    if (verbose) message("Finished (", round(jsonReq$query_time_s, 3), "s)")

    return(formatSoilGrid(jsonReq))
  }


  getSoil_helper<-function(.env,.lat,.lon,.properties='clay',.sleep=10) {

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

    # Initialize a list to store the results
    n.env = length(.env)
    my_result <- vector("list", n.env)
    pb <- progress(max = n.env, style = 4)


    for (i in 1:n.env) {
      success <- FALSE
      while(!success) {
        tryCatch({
          my_result[[i]] <- suppressWarnings(getSoilData(lon = .lon[i], lat = .lat[i], properties = .properties))
          names(my_result)[i] <- .env[i]  # Changed from .env[i] to envs[i]
          success <- TRUE
          # message(paste0('Doing....',i,'/',n.env,'\n'))
        }, error = function(e) {
          message("Error encountered on iteration ", i, ": ", e$message, "\n")
          message("Waiting ", .sleep, "s for a new query to the API.")
          Sys.sleep(.sleep)
        })
      }
      msg <- paste0("Env ", .env[i], " (", i, "/",  n.env, ") ", "downloaded")
      run_progress(pb, actual = i, text = msg)
    }

    # renaming it

    # Return the list of results
    return(my_result)
  }

  # split a vector into chunks with a given length
  split_chunk <- function(vec, length){
    split(vec, ceiling(seq_along(vec) / length))
  }




  message("------------------------------------------------")
  message('get_soil() - Pulling Soil Features from SoilGrids')
  message('Connecting to the API client')
  message('https://soilgrids.org')
  message("------------------------------------------------  \n")



  if(is.null(home.path)) home.path = getwd()
  # if(is.null(   output.path))    output.path = getwd()

  if (is.null(variables.names)) {
    variables.names = c("bdod",
                        "cec",
                        "clay",
                        "nitrogen",
                        "phh2o",
                        "sand",   "wv0010", "wv0033", "wv1500",
                        "soc")

  }


  envs_to_pull <- unique(env.id)
  startTime <- Sys.time()
  message(paste0('Start at...........................', format(startTime, "%a %b %d %X %Y"),'\n'))
  message(paste0('Number of Environmental Units ........',length( envs_to_pull)))
  # message(paste0('Parallelization.....................[ ',ifelse(isTRUE(parallel),'x',''),' ]'))


  # sequential strategy
  results <- getSoil_helper(.lat = lat,.lon = lon,.env = env.id,.sleep = sleep,.properties = variables.names )
  message(paste0('Soil Grids: Done! Thank you, Poggio et al 2021 ! \n'))
  results <- plyr::ldply(  results,.id = 'environmental_unit')
  results <- .features_id_creation(.data =   results)
  results <- reshape2::melt(results ,id.var=c('environmental_unit','property','lat','lon','unit','uncertainty','depth','feature'))
  results <- reshape2::dcast(results,environmental_unit~feature+variable,value.var = 'value')


  # results_joint <- reshape2::dcast( plyr::ldply(results),env~feature,value.var = 'mean')


  # results <- reshape2::dcast( results ,env~feature,value.var = 'mean') #plyr::ldply(results)

  # property | unit | depth | stat (mean etc)


  #  output_soil_mean <- temp #reshape2::dcast(   temp,env~feature,value.var = 'mean')
  # output_soil_mean <- reshape2::dcast(   temp,env~feature,value.var = 'mean')

  endTime <- Sys.time()
  #  message(paste0('Environmental Units downloaded ......',  length(unique(  output_soil_mean$env)),'/',length( envs_to_pull)))
  #  message(paste0('Soil Features downloaded.............',  ncol(  output_soil_mean )-1),'\n')
  # message(paste0('Environmental data...................',  length(unique(  output_soil_mean$env))*(ncol(  output_soil_mean )-1)),'\n')

  message(paste0('Done!..............................', format(endTime, "%a %b %d %X %Y")))
  message(paste0('Total Time.........................',  round(difftime(endTime,startTime,units = 'secs'),3)))
  message("------------------------------------------------  \n")


  #  return( list(raw_data = results, table_data = results_joint ) )
  return(results)

}

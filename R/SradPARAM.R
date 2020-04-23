#'@title  Support functions to estimate basic radiation parameters
#'
#' @description Core of functions to estimate atmospheric parameters related to solar radiation. Fore more details about the equations see de FAO-evapotranspiration publication http://www.fao.org/3/X0490E/x0490e07.htm#radiation
#' @author Germano Costa Neto
#' @param K_E list of envirotype-related kernels (n x n genotypes-environment).
#' If NULL, benchmarck genomic kernels are built.
#' @param K_G list of genomic enabled kernels (p x p genotypes)
#' @param Y data.frame contaning the following colunms: environemnt, genotype, trait value
#' @param model model structure for genomic predicion. It can be c('MM','MDs','E-MM','E-MDs'),
#' which MM (main effect model or Y=fixed + G) amd MDs (Y=fixed+G+GxE)
#' @param reaction boolean, inclusion of a reaction norm based GxE kernel (default = FALSE)
#' @param intercept.random boolean, inclusion of a genomic random intercept (default = FALSE). For more details, see BGGE package vignette.
#' @importFrom BGGE getK
#' @importFrom stats median

# http://www.fao.org/3/X0490E/x0490e07.htm#radiation
Param_Radiation <-function(df,DOY=NULL, LAT=NULL,merge=FALSE){

  if(is.null(DOY)) DOY <-'DOY'
  if(is.null(LAT)) LAT <- 'LAT'

  DOY <- df[,DOY]
  LAT <- df[,LAT]

  Ra <- function(J,lat){
    rlat = deg2rad(lat)
    fi = .409*sin((2*pi/365)*J-1.39)
    dr = 1+ .033*cos((2*pi/365)*J)
    ws = acos(-tan(rlat)*tan(fi))
    Ra = (1440/pi)*.0820*dr*(ws*sin(rlat)*sin(fi)+cos(rlat)*cos(fi)*sin(ws))
    N = (24/pi)*ws
    cat('------------------------------------------------ \n')
    cat('Extraterrestrial radiation (RTA, MJ/m^2/day)  \n')
    cat('Daylight hours (N, hours) \n')

    return(data.frame(Ra=Ra,N=N))
  }


  RadN <-Ra(J = DOY,lat = LAT)
  LWD <- df$ALLSKY_SFC_LW_DWN
  DWN <- df$ALLSKY_SFC_SW_DWN
  LWD[LWD == -99] <- NA
  DWN[DWN == -99] <- NA
  LWD[is.na(LWD)] <- median(LWD,na.rm=T)
  DWN[is.na(DWN)] <- median(DWN,na.rm=T)

  Srad <- LWD-DWN
  n <- RadN$N*(Srad/RadN$Ra)

  cat('Actual duration of sunshine (n, hours) \n')
  cat('Solar Radiation (SRAD, MJ/m^2/day) \n')
  cat('------------------------------------------------ \n')
  cat('\n')

  if(isFALSE(merge)) return(data.frame(n=n,N=RadN$N,RTA=RadN$Ra, SRAD=Srad))
  if(isTRUE(merge)) return(data.frame(df,data.frame(n=n,N=RadN$N,RTA=RadN$Ra, SRAD=Srad)))

}








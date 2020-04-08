#'@title  Support functions to estimate atmospheric parameters
#'
#' @description Core of functions to estimate atmospheric parameters related to evapotranspiration. Fore more details about the equations see de FAO-evapotranspiration publication http://www.fao.org/3/X0490E/x0490e07.htm#atmospheric%20parameters
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



# http://www.fao.org/3/X0490E/x0490e07.htm#atmospheric%20parameters

Param_Atmospheric <- function(PREC=NULL,Tdew=NULL,
                            Tmin=NULL,Tmax=NULL,Alt=NULL,RH=NULL,
                            Rad=NULL,G=NULL,alpha=1.26,
                            df,merge=FALSE){


  teten <- function(Temp) return(.61078*exp((17.27*Temp)/(Temp+237.3)))
  psyco <- function(atm) return((((1.013*1E-3)*atm)/.622*2.45)*.665*atm*1E-3)
  AtmP <-  function(elevation=600) return(101.3*((293-.0065*elevation)/293)^5.26)
  slope.vapor <- function(Tmed) return(4098*(.6108*exp((17.27*Tmed)/(Tmed+237.3)))/(Tmed+237.2)^2)

  # Pristley-Taylor Equation
  EToPT<-function(alfa=1.26,Srad,G=NULL,slope,psyc){
    if(is.null(G)){G=0}
    W = slope/(slope+psyc)
    return(alfa*W*(Srad-G)*.408)}



  if(is.null(PREC)) PREC <-'PRECTOT'; PREC <- df[,PREC]
  if(is.null(Tdew)) Tdew <-'T2MDEW';Tdew <- df[,Tdew]
  if(is.null(Tmin)) Tmin <-'T2M_MIN';Tmin<- df[,Tmin]
  if(is.null(Tmax)) Tmax <-'T2M_MAX';Tmax <- df[,Tmax]
  if(is.null(Alt))  Alt <-'ALT'
  if(isTRUE(Alt  %in% names(df))) Alt <- df[,Alt]
  if(isFALSE(Alt  %in% names(df))){Alt <- 600; cat('Missing ALT value. We adopted 600m. Please use the Extract_GIS funciton to collect ALT from SRTM files \n')}

  if(is.null(RH))  RH <-'RH2M'; RH<- df[,RH]
  if(is.null(Rad)) Rad <- 'SRAD';Rad <- df[,Rad]


  Tmed <- (Tmin+Tmax)/2
  ATP <- AtmP(Alt)
  Es <- teten(Tdew)
  Eamin <- teten(Tmin)
  Eamax <- teten(Tmax)
  VPD <- ((Eamin +Eamax)/2)-Es
  #VPD <- Es*(1-RH/100)
  Psy <- psyco(atm = ATP)
  Slope <- 4098*(.6108*exp((17.27*Tmed)/(Tmed+237.3)))/(Tmed+237.2)^2
  ETo <- EToPT(alfa=alpha,Srad=Rad,G=G,slope=Slope,psyc=Psy)
  PETo <- PREC-ETo

  cat('------------------------------------------------ \n')
  cat('Slope of saturation vapour pressure curve (SPV, in kPa.Celsius) \n')
  cat('Vapour pressure deficit (VPD, in kPa) \n')
  cat('Potential Evapotranspiration (ETP, in mm.day) \n')
  cat('Deficit by Precipitation - ETP (PETP, in mm.day) \n')
  cat('------------------------------------------------ \n')
  cat('\n')

  if(isFALSE(merge)) return(data.frame(VPD=VPD,SPV=Slope,ETP=ETo,PETP=PETo))
  if(isTRUE(merge)) return(data.frame(df,data.frame(VPD=VPD,SPV=Slope,ETP=ETo,PETP=PETo)))


}








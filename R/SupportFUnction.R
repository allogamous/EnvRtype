#' @importFrom stats median

# http://www.fao.org/3/X0490E/x0490e07.htm#atmospheric%20parameters

AtmosphericParam<- function(PREC=NULL,Tdew=NULL,
                            Tmin=NULL,Tmax=NULL,Alt=600,RH=NULL,
                            Rad=NULL,G=NULL,alpha=1.26,
                            df,merge=FALSE){
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

  if(!isTRUE(merge)) return(data.frame(VPD=VPD,SPV=Slope,ETP=ETo,PETP=PETo))
  if(isTRUE(merge)) return(data.frame(df,data.frame(VPD=VPD,SPV=Slope,ETP=ETo,PETP=PETo)))


}

deg2rad <- function(deg) (deg * pi) / (180)


teten <- function(Temp) return(.61078*exp((17.27*Temp)/(Temp+237.3)))

psyco <- function(atm) return((((1.013*1E-3)*atm)/.622*2.45)*.665*atm*1E-3)
AtmP <-  function(elevation=600) return(101.3*((293-.0065*elevation)/293)^5.26)
slope.vapor <- function(Tmed) return(4098*(.6108*exp((17.27*Tmed)/(Tmed+237.3)))/(Tmed+237.2)^2)

# Pristley-Taylor Equation
EToPT<-function(alfa=1.26,Srad,G=NULL,slope,psyc){
  if(is.null(G)){G=0}
  W = slope/(slope+psyc)
  return(alfa*W*(Srad-G)*.408)}

Ra <- function(J,lat){
  rlat = deg2rad(lat)
  fi = .409*sin((2*pi/365)*J-1.39)
  dr = 1+ .033*cos((2*pi/365)*J)
  ws = acos(-tan(rlat)*tan(fi))
  Ra = (1440/pi)*.0820*dr*(ws*sin(rlat)*sin(fi)+cos(rlat)*cos(fi)*sin(ws))
  N = (24/pi)*ws
  cat('------------------------------------------------ \n')
  cat('Extraterrestrial radiation (Ra, MJ/m^2/day)  \n')
  cat('Daylight hours (N, hours) \n')

  return(data.frame(Ra=Ra,N=N))
}


# http://www.fao.org/3/X0490E/x0490e07.htm#radiation
SradParam<-function(df,DOY, LAT,merge=FALSE){


  RadN <-Ra(J = DOY,lat = LAT)
  LWD <- df$ALLSKY_SFC_LW_DWN
  DWN <- df$ALLSKY_SFC_SW_DWN
  LWD[LWD == -99] <- NA
  DWN[DWN == -99] <- NA
  LWD[is.na(LWD)] <- median(LWD,na.rm=TRUE)
  DWN[is.na(DWN)] <- median(DWN,na.rm=TRUE)

  Srad <- LWD-DWN
  n <- RadN$N*(Srad/RadN$Ra)

  cat('Actual duration of sunshine (n, hours) \n')
  cat('Solar Radiation (Srad, MJ/m^2/day) \n')
  cat('------------------------------------------------ \n')
  cat('\n')

  if(!isTRUE(merge)) return(data.frame(n=n,N=RadN$N,Ra=RadN$Ra, Srad=Srad))
  if(isTRUE(merge)) return(data.frame(df,data.frame(n=n,N=RadN$N,Ra=RadN$Ra, Srad=Srad)))

}



F.RUE.Temperature = function(Tbase1=9,Tbase2=45,Tmed,Topt1=26,Topt2=32){

  n = length(Tmed)
  FTAR = rep(1,n)
  for(i in 1:n){
    if(Tmed[i] < Topt1){FTAR[i] = (Tmed[i]-Tbase1)/(Topt1-Tbase1)}
    if(Tmed[i] > Topt2){FTAR[i] = (Tmed-Topt2)/(Tbase2-Topt2)}
  }
  return(FTAR)
}

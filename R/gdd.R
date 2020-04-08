#'@title Support functions to estimate thermal parameters
#'
#'
#' @description Calculation Of Heat According To The Growing Degree Day Model. From maximum (tmax) and minimum (tmin) temperature values, and according to the cardinal values of base temperature (Tbase) computes the daily thermal sum in ° C / day
#' @author Germano Costa Neto
#' @param Tmax numeric (vetor). Maximum air temperature in degree Celsius
#' @param Tmin numeric (vetor). Minimum air temperature in degree Celsius
#' @param Tbase1 numeric (vetor). cardinal value for temperature base for phenological development
#' @param Tbase2 numeric (vetor). maximum cardinal value for temperature base for phenological development


Param_Temperature <- function(df,Tmax=NULL, Tmin=NULL,
                              Tbase1=9,Tbase2=45,Tmed,Topt1=26,Topt2=32,merge=T) {

  if(is.null(Tmin)) Tmin <-'T2M_MIN';Tmin<- df[,Tmin]
  if(is.null(Tmax)) Tmax <-'T2M_MAX';Tmax <- df[,Tmax]

  Tmed <- (Tmin+Tmax)/2


  adjust_for_Tbase <- function(x, Tbase1) {
    ifelse(test = x < Tbase1, yes = Tbase1, no = x)
  }
  adjust_for_Tbase2 <- function(x, Tbase2) {
    ifelse(test = x > Tbase2, yes = Tbase2, no = x)
  }

  F.RUE.Temperature = function(Tbase1,Tbase2,Tmed,Topt1,Topt2){

    n = length(Tmed)
    FTAR = rep(1,n)
    for(i in 1:n){
      if(Tmed[i] < Topt1){FTAR[i] = (Tmed[i]-Tbase1)/(Topt1-Tbase1)}
      if(Tmed[i] > Topt2){FTAR[i] = (Tmed[i]-Topt2)/(Tbase2-Topt2)}
    }
    return(FTAR)
  }

  FRUE <- F.RUE.Temperature(Tbase1=Tbase1,Tmed=Tmed,Tbase2=Tbase2,Topt1=Topt1,Topt2=Topt2)

  tmax_adjusted <- adjust_for_Tbase(Tmax, Tbase1)
  tmin_adjusted <- adjust_for_Tbase(Tmin, Tbase1)

  tmax_adjusted <- adjust_for_Tbase2(tmax_adjusted, Tbase2)
  tmin_adjusted <- adjust_for_Tbase2(tmin_adjusted, Tbase2)

  gdd_temp <- (tmax_adjusted + tmin_adjusted) / 2 - Tbase1

  cat('------------------------------------------------ \n')
  cat('Growing Degree Day (GDD,oC/day) \n')
  cat('Effect of temperature on radiation use efficiency (FRUE, from 0 to 1)  \n')
  cat('Daily Temperature Range, (T2M_RANGE, oC day)\n')
  if(isFALSE(merge))   return(data.frame(GDD=gdd_temp,FRUE = FRUE, T2M_RANGE=Tmax-Tmin ))
  if(isTRUE(merge))   return(data.frame(df,GDD=gdd_temp,FRUE = FRUE, T2M_RANGE=Tmax-Tmin ))

}

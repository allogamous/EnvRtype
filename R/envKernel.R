#'@title Easily Building of Environmental Relatedness Kernels
#'
#'
#' @description Returns environmental kinships for reaction norm models. Output is a list containing the objects varCov kinship for environmental variables and envCov kinshp for environmental relatedness.
#' @author Germano Costa Neto
#' @param weather.data matrix. Data from environmental variables (or markers) per environment (or combinations of genotype-environment).
#' @param is.scaled boolean. If environmental data is mean-centered and scaled (default = TRUE), assuming x~N(0,1).
#' @param sd.tol numeric. Maximum standard deviation value for quality control. Coluns above this value are eliminated.
#' @param tol numeric. Value of tolerance (default = 0.001).
#' @param Y data.frame. Phenotypic data set containing environment id, genotype id and trait value.
#' @param merge boolean. if TRUE, the environmental covariables are merged with Y to build a n x n dimension env.kernel.
#' @param env.id character. Identification of experiment.
#' @param bydiag boolean. If TRUE, the parametrization by WW'/diag(WW') is applied. If FALSE, WW'/ncol(W).
#' @param gaussian boolean. If TRUE, uses the gaussian kernel parametrization for W, where envCov = exp(-h*d/q).
#' @param h.gaussian numeric. If gaussian = TRUE, returns the h parameter for exp(-h*d/q).
#' @importFrom stats sd dist
#' @export


EnvKernel <-function(weather.data,Y=NULL, is.scaled=T, sd.tol = 1,
                     tol=1E-3,bydiag=FALSE,merge=FALSE,
                     env.id=NULL,gaussian=FALSE, h.gaussian=NULL){

  nr<-nrow(weather.data)
  nc <-ncol(weather.data)

  if(!is.matrix(weather.data)){stop('weather.data must be a matrix')}
  if(isFALSE(is.scaled)){
    Amean <- weather.data-apply(weather.data,2,mean)+tol
    sdA   <- apply(Amean,2,sd)
    A. <- Amean/sdA
    removed <- names(sdA[sdA < sd.tol])
    weather.data <- A.[,!colnames(A.) %in% removed]
    t <- ncol(weather.data)
    r <- length(removed)
    cat(paste0('------------------------------------------------','\n'))
    cat(paste0(' Removed envirotype markers:','\n'))
    cat(paste0(r,' from ',t,'\n'))
    cat(paste0(removed,'\n'))
    cat(paste0('------------------------------------------------','\n'))

  }

  if(isTRUE(merge)){
    if(is.null(env.id)) env.id <- 'env'
    weather.data <- envK(weather.data = weather.data,df.pheno=Y,env.id=env.id)
  }
  if(isTRUE(gaussian)){
    O <- gaussian(x = weather.data,h=h.gaussian)
    H <- gaussian(x = t(weather.data),h=h.gaussian)
    return(list(varCov=H,envCov=O))

  }
  if(isFALSE(gaussian)){
    O <- tcrossprod(weather.data)#/ncol(weather.data)  # env.relatedness kernel from covariates
    H <- crossprod(weather.data)#/nrow(weather.data)   # covariable relatedness kernel from covariates
    if(isTRUE(bydiag)){
      O <- O/(sum(diag(O))/nc) + diag(1e-2, nrow(O))
      H <- H/(sum(diag(H))/nr) + diag(1e-2, nrow(H))
    }
    if(isFALSE(bydiag)){
      O <- O/nc + diag(1e-2, nrow(O))
      H <- H/nr + diag(1e-2, nrow(H))
    }


  }
  return(list(varCov=H,envCov=O))
}



gaussian <- function(x,h=NULL){
  d<-as.matrix(dist(x,upper = T,diag = T))^2
  q <- median(d)
  if(is.null(h)) h <- 1

  return(exp(-h*d/q))
}

envK = function(weather.data,df.pheno,skip=3,env.id){
  df.pheno <-data.frame(df.pheno)
  weather.data <-data.frame(weather.data)
  weather.data$env <- as.factor(rownames(weather.data))
  W <- as.matrix(merge(df.pheno,weather.data, by=env.id)[,-c(1:skip)])
  return(W)

}

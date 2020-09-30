#'@title Easily Building of Environmental Relatedness Kernels
#'
#'
#' @description Returns environmental kinships for reaction norm models. Output is a list containing the objects varCov kinship for environmental variables and envCov kinshp for environmental relatedness.
#' @author Germano Costa Neto
#' @param env.data matrix. Data from environmental variables (or markers) per environment (or combinations of genotype-environment).
#' @param is.scaled boolean. If environmental data is mean-centered and scaled (default = TRUE), assuming x~N(0,1).
#' @param sd.tol numeric. Maximum standard deviation value for quality control. Coluns above this value are eliminated.
#' @param tol numeric. Value of tolerance (default = 0.001).
#' @param Y data.frame. Phenotypic data set containing environment id, genotype id and trait value.
#' @param merge boolean. if TRUE, the environmental covariables are merged with Y to build a n x n dimension env.kernel.
#' @param env.id character. Identification of experiment.
#' @param Z_E matrix. NULL by default. is the model.matrix for environments (if merge = TRUE)
#' @param gaussian boolean. If TRUE, uses the gaussian kernel parametrization for W, where envCov = exp(-h*d/q).
#' @param h.gaussian numeric. If gaussian = TRUE, returns the h parameter for exp(-h*d/q).
#'
#' @return A list with environmental kinships for reaction norm models. Two matrices are produced. varCov with the distance for environmental covariables, and envCov with distances for genotypes.
#'
#' @details
#' TODO
#'
#' @examples
#' ### Loading the genomic, phenotype and weather data
#' data('maizeYield'); data("maizeWTH")
#'
#' ### getting the W matrix from weather data
#' W.cov <- W.matrix(env.data = maizeWTH)
#'
#' ### Parametrization by K_W = WW'/ncol(W)
#' EnvKernel(env.data = W.cov,
#'           Y = maizeYield,
#'           merge = FALSE,
#'           env.id = 'env',
#'           gaussian = FALSE)
#'
#' ### Parametrization by K_W = WW'/diag( WW')
#' EnvKernel(env.data = W.cov,
#'           Y = maizeYield,
#'           merge = FALSE,
#'           env.id = 'env',
#'           bydiag = TRUE)
#'
#' @seealso W.matrix
#'
#' @importFrom stats sd dist
#' @export


EnvKernel <-function(env.data,Y=NULL, is.scaled=TRUE, sd.tol = 1,
                     tol=1E-3, merge=FALSE,Z_E = NULL,
                     env.id='env',gaussian=FALSE, h.gaussian=NULL){
  
  nr<-nrow(env.data)
  nc <-ncol(env.data)
  
  GB_Kernel <-function(X,is.center=FALSE){
    if(isFALSE(is.center)) X = scale(x = X,center = T,scale = F)
    XXl <- X %*% t(X)
    K_G <- XXl/(sum(diag(XXl))/nrow(X)) + diag(1e-6, nrow(XXl))
    return(K_G)
  }
  
  if(!is.matrix(env.data)){stop('env.data must be a matrix')}
  if(isFALSE(is.scaled)){
    Amean <- env.data-apply(env.data,2,mean)+tol
    sdA   <- apply(Amean,2,sd)
    A. <- Amean/sdA
    removed <- names(sdA[sdA < sd.tol])
    env.data <- A.[,!colnames(A.) %in% removed]
    t <- ncol(env.data)
    r <- length(removed)
    cat(paste0('------------------------------------------------','\n'))
    cat(paste0(' Removed envirotype markers:','\n'))
    cat(paste0(r,' from ',t,'\n'))
    cat(paste0(removed,'\n'))
    cat(paste0('------------------------------------------------','\n'))
    
  }
  
  if(isTRUE(merge)){
    if(is.null(Z_E)){
      .DF = data.frame(Y[,env.id])
      names(.DF) ='env'
      Z_E = model.matrix(~0+env,.DF)
    }
    env.data = Z_E %*% env.data
  } 
  
  if(isTRUE(gaussian)){
    O <- gaussian(x = env.data,h=h.gaussian)
    H <- gaussian(x = t(env.data),h=h.gaussian)
  }
  if(isFALSE(gaussian)){
    
    O = GB_Kernel(env.data)
    H = GB_Kernel(t(env.data))
    
  }
  
 
  return(list(varCov=H,envCov=O))
}

gaussian <- function(x,h=NULL){
  d<-as.matrix(dist(x,upper = TRUE,diag = TRUE))^2
  q <- median(d)
  if(is.null(h)) h <- 1

  return(exp(-h*d/q))
}

envK = function(env.data,df.pheno,skip=3,env.id){
  df.pheno <-data.frame(df.pheno)
  env.data <-data.frame(env.data)
  env.data$env <- as.factor(rownames(env.data))
  W <- as.matrix(merge(df.pheno,env.data, by=env.id)[,-c(1:skip)])
  return(W)

}

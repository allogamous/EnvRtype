#'@title Easily Building of Environmental Relatedness Kernels
#'
#'
#' @description Returns environmental kinships for reaction norm models. Output is a list containing the objects varCov kinship for environmental variables and envCov kinshp for environmental relatedness
#' @author Germano Costa Neto
#' @param df.cov matrix. Data from environmental variables (or markers) per environment (or combinations of genotype-environment)
#' @param is.scaled boolean. If environmental data is mean-centered and scaled (default = TRUE), assuming x~N(0,1)
#' @param sd.tol numeric. Maximum standard deviation value for quality control. Coluns above this value are eliminated


EnvKernel <-function(df.cov,Y=NULL, is.scaled=T, sd.tol = 1,
                     tol=1E-3,bydiag=FALSE,merge=FALSE,
                     env.id=NULL){

  envK = function(df.cov,df.pheno,skip=3,env.id){
    df.pheno <-data.frame(df.pheno)
    df.cov <-data.frame(df.cov)
    df.cov$env <- as.factor(rownames(df.cov))
    W <- as.matrix(merge(df.pheno,df.cov, by=env.id)[,-c(1:skip)])
    return(W)
  }


  if(!is.matrix(df.cov)){stop('df.cov must be a matrix')}
  if(isFALSE(is.scaled)){
    Amean <- df.cov-apply(df.cov,2,mean)+tol
    sdA   <- apply(Amean,2,sd)
    A. <- Amean/sdA
    removed <- names(sdA[sdA < sd.tol])
    df.cov <- A.[,!colnames(A.) %in% removed]
    t <- ncol(df.cov)
    r <- length(removed)
    cat(paste0('------------------------------------------------','\n'))
    cat(paste0(' Removed envirotype markers:','\n'))
    cat(paste0(r,' from ',t,'\n'))
    cat(paste0(removed,'\n'))
    cat(paste0('------------------------------------------------','\n'))

  }

  if(isTRUE(merge)){
    if(is.null(env.id)) env.id <- 'env'
    df.cov <- envK(df.cov = df.cov,df.pheno=Y,env.id=env.id)
  }

  O <- tcrossprod(df.cov)#/ncol(df.cov)  # env.relatedness kernel from covariates
  H <- crossprod(df.cov)#/nrow(df.cov)   # covariable relatedness kernel from covariates
  O <- O + diag(1E-4,nrow=nrow(O),ncol=ncol(O))
  H <- H + diag(1E-4,nrow=nrow(H),ncol=ncol(H))
  if(isTRUE(bydiag)){
    O <- O/diag(O)
    H <- H/diag(H)
  }
  if(isFALSE(bydiag)){
    O <- O/ncol(df.cov)
    H <- H/nrow(df.cov)
  }


  return(list(varCov=H,envCov=O))
}


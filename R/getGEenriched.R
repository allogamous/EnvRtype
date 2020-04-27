#'@title  Envirotype-informed kernels for statistical models
#'
#' @description Get multiple genomic and/or envirotype-informed kernels for bayesian genomic prediciton.
#' @author Germano Costa Neto
#' @param K_E list. Contains nmatrices of envirotype-related kernels (n x n genotypes-environment). If NULL, benchmarck genomic kernels are built.
#' @param K_G list. Constains matrices of genomic enabled kernels (p x p genotypes). See BGGE::getK for more information.
#' @param Y data.frame. Should contain the following colunms: environemnt, genotype, phenotype.
#' @param model character. Model structure for genomic predicion. It can be \code{c('MM','MDs','E-MM','E-MDs')}, in which MM (main effect model \eqn{Y=fixed + G}) and MDs (\eqn{Y=fixed+G+GxE}).
#' @param reaction boolean. Indicates the inclusion of a reaction norm based GxE kernel (default = FALSE).
#' @param intercept.random boolean. Indicates the inclusion of a genomic random intercept (default = FALSE). For more details, see BGGE package vignette.
#' @param size_E character. \code{size_E=c('full','environment')}. In the first, 'full' means taht the environmental relationship kernel has the dimensions of n x n observations, which n = pq (p genotypes, q environments). If 'environment' the size of E-kernel is q x q.
#'
#' @return
#' A list of kernels (relationship matrices) to be used in genomic models.
#'
#' @details
#' TODO Define models.
#'
#' @examples
#' ### Loading the genomic, phenotype and weather data
#' data('maizeG'); data('maizeYield'); data("maizeWTH")
#'
#' ### Y = fixed + G
#' MM <- get_kernel(K_G = list(G = as.matrix(maizeG)),
#'                  Y = maizeYield, model = 'MM')
#' ### Y = fixed + G + GE
#' MDs <- get_kernel(K_G = list(G = as.matrix(maizeG)),
#'                   Y = maizeYield, model = 'MDs')
#'
#' ### Enriching models with weather data
#' W.cov <- W.matrix(weather.data = maizeWTH)
#' H <- EnvKernel(weather.data = W.cov, Y = maizeYield, merge = TRUE, env.id = 'env')
#'
#' EMM <- get_kernel(K_G = list(G = as.matrix(maizeG)),
#'                   Y = maizeYield,K_E = list(W = H$envCov),
#'                   model = 'EMM') # or model = MM
#'
#' ### Y = fixed + G + W + GE
#' EMDs <- get_kernel(K_G = list(G = as.matrix(maizeG)),
#'                    Y = maizeYield,
#'                    K_E = list(W = H$envCov),
#'                    model = 'MDs') # or model = MDs
#'
#' ### Y = fixed + W + G + GW
#' RN <- get_kernel(K_G = list(G = as.matrix(maizeG)),
#'                  Y = maizeYield,
#'                  K_E = list(W = H$envCov),
#'                  model = 'RNMM')
#'
#' ### Y = fixed + W + G + GW + GE
#' fullRN <- get_kernel(K_G = list(G = as.matrix(maizeG)),
#'                      Y = maizeYield,
#'                      K_E = list(W = H$envCov),
#'                      model = 'RNMDs')
#'
#' @seealso
#' BGGE::getk W.matrix
#'
#' @importFrom BGGE getK
#' @importFrom stats model.matrix
#'
#' @export

get_kernel <-function(K_E = NULL,                    #' environmental kernel
                      K_G,                           #' genotypic kernel (p x p genotypes)
                      Y,                             #' phenotypic dataframe
                      model = NULL,                  #' family model c('MM','MDs','EMM','EMDs','RNMM','RNMDs'),
                      intercept.random = FALSE,      #' insert genomic random intercept)
                      reaction = FALSE,
                      size_E = NULL                  #c('full','environment'),
                      ){
  #----------------------------------------------------------------------------
  # Start Step
  #  Y <- data.frame(env=Y[,env.id],gid=Y[,gen.id],value=Y[,trait.id])
  if (is.null(K_G))   stop('Missing the list of genomic kernels')
  if (!requireNamespace('BGGE')) utils::install.packages("BGGE")
  # if(!any(model %in% c("MM","MDs",'E-MM','E-MDs'))) stop("Model not specified. Choose between MM or MDs")
  if(is.null(model)) model <- 'MM'
  if(model == 'MM'){reaction <- FALSE; model_b <- 'MM';K_E=NULL}
  if(model == 'MDs'){reaction <- FALSE; model_b <-'MDs';K_E=NULL}
  if(model == 'EMM'){reaction <- FALSE; model_b <- 'MM'}
  if(model == 'EMDs'){reaction <- FALSE; model_b <-'MDs'}
  if(model == 'RNMM'){reaction <- TRUE; model_b <- 'MM'}
  if(model == 'RNMDs'){reaction <- TRUE; model_b <- 'MDs'}

  #----------------------------------------------------------------------------
  # getting genomic kernels (see BGGE)
  #----------------------------------------------------------------------------
  Y <- droplevels(Y)
  K = BGGE::getK(Y = Y, setKernel = K_G, model = model_b,intercept.random = intercept.random);
  names(K)   <- paste0('KG_',names(K))
  #----------------------------------------------------------------------------
  # If K_E is null, return benchmark genomic model
  #----------------------------------------------------------------------------
  if(is.null(K_E)){
    if(isFALSE(reaction)){
      cat("----------------------------------------------------- \n")
      cat('ATTENTION \n')
      cat('No K_E kernel was provided \n')
      cat('Environment effects assumed as fixed \n')
      cat("----------------------------------------------------- \n")
      cat(paste0('Model: ',model_b,'\n'))
      cat(paste0('Reaction Norm for E effects: ',FALSE,'\n'))
      cat(paste0('Reaction Norm for GxE effects: ',reaction,'\n'))
      cat(paste0('Intercept random: ',intercept.random,'\n'))
      cat(paste0("Kernels used: ",length(K),'\n'))
      cat("----------------------------------------------------- \n")
      return(K)
    }
    return(K)

  }
  #----------------------------------------------------------------------------
  # Envirotype-enriched models (for E effects)
  #----------------------------------------------------------------------------
  if(is.null(size_E)) size_E <- 'full'
  if(size_E == 'environment') for(q in 1:length(K_E)) K_E[[q]] <- EnvKernel(weather.data = K_E[[q]],Y = Y,merge = TRUE,env.id = 'env')$envCov


  h <- length(K_E);
  n <- length(K);

  K_e <- c()
  for(q in 1:h) K_e[[q]] = list(Kernel = K_E[[q]], Type = "D")
  names(K_e) <- paste0('KE_',names(K_E))


  K_f <- Map(c,c(K,K_e))

  #----------------------------------------------------------------------------
  # Envirotype-enriched models (for GE+E effects)
  #----------------------------------------------------------------------------
  if(isTRUE(reaction)){
    Zg <- stats::model.matrix(~0+gid,Y)
    ng <- length(K_G)
    Ng<-names(K_G)
    for(i in 1:ng) K_G[[i]] <-tcrossprod(Zg%*%K_G[[i]])
    ne <- length(K_E)
    A<-c()
    nome<-c()
    Ne<-names(K_E)
    for(g in 1:ng){for(e in 1:ne) {A <- cbind(A,list(K_G[[g]]*K_E[[e]])); nome <- c(nome,paste0('KGE_',Ng[g],Ne[e]))}}
    K_GE <- c()
    for(ge in 1:length(A)) K_GE[[ge]] <- list(Kernel=A[[ge]],Type='D')
    names(K_GE) <- nome
    K_f <- Map(c,c(K,K_e,K_GE))
  }

  if(isTRUE(intercept.random)) K_f<-K_f[-grep(names(K_f),pattern = 'GE_Gi')]
  #----------------------------------------------------------------------------
  # Reporting status
  #----------------------------------------------------------------------------

  cat("----------------------------------------------------- \n")
  cat(paste0('Model: ',model_b,'\n'))
  cat(paste0('Reaction Norm for E effects: ',TRUE,'\n'))
  cat(paste0('Reaction Norm for GxE effects: ',reaction,'\n'))
  cat(paste0('Intercept random: ',intercept.random,'\n'))
  cat(paste0("Total number of kernels: ",length(K_f),'\n'))
  cat("----------------------------------------------------- \n")
  return(K_f)
}


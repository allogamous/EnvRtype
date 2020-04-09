#'@title  Envirotype-informed kernels for statistical models
#'
#' @description Get multiple genomic and/or envirotype-informed kernels for bayesian genomic prediciton
#' @author Germano Costa Neto
#' @param K_E list of envirotype-related kernels (n x n genotypes-environment).
#' If NULL, benchmarck genomic kernels are built.
#' @param K_G list of genomic enabled kernels (p x p genotypes)
#' @param Y data.frame contaning the following colunms: environemnt, genotype, trait value
#' @param model model structure for genomic predicion. It can be c('MM','MDs','E-MM','E-MDs'),
#' which MM (main effect model or Y=fixed + G) amd MDs (Y=fixed+G+GxE)
#' @param reaction boolean, inclusion of a reaction norm based GxE kernel (default = FALSE)
#' @param intercept.random boolean, inclusion of a genomic random intercept (default = FALSE). For more details, see BGGE package vignette.
#' @param size_E character. size_E=c('full','environment'). In the first, 'full' means taht the environmental relationship kernel has the dimensions of n x n observations, which n = pq (p genotypes, q environments). If 'environment' the size of E-kernel is q x q.
#' @importFrom BGGE getK


#'----------------------------------------------------------------------------------------
get_kernel <-function(K_E = NULL,                    #' environmental kernel
                        K_G,                           #' genotypic kernel (p x p genotypes)
                        Y,                            #' phenotypic dataframe
                        model = NULL,  #' family model c('MM','MDs','E-MM','E-MDs'),
                        reaction=FALSE,
                        intercept.random = FALSE,      #' insert genomic random intercept)
                        size_E = NULL#c('full','environment'),
                        ){
  #'----------------------------------------------------------------------------
  #' Start Step
#  Y <- data.frame(env=Y[,env.id],gid=Y[,gen.id],value=Y[,trait.id])
  if (is.null(K_G))   stop('Missing the list of genomic kernels')
  if (!require(BGGE)) install.packages("BGGE");require(BGGE)
  if(!any(model %in% c("MM","MDs",'E-MM','E-MDs'))) stop("Model not specified. Choose between MM or MDs")
  if(is.null(model)) model <- 'MM'
  if(model == 'MM' | model =='E-MM')     model_b <- 'MM'
  if(model == 'MDs'| model == 'E-MDs')   model_b <- 'MDs'

  #'----------------------------------------------------------------------------
  #' getting genomic kernels (see BGGE)
  #'----------------------------------------------------------------------------
  Y <- droplevels(Y)
  K = getK(Y = Y, setKernel = K_G, model = model_b,intercept.random = intercept.random);
  names(K)   <- paste0('KG_',names(K))
  #'----------------------------------------------------------------------------
  #' If K_E is null, return benchmark genomic model
  #'----------------------------------------------------------------------------
  if(is.null(K_E)){
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
  #'----------------------------------------------------------------------------
  #' Envirotype-enriched models (for E effects)
  #'----------------------------------------------------------------------------
  if(is.null(size_E)) size_E <- 'full'
  if(size_E == 'environment') for(q in 1:length(K_E)) K_E[[q]] <- EnvKernel(df.cov = K_E[[q]],Y = Y,merge = T,env.id = 'env')$envCov


  h <- length(K_E);
  n <- length(K);

  K_e <- c()
  for(q in 1:h) K_e[[q]] = list(Kernel = K_E[[q]], Type = "D")
  names(K_e) <- paste0('KE_',names(K_E))


  K_f <- Map(c,c(K,K_e))

  #'----------------------------------------------------------------------------
  #' Envirotype-enriched models (for GE+E effects)
  #'----------------------------------------------------------------------------
  if(isTRUE(reaction)){
    Zg <- model.matrix(~0+gid,Y)
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
  #'----------------------------------------------------------------------------
  #' Reporting status
  #'----------------------------------------------------------------------------

  cat("----------------------------------------------------- \n")
  cat(paste0('Model: ',model_b,'\n'))
  cat(paste0('Reaction Norm for E effects: ',TRUE,'\n'))
  cat(paste0('Reaction Norm for GxE effects: ',reaction,'\n'))
  cat(paste0('Intercept random: ',intercept.random,'\n'))
  cat(paste0("Total number of kernels: ",length(K_f),'\n'))
  cat("----------------------------------------------------- \n")
  return(K_f)
}

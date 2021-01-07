#'@title  Envirotype-informed kernels for statistical models
#'
#' @description Get multiple genomic and/or envirotype-informed kernels for bayesian genomic prediciton.
#' @author Germano Costa Neto
#' @param K_E list. Contains nmatrices of envirotype-related kernels (n x n genotypes-environment). If NULL, benchmarck genomic kernels are built.
#' @param K_G list. Constains matrices of genomic enabled kernels (p x p genotypes). See BGGE::getK for more information.
#' @param gid character. denotes the name of the column respectively to genotypes
#' @param env character. denotes the name of the column respectively to environments
#' @param y character. denotes the name of the column respectively to phenotype values
#' @param data data.frame. Should contain the following colunms: environemnt, genotype, phenotype.
#' @param model character. Model structure for genomic predicion. It can be \code{c('MM','MDs','E-MM','E-MDs')}, in which MM (main effect model \eqn{Y=fixed + G}) and MDs (\eqn{Y=fixed+G+GxE}).
#' @param reaction boolean. Indicates the inclusion of a reaction norm based GxE kernel (default = FALSE).
#' @param intercept.random boolean. Indicates the inclusion of a genomic random intercept (default = FALSE). For more details, see BGGE package vignette.
#' @param dimension_KE character. \code{size_E=c('q','n')}. In the first, 'q' means taht the environmental relationship kernel has the dimensions of q x q observations,in which q is the number of environments. If 'n', the relationship kernel has the dimension n=pq, in which p is the number of genotypes
#' @param ne numeric. denotes the number of environments (q)
#' @param ng numeric. denotes the number of genotypes (p)
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
#'data("maizeYield") # toy set of phenotype data (grain yield per environment)
#'data("maizeG"    ) # toy set of genomic relationship for additive effects
#'data("maizeWTH")   # toy set of environmental data
#'y   = "value"      # name of the vector of phenotypes
#'gid = "gid"        # name of the vector of genotypes
#'env = "env"        # name of the vector of environments
#'ECs  = W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100, ], var.id = c("FRUE",'PETP',"SRAD","T2M_MAX"),statistic = 'mean')
#'## KG and KE might be a list of kernels
#'KE = list(W = env_kernel(env.data = ECs)[[2]])
#'KG = list(G=maizeG);
#'## Creating kernel models with get_kernel
#'## y = fixed + Genomic + error (main effect model, MM)
#'MM    = get_kernel(K_G = KG, y = y,gid = gid,env = env, data = maizeYield,model = "MM")
#'## y = fixed + Genomic + Genomic x Environment + error (MM plus a single GE deviation, MDs model, assuming GE as a block diagonal genomic effects)
#'MDs   = get_kernel(K_G = KG, y = y,gid = gid,env = env,  data = maizeYield, model = "MDs")
#'EMM   = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "EMM")
#'EMDs  = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "EMDs")
#'RMMM  = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMM")
#'RNMDs = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs")
#'## Examples of Models without any genetic relatedness, which G = I a identity for genotypes
#'MDs   = get_kernel(K_G = NULL, y = y,gid = gid,env = env,  data = maizeYield, model = "MDs")
#'EMDs  = get_kernel(K_G = NULL, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "EMDs")
#'
#' @seealso
#' env_typing W_matrix kernel_model
#'
#' @importFrom stats model.matrix
#'
#' @export


get_kernel <-function(K_E = NULL,                    #' environmental kernel ()
                      K_G = NULL,                    #' genotypic kernel (p x p genotypes)
                      data = NULL,                    #' phenotypic dataframe named after env, gid and trait
                      model = NULL,                  #' family model c('MM','MDs','EMM','EMDs','RNMM','RNMDs'),
                      intercept.random = FALSE,      #' insert genomic random intercept)
                      reaction = FALSE,              #' include reaction-norms (see model arguments)
                      dimension_KE = NULL,           #' k environments or n observations (n = pq)
                      ne = NULL,                     #' number of environments (calculated by default)
                      ng = NULL,                     #' number of genotypes (calculated by default)
                      env = 'env',
                      gid = 'gid',
                      y   = 'value'
                      ){
  #----------------------------------------------------------------------------
  # Start Step
  #
  #if (is.null(K_G))   stop('Missing the list of genomic kernels')

  if(is.null(model)) model <- 'MM'
  if(model == 'MM'   ){reaction <- FALSE; model_b <- 'MM';K_E=NULL}
  if(model == 'MDs'  ){reaction <- FALSE; model_b <-'MDs';K_E=NULL}
  if(model == 'EMM'  ){reaction <- FALSE; model_b <- 'MM'}
  if(model == 'EMDs' ){reaction <- FALSE; model_b <-'MDs'}
  if(model == 'RNMM' ){reaction <- TRUE; model_b <- 'MM'}
  if(model == 'RNMDs'){reaction <- TRUE; model_b <- 'MDs'}

  #----------------------------------------------------------------------------
  # getting genomic kernels
  #----------------------------------------------------------------------------
#  names(Y)[1:2] = c('env','gid')
 # Y <- droplevels(Y)

  Y <- data.frame(env=data[,env],gid=data[,gid],value=data[,y])

  if(is.null(ne)) ne = length(unique(Y$env))

  Zg <- stats::model.matrix(~0+gid,Y)
  colnames(Zg) = gsub(colnames(Zg),pattern = 'gid',replacement = '')
  ng <- length(unique(Y$gid))

  if(is.null(K_G)){
    cat("----------------------------------------------------- \n")
    cat('ATTENTION \n')
    cat(paste0('K_G = NULL, no K_G was provided','\n'))
    cat(paste0('Genomic effects assumed as an identity matrix Ig','\n'))
    K_G <- list(G = crossprod(Zg)/ne)
    }

  getK <- function(Y, X, kernel = c("GK", "GB"), setKernel = NULL, bandwidth = 1, model = c("SM", "MM", "MDs", "MDe"), quantil = 0.5,
                   intercept.random = FALSE)
  {
    #Force to data.frame
    Y <- as.data.frame(Y)

    Y[colnames(Y)[1:2]] <- lapply(Y[colnames(Y)[1:2]], factor)

    subjects <- levels(Y[,2])
    env <- levels(Y[,1])
    nEnv <- length(env)

    # check for repeated genotypes
    if(any(table(Y[,c(1:2)]) > 1))
      warning("There are repeated genotypes in some environment. They were kept")

    switch(model,
           'SM' = {
             if (nEnv > 1)
               stop("Single model choosen, but more than one environment is in the phenotype file")
             Zg <- model.matrix(~factor(Y[,2L]) - 1)
           },
           'Cov' = {
             Zg <- model.matrix(~factor(subjects) - 1)
           },{
             Ze <- model.matrix(~factor(Y[,1L]) - 1)
             Zg <- model.matrix(~factor(Y[,2L]) - 1)
           })

    if(is.null(setKernel)){
      if(is.null(rownames(X)))
        stop("Genotype names are missing")

      if (!all(subjects %in% rownames(X)))
        stop("Not all genotypes presents in the phenotypic file are in marker matrix")

      X <- X[subjects,]

      switch(kernel,
             'GB' = {
               # case 'G-BLUP'...
               ker.tmp <- tcrossprod(X) / ncol(X)
               #G <- list(Zg %*% tcrossprod(ker.tmp, Zg))
               G <- list(list(Kernel = Zg %*% tcrossprod(ker.tmp, Zg), Type = "D"))
             },
             'GK' = {
               # case 'GK'...
               D <- (as.matrix(dist(X))) ^ 2

               G <- list()
               for(i in 1:length(bandwidth)){
                 ker.tmp <- exp(-bandwidth[i] * D / quantile(D, quantil))
                 #G[[i]] <- Zg %*% tcrossprod(ker.tmp, Zg)
                 G[[i]] <- list(Kernel = Zg %*% tcrossprod(ker.tmp, Zg), Type = "D")
               }
             },
             {
               stop("kernel selected is not available. Please choose one method available or make available other kernel through argument K")
             })

      names(G) <- seq(length(G))

    }else{
      ## check kernels
      nullNames <- sapply(setKernel, function(x) any(sapply(dimnames(x), is.null)))
      if(any(nullNames))
        stop("Genotype names are missing in some kernel")

      # Condition to check if all genotype names are compatible
      equalNames <- sapply(setKernel, function(x) mapply(function(z, y) all(z %in% y), z=list(subjects), y=dimnames(x)) )
      if(!all(equalNames))
        stop("Not all genotypes presents in phenotypic file are in the kernel matrix.
             Please check dimnames")

      K <- lapply(setKernel, function(x) x[subjects, subjects]) # reordering kernel

      ker.tmp <- K
      #G <- list(Zg %*% tcrossprod(ker.tmp, Zg))
      G <- lapply(ker.tmp, function(x) list(Kernel = Zg %*% tcrossprod(x, Zg), Type = "D") )


      # Setting names
      if(is.null(names(K))){
        names(G) <- seq(length(G))
      }else{
        names(G) <- names(setKernel)
      }
    }

    tmp.names <- names(G)
    names(G) <- if(length(G) > 1) paste("G", tmp.names, sep ="_") else "G"


    switch(model,

           'SM' = {
             out <- G
           },

           'MM' = {
             out <- G
           },

           'MDs' = {
             E <- tcrossprod(Ze)
             #GE <- Map('*', G, list(E))
             GE <- lapply(G, function(x) list(Kernel = x$Kernel * E, Type = "BD"))
             names(GE) <- if(length(G) > 1) paste("GE", tmp.names, sep ="_") else "GE"
             out <- c(G, GE)
           },

           'MDe' = {
             ZEE <- matrix(data = 0, nrow = nrow(Ze), ncol = ncol(Ze))

             out.tmp <- list()

             for(j in 1:length(G)){
               out.tmp <- c(out.tmp, lapply(1:nEnv, function(i){
                 ZEE[,i] <- Ze[,i]
                 ZEEZ <- ZEE %*% t(Ze)
                 #K3 <- G[[j]] * ZEEZ
                 K3 <- list(Kernel = G[[j]]$Kernel * ZEEZ, Type = "BD")
                 return(K3)
               }))
             }
             if(length(G) > 1){
               names(out.tmp) <- paste(rep(env, length(G)), rep(tmp.names, each = nEnv), sep = "_" )
             }else{
               names(out.tmp) <-  env
             }
             out <- c(G, out.tmp)
           }, #DEFAULT CASE
           {
             stop("Model selected is not available ")
           })

    if(intercept.random){
      Gi <- list(Kernel = Zg %*% tcrossprod(diag(length(subjects)), Zg), Type = "D")
      out <- c(out, list(Gi = Gi))
    }

    return(out)
  }

  K = getK(Y = Y, setKernel = K_G, model = model_b,intercept.random = intercept.random);
  names(K)   <- paste0('KG_',names(K))


  #----------------------------------------------------------------------------
  # If K_E is null, return benchmark genomic model
  #----------------------------------------------------------------------------
  if(is.null(K_E)){
    if(isFALSE(reaction)){
      cat("----------------------------------------------------- \n")
      cat('ATTENTION \n')
      cat('K_E = NULL, no K_E kernel was provided \n')
      #cat('Environment effects assumed as fixed \n')
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
  if(is.null(dimension_KE)) dimension_KE <- 'q'
  # main envirotype effect

  if(dimension_KE == 'q'){
    K_Em = list()
    for(q in 1:length(K_E)) K_Em[[q]] <- K_E[[q]] %x% matrix(1,ncol=ng,nrow = ng)

    h <- length(K_E);
    n <- length(K);
  }
  if(dimension_KE =='n') K_Em <- K_E

  K_e <- c()
  for(q in 1:h) K_e[[q]] = list(Kernel = K_Em[[q]], Type = "D")
  names(K_e) <- paste0('KE_',names(K_E))


  K_f <- Map(c,c(K,K_e))

  #----------------------------------------------------------------------------
  # Envirotype-enriched models (for GE+E effects)
  #----------------------------------------------------------------------------
  if(isTRUE(reaction)){
    if(dimension_KE == 'n'){
      Ng<-names(K_G)
      for(i in 1:ng) K_G[[i]] <- matrix(1,ncol=ne,nrow=ne) %x% K_G[[i]]#tcrossprod(Zg%*%K_G[[i]])
      ne <- length(K_E)
      A<-c()
      nome<-c()
      Ne = names(K_E)
      ng = length(K_G)
      for(g in 1:ng){for(e in 1:ne) {A <- cbind(A,list(K_G[[g]]*K_E[[e]])); nome <- c(nome,paste0('KGE_',Ng[g],Ne[e]))}}
      K_GE <- c()
      for(ge in 1:length(A)) K_GE[[ge]] <- list(Kernel=A[[ge]],Type='D')
      names(K_GE) <- nome
      K_f <- Map(c,c(K,K_e,K_GE))
    }
    if(dimension_KE == 'q'){
      Ng<-names(K_G)
      #   for(i in 1:ng) K_G[[i]] <- matrix(1,ncol=ne,nrow=ne) %x% K_G[[i]]#tcrossprod(Zg%*%K_G[[i]])
      #  ne <- length(K_E)
      A<-c()
      nome<-c()
      Ne<-names(K_E)
      ne = length(K_E)
      ng = length(K_G)
      for(g in 1:ng){for(e in 1:ne) {A <- cbind(A,list(K_E[[e]]%x%K_G[[g]])); nome <- c(nome,paste0('KGE_',Ng[g],Ne[e]))}}
      K_GE <- c()
      for(ge in 1:length(A)) K_GE[[ge]] <- list(Kernel=A[[ge]],Type='D')
      names(K_GE) <- nome
      K_f <- Map(c,c(K,K_e,K_GE))
    }

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

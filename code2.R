#home.dir <- getwd()
#data.dir <- 'C:/Users/pc/Documents/Germano/Crossa_baseline_WGP/data'
#setwd(data.dir)
################################################################################
################################################################################
# GBLUP KERNELS

seed    = 71621
ite     = 3000
burn    = 600
thin    = 2
nFold   = 30


# Penotypic data
pheno <- droplevels(readRDS('GY594'))

EH <- droplevels(readRDS('EH594'))
PH <- droplevels(readRDS('PH594'))
GY <- droplevels(readRDS('GY594'))

phenoGE <- list(EH=EH,PH=PH,GY=GY)

# Molecular data
M <- readRDS('M594')


# Env data 1: Medias por ambiente
W1 <- readRDS(file = "W1")

# Env data 2: Medias por ambiente x intervalo
W2 <- readRDS(file = "W2")

# Env data 3: medias por genotipo-ambiente
K_W3 <- readRDS('KW3')

# Env data 3: medias por genotipo-ambiente ots (GK)
#K_W4 <- readRDS('W_gid_GK')


EnvKernel <-function(env.data,Y=NULL, is.scaled=TRUE, sd.tol = 1,
                     tol=1E-3, bydiag=FALSE, merge=FALSE,
                     env.id=NULL,gaussian=FALSE, h.gaussian=NULL){
  
  nr<-nrow(env.data)
  nc <-ncol(env.data)
  
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
    if(is.null(env.id)) env.id <- 'env'
    env.data <- envK(env.data = env.data,df.pheno=Y,env.id=env.id)
  }
  if(isTRUE(gaussian)){
    O <- gaussian(x = env.data,h=h.gaussian)
    H <- gaussian(x = t(env.data),h=h.gaussian)
    return(list(varCov=H,envCov=O))
    
  }
  if(isFALSE(gaussian)){
    O <- tcrossprod(env.data)#/ncol(env.data)  # env.relatedness kernel from covariates
    H <- crossprod(env.data)#/nrow(env.data)   # covariable relatedness kernel from covariates
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

Vcomp.BGGE<-function(model,digits=4){
  K = model$K
  size = length(K)
  comps = data.frame(matrix(NA,ncol=3,nrow=size))
  VarE =  data.frame(matrix(NA,ncol=3,nrow=1))
  names(comps) = names(VarE) = c("K","Var","SD.var")
  for(k in 1:size){
    comps [k,1] = names(K)[k]
    comps [k,2] = round(K[[k]]$varu,digits   )
    comps [k,3] = round(K[[k]]$varu.sd,digits)
  }
  VarE  [1,1] = "Residual"
  VarE  [1,2] = round(model$varE, digits    )
  VarE  [1,3] = round(model$varE.sd,digits   )
  comps = rbind(comps,VarE)
  return(comps)
}

#require(devtools)
#devtools::install_github('allogamous/EnvRtype')
#require(EnvRtype)
#install_github('cran/foreach')
require(foreach)

BGGE <- function(y, K, XF = NULL, ne, ite = 1000, burn = 200, thin = 3, verbose = FALSE, tol = 1e-10, R2 = 0.5) {
  
  ### PART I  - Conditional distributions functions and eigen descomposition ####
  # Conditional Distribution of tranformed genetic effects b (U'u)
  #' @import stats
  dcondb <- function(n, media, vari) {
    sd <- sqrt(vari)
    return(rnorm(n, media, sd))
  }
  
  # Conditional distribution of compound variance of genetic effects (sigb)
  dcondsigb <- function(b, deltav, n, nu, Sc) {
    z <- sum(b ^ 2 * deltav)
    return(1 / rgamma(1, (n + nu) / 2, (z + nu * Sc) / 2))
  } 
  
  # Conditional distribution of residual compound variance 
  dcondsigsq <- function(Aux, n, nu, Sce) {
    return(1 / rgamma(1, (n + nu) / 2, crossprod(Aux) / 2 + Sce / 2))
  }
  
  # Conditional fixed effects
  rmvnor <- function(n,media,sigma){
    z <- rnorm(n)
    return( media + crossprod(chol(sigma), z))
  }
  
  # Function for eigen descompositions
  eig <- function(K, tol) {
    ei <- eigen(K)
    fil <- which(ei$values > tol)
    return(list(ei$values[fil], ei$vectors[, fil]))
  }
  
  # Set spectral decomposition
  setDEC <- function(K, tol, ne) {
    
    sK <- vector("list", length = length(K))
    typeM <- sapply(K, function(x) x$Type)
    
    if (!all(typeM %in% c("BD", "D")))
      stop("Matrix should be of types BD or D")
    
    if (missing(ne)) {
      if (any(typeM == "BD"))
        stop("For type BD, number of subjects in each sub matrices should be provided")
    } else {
      if (length(ne) <= 1 & any(typeM == "BD"))
        stop("ne invalid. For type BD, number of subjects in each sub matrices should be provided")
      
      nsubK <- length(ne)
      
      if (nsubK > 1) {
        posf <- cumsum(ne)
        posi <- cumsum(c(1,ne[-length(ne)]))
      }
    }
    
    
    for (i in 1:length(K)) {
      if (K[[i]]$Type == "D") {
        tmp <- list()
        ei <- eig(K[[i]]$Kernel, tol)
        tmp$s <- ei[[1]]
        tmp$U <- ei[[2]]
        tmp$tU <- t(ei[[2]])
        tmp$nr <- length(ei[[1]])
        tmp$type <- "D"
        tmp$pos <- NA
        sK[[i]] <- list(tmp)
      }
      
      if (K[[i]]$Type == "BD") {
        cont <- 0
        tmp <- list()
        for (j in 1:nsubK) {
          Ktemp <- K[[i]]$Kernel[(posi[j]:posf[j]), (posi[j]:posf[j])]
          ei <- eig(Ktemp, tol)
          if (length(ei[[1]]) != 0) {
            cont <- cont + 1
            tmp[[cont]] <- list()
            tmp[[cont]]$s <- ei[[1]]
            tmp[[cont]]$U <- ei[[2]]
            tmp[[cont]]$tU <- t(ei[[2]])
            tmp[[cont]]$nr <- length(ei[[1]])
            tmp[[cont]]$type <- "BD"
            tmp[[cont]]$pos <- c(posi[j], posf[j])
          }
        }
        
        if (length(tmp) > 1) {
          sK[[i]] <- tmp
        } else {
          sK[[i]] <- list(tmp[[1]])
        }
      }
    }
    return(sK)
  }
  
  # verbose part I
  if (as.numeric(verbose) != 0) {
    cat("Setting parameters...", "\n", "\n")
  }
  
  ### PART II  Preparation for Gibbs sampler ######
  # Identification of NA's values
  y <- as.numeric(y)
  yNA <- is.na(y)
  whichNa <- which(yNA)
  nNa <- length(whichNa)
  mu <- mean(y, na.rm = TRUE)
  n <- length(y)
  
  if (nNa > 0) { 
    y[whichNa] <- mu
  }
  
  # name of each kernel (important to following procedures)
  if (is.null(names(K))) {
    names(K) <- paste("K", seq(length(K)), sep = "")
  }
  
  # initial values of fixed effects 
  if (!is.null(XF)) {
    Bet <- solve(crossprod(XF), crossprod(XF, y))
    nBet <- length(Bet)
    tXX <- solve(crossprod(XF))
  }
  
  # Eigen descomposition for nk Kernels
  nk <- length(K)
  nr <- numeric(length(K))
  typeM <- sapply(K, function(x) x$Type)
  
  if (!all(typeM %in% c("BD", "D")))
    stop("Matrix should be of types BD or D")
  
  if (length(ne) == 1 & any(typeM == "BD"))
    stop("Type BD should be used only for block diagonal matrices")
  
  Ei <- setDEC(K = K, ne = ne, tol = tol)
  
  # Initial values for Monte Carlo Markov Chains (MCMC)
  
  nCum <- sum(seq(1, ite) %% thin == 0)
  
  chain <- vector("list", length = 3)
  names(chain) <- c("mu", "varE", "K")
  chain$varE <- numeric(nCum)
  chain$mu <- numeric(nCum)
  
  chain$K <- vector("list", length = nk)
  names(chain$K) <- names(K)
  chain$K[seq(nk)] <- list(list(varU = numeric(nCum)))
  
  cpred <- vector('list', length = nk)
  names(cpred) <- names(K)
  cpred[seq(nk)] <- list(U = matrix(NA_integer_, nrow = nCum, ncol = n))
  
  
  nu <- 3
  Sce <- (nu + 2) * (1 - R2) * var(y, na.rm = TRUE)
  Sc <- numeric(length(K))
  for (i in 1:length(K)) {
    Sc[i] <- (nu + 2) * R2 * var(y, na.rm = T)/mean(diag(K[[i]]$Kernel))
  }
  tau <- 0.01
  u <- list()
  for (j in 1:nk) {
    u[[j]] <- rnorm(n, 0, 1 / (2 * n))
  }
  sigsq <- var(y)
  #Sc <- rep(0.01, nk)
  sigb <- rep(0.2, nk)
  
  temp <- y - mu
  
  if (!is.null(XF)) {
    B.mcmc <- matrix(0, nrow = nCum, ncol = nBet)
    temp <- temp - XF %*% Bet
  }
  
  temp <- temp - Reduce('+', u)
  nSel <- 0
  i <- 1
  
  ### PART III  Fitted model with training data ####
  
  # Iterations of Gibbs sampler
  while (i <= ite) {
    time.init <- proc.time()[3]
    
    # Conditional of mu
    temp <- temp + mu
    mu <- rnorm(1, mean(temp), sqrt(sigsq/n))
    #mu.mcmc[i] <- mu
    temp <- temp - mu
    
    # Conditional of fixed effects
    if (!is.null(XF)) {
      temp <- temp + XF %*% Bet
      vari <- sigsq * tXX
      media <- tXX %*% crossprod(XF, temp)
      Bet <- rmvnor(nBet, media, vari)
      temp <- temp - XF %*% Bet
    }
    
    # Conditionals x Kernel
    for (j in 1:nk) {
      # Sampling genetic effects
      
      if (typeM[j] == "D") {
        temp <- temp + u[[j]]
        d <- crossprod(Ei[[j]][[1]]$U, temp)
        s <- Ei[[j]][[1]]$s
        deltav <- 1 / s
        lambda <- sigb[j]
        vari <- s * lambda / (1 + s * lambda * tau)
        media <- tau * vari * d
        nr <- Ei[[j]][[1]]$nr
        b <- dcondb(nr, media, vari)
        u[[j]] <- crossprod(Ei[[j]][[1]]$tU ,b)
        temp <- temp - u[[j]]
      }
      
      if (typeM[j] == "BD") {
        nsk <- length(Ei[[j]])
        
        if (length(nsk > 1)) {
          temp <- temp + u[[j]]
          d <- NULL
          s <- NULL
          neiv <- numeric(nsk)
          pos <- matrix(NA, ncol = 2, nrow = nsk)
          for (k in 1:nsk) {
            pos[k,] <- Ei[[j]][[k]]$pos
            d <- c(d, crossprod(Ei[[j]][[k]]$U, temp[pos[k, 1]:pos[k, 2]]))
            neiv[k] <- length(Ei[[j]][[k]]$s)
            s <- c(s, Ei[[j]][[k]]$s)
          }
          deltav <- 1/s
          lambda <- sigb[j]
          vari <- s * lambda / (1 + s * lambda * tau)
          media <- tau*vari*d
          nr <- length(s)
          b <- dcondb(nr, media, vari)
          
          posf <- cumsum(neiv)
          posi <- cumsum(c(1, neiv[-length(neiv)]))
          utmp <- numeric(n)
          for (k in 1:nsk) {
            utmp[pos[k, 1]:pos[k, 2] ] <- crossprod(Ei[[j]][[k]]$tU, b[posi[k]:posf[k] ])
          }
          u[[j]] <- utmp
          temp <- temp - u[[j]]
          
        }else{
          temp <- temp + u[[j]]
          pos <- Ei[[j]]$pos
          d <- crossprod(Ei[[j]][[1]]$U, temp[pos[1]:pos[2]])
          s <- Ei[[j]][[1]]$s
          deltav <- 1/s
          lambda <- sigb[j]
          vari <- s * lambda / (1 + s * lambda * tau)
          media <- tau*vari*d
          nr <- Ei[[j]][[1]]$nr
          b <- dcondb(nr, media, vari)
          utmp <- numeric(n)
          utmp[pos[1]:pos[2]] <-  crossprod(Ei[[j]][[1]]$tU, b)
          u[[j]] <- utmp
          temp <- temp - u[[j]]
        }
      }
      
      # Sampling scale hyperparameters and variance of genetic effects
      sigb[j] <- dcondsigb(b, deltav, nr, nu, Sc[j])
    }
    
    # Sampling residual variance 
    res <- temp
    #Sce <- dcondSc(nu, sigsq)
    sigsq <- dcondsigsq(res, n, nu, Sce)
    tau <- 1 / sigsq
    
    # Predicting missing values
    if (nNa > 0) {
      uhat <- Reduce('+', u)
      
      if (!is.null(XF)) {
        aux <- XF[yNA,] %*% Bet
      }else{
        aux <- 0
      }
      
      y[whichNa] <- aux + mu + uhat[whichNa] + rnorm(n = nNa, sd = sqrt(sigsq))
      temp[whichNa] <- y[whichNa] - uhat[whichNa] - aux - mu
    }
    
    # Separating what is for the chain
    if (i %% thin == 0) {
      nSel <- nSel + 1
      chain$varE[nSel] <- sigsq
      chain$mu[nSel] <- mu
      if (!is.null(XF)) {
        B.mcmc[nSel,] <- Bet
      }
      for (j in 1:nk) {
        cpred[[j]][nSel,] <- u[[j]]
        chain$K[[j]]$varU[nSel] <- sigb[j]
      }
    }
    
    # Verbose 
    if (as.numeric(verbose) != 0 & i %% as.numeric(verbose) == 0) {
      time.end <- proc.time()[3]
      cat("Iter: ", i, "time: ", round(time.end - time.init, 3),"\n")
    }
    i <- i + 1
  }
  
  
  ###### PART IV  Output ######
  #Sampling
  draw <- seq(ite)[seq(ite) %% thin == 0] > burn
  
  mu.est <- mean(chain$mu[draw])
  yHat <- mu.est
  
  if (!is.null(XF)) {
    B <- colMeans(B.mcmc[draw,])
    yHat <- yHat + XF %*% B
  }
  
  u.est <- sapply(cpred, FUN = function(x) colMeans(x[draw, ]) )
  yHat <- yHat + rowSums(u.est)
  
  out <- list()
  out$yHat <- yHat
  out$varE <- mean(chain$varE[draw])
  out$varE.sd <- sd(chain$varE[draw])
  
  out$K <- vector('list', length = nk)
  names(out$K) <- names(cpred)
  for (i in 1:nk) {
    out$K[[i]]$u <- colMeans(cpred[[i]][draw, ])
    out$K[[i]]$u.sd <- apply(cpred[[i]][draw, ], MARGIN = 2, sd)
    out$K[[i]]$varu <- mean(chain$K[[i]]$varU[draw])
    out$K[[i]]$varu.sd <- sd(chain$K[[i]]$varU[draw])
  }
  
  out$chain <- chain
  out$ite <- ite
  out$burn <- burn
  out$thin <- thin
  out$model <- K$model
  out$kernel <- K$kernel
  out$y <- y
  class(out) <- "BGGE"
  return(out)
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
  if(size_E == 'environment') for(q in 1:length(K_E)) K_E[[q]] <- EnvKernel(env.data = K_E[[q]],Y = Y,merge = TRUE,env.id = 'env')$envCov
  
  
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

#install_github('cran/dofuture')
require(doFuture)

source('https://raw.githubusercontent.com/gcostaneto/KernelMethods/master/GBLUP_Kernel.R')
#require(EnvRtype)

### genomic kernels #####
pheno <- droplevels(pheno[pheno$gid %in% rownames(M),])
GB <- GB_kernel(M = M)

dim(GB$Ga)

K_W1 <- EnvKernel(env.data = W1,Y = pheno,bydiag = T,merge = T)[[2]]
K_W2 <- EnvKernel(env.data = W2,Y = pheno,bydiag = T,merge = T)[[2]]

# Modelo I: basico y = b + A
MM  <- get_kernel(K_G = list(A=GB$Ga),model = 'MM',Y = pheno)
# Modelo I: basico y = b + A + AE
MDs <- get_kernel(K_G = list(A=GB$Ga),model = 'MDs',Y = pheno)

E1MM <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W1),model = 'EMM',Y = pheno)
E2MM <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W2),model = 'EMM',Y = pheno)
E3MM <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W3),model = 'EMM',Y = pheno)
# <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W4),model = 'EMM',Y = pheno)

E1MD <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W1),model = 'EMDs',Y = pheno)
E2MD <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W2),model = 'EMDs',Y = pheno)
E3MD <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W3),model = 'EMDs',Y = pheno)


E1RMD <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W1),model = 'RNMDs',Y = pheno)
E2RMD <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W2),model = 'RNMDs',Y = pheno)
E3RMD <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W3),model = 'RNMDs',Y = pheno)
#E4MD <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W4),model = 'EMDs',Y = pheno)

E1RMM <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W1),model = 'RNMM',Y = pheno)
E2RMM <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W2),model = 'RNMM',Y = pheno)
E3RMM <-get_kernel(K_G = list(A=GB$Ga),K_E = list(W=K_W3),model = 'RNMM',Y = pheno)


Kernels <- list(MM,MDs,E1MM,E2MM,E3MM,
                E1MD,E2MD,E3MD,
                E1RMD,E2RMD,E3RMD,
                E1RMM,E2RMM,E3RMM)

models <-c('MM','MDs','EMM1','EMM2','EMM3',
           'EMD1','EMD2','EMD3',
           'E1RMD','E2RMD','E3RMD',
           'E1RMM','E2RMM','E3RMM')

traits <- c('EH','PH','GY')

#Kernels <- list(M1usp,M1hel,M5usp,M5hel)


rm(M1usp,M1hel,M2usp,M2hel, M3usp,M3hel,M4usp,M4hel)

cat(paste0('Kernels: done  ',Sys.time(),'\n'))
#### CV SCHEMES ###########
source('https://raw.githubusercontent.com/gcostaneto/SelectivePhenotyping/master/cvrandom.R')
TS <- Sampling.CV0(gids = pheno$gid,ngids = nlevels(pheno$gid),envs = pheno$env,f = .7,out.env = 1,seed =seed,rep=nFold)


require(doParallel)

cl <- makeCluster(13)
registerDoParallel(cl)


results <-  foreach(t = 1:length(traits), .combine = "rbind")%:%
  foreach(f = 1:nFold, .combine = "rbind")%:%
  foreach(m = 1:length(models), .combine = "rbind")%dopar% {
    
    Y <- phenoGE[[t]]
    tr <- TS[[f]][[4]]
    Y$value[-tr] <- NA
    # Ze <- model.matrix(~0+env,Y)
    ne <- as.vector(table(Y$env))
    fit <- BGGE(y = Y$value,
                K = Kernels[[m]],
                # XF = Ze,
                tol = 1e-20,
                ne = ne,
                ite = ite,
                burn = burn,
                thin = thin,
                verbose = FALSE)
    
    output <- data.frame(obs=phenoGE[[t]]$value,pred=fit$yHat,
                         gid=phenoGE[[t]]$gid, env=phenoGE[[t]]$env,
                         trait = traits[t],
                         Model= models[m],rep=f,pop=NA)
    
    output$pop[TS[[f]][[1]]] <- 'NewG'
    output$pop[TS[[f]][[2]]] <- 'NewE'
    output$pop[TS[[f]][[3]]] <- 'NewGE'
    output$pop[TS[[f]][[4]]] <- 'training'
    
    Vcomp = data.frame(Vcomp.BGGE(model = fit,5),rep=f,Model=as.character(models[m]))
    write.table(x=Vcomp,file='Vcomps0.txt',sep=',',append = T,row.names=F)
    
    df<-data.frame(Trait = 'GY', Model = models[m],rep=f,
                   rTr=cor(phenoGE[[t]]$value[tr ], fit$yHat[tr ],use = 'complete.obs'),
                   rTs=cor(phenoGE[[t]]$value[-tr], fit$yHat[-tr],use = 'complete.obs'))
    
    write.table(x = df,file = 'PA.txt',sep=',',append = T,row.names=F)
    return(output)
  }
stopCluster(cl)

#save(   object = results, file = 'results.Rdata')
saveRDS(object = results, file = 'maizeResults' )
Sys.time()





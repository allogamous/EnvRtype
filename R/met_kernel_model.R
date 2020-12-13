#'@title  Kernel Models for Predicting Phenotypes across Multi-Environment Conditions
#'
#'
#' @description Runs Bayesian Linear-Mixed Models for Multiple Environments using kernels from get_kernel function
#'
#' @author Germano Costa Neto (Adapted from Granato et al. 2018, BGGE package)

#' @param phenotypes Vector of data. Should be numeric and NAs are allowed.
#' @param gid Vector of characters for genotypes identification.
#' @param env Vector of characters for environments identification.
#' @param random list A two-level list Specify the regression kernels (co-variance matrix). The former is the \code{Kernel},
#' where is included the regression kernels. The later is the \code{Type}, specifying if the matrix is either \code{D} Dense or
#' \code{BD} Block Diagonal. A number of regression kernels or random effects to be fitted are specified in this list.
#' @param fixed matrix Design matrix (\eqn{n \times p}) for fixed effects
#' @param iterations numeric Number of iterations.
#' @param burnin numeric Number of iterations to be discarded as burn-in.
#' @param thining numeric Thinin interval.
#' @param digits numeric. Digits for round variance components.
#' @param verbose Should iteration history be printed on console? If TRUE or 1 then it is printed,
#' otherwise, if another number $n$ is choosen the history is printed every $n$ times. The default is \code{FALSE}
#' @param tol a numeric tolerance level. Eigenvalues lower than \code{tol} are discarded. Default is \code{1e-10}.
#' @param R2 the proportion of variance expected to be explained by the regression.
#'
#'
#' @details
#' TODO
#'

#' @return
#'  A list with estimated posterior means of residual and genetic variance component for each term in the linear model and the genetic value predicted. Also the
#'  values along with the chains are released.


kernel_model <- function(phenotypes, random = NULL, fixed = NULL,env, gid, verbose=FALSE, iterations=1E3, burnin=2E2, thining=10, tol=1e-20, R2=0.5, digits=4 ){

  cat(paste0('--------------------------------------------------------','\n'))
  cat(paste0('Running BGGE (Bayesian Genotype + Genotype x Environment)','\n'))
  cat(paste0('More Detail in Granato et al (2018) G3','\n'))
  cat(paste0('--------------------------------------------------------','\n'))

  if(is.null(random)) stop('Missing the list of kernels for random effects')
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


  Vcomp.BGGE<-function(model,env,gid,digits=digits,alfa=.10){

    t = length(unique(gid))
    e = length(unique(env))
    #GLres = p*q - (p-1) - (q-1)
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


    comps$Type <- comps$K


    comps$Type[grep(comps$K,pattern = 'KGE_')] = 'GxE'
    comps$Type[grep(comps$K,pattern = 'KG_')] = 'Genotype (G)'
    comps$Type[grep(comps$K,pattern = 'E')] = 'GxE'
    comps$Type[grep(comps$K,pattern = 'KE_')] = 'Environment (E)'

    comps$CI_upper = NA
    comps$CI_lower = NA

      ENV = which(comps$Type %in% 'Environment (E)')
      GID = which(comps$Type %in% 'Genotype (G)')
      GE = which(comps$Type %in% 'GxE')
      R = which(comps$Type %in% 'Residual')
    comps$CI_upper[ENV] = (n-e)  *comps$Var[ENV]/qchisq((alfa/2), n-e)
    comps$CI_upper[GID] = (n-t)  *comps$Var[GID]/qchisq((alfa/2), n-t)
    comps$CI_upper[GE ] = (n-t-e)*comps$Var[GE]/qchisq((alfa/2), n-t-e)
    comps$CI_upper[R  ] = (n-t-e)*comps$Var[R]/qchisq((alfa/2), n-t-e)
    comps$CI_lower[ENV] = (n-e)  *comps$Var[ENV]/qchisq((1-alfa/2), n-e)
    comps$CI_lower[GID] = (n-t)  *comps$Var[GID]/qchisq((1-alfa/2), n-t)
    comps$CI_lower[GE ] = (n-t-e)*comps$Var[GE]/qchisq((1-alfa/2), n-t-e)
    comps$CI_lower[R  ] = (n-t-e)*comps$Var[R]/qchisq((1-alfa/2), n-t-e)

    comps$CI_upper = round(comps$CI_upper,digits)
    comps$CI_lower = round(comps$CI_lower,digits)

    comps <- comps[,c(4,1:2,6,5,3)]


    return(comps)
  }


  ne <- as.vector(table(env))
  start=Sys.time()
  fit <- BGGE(y = phenotypes,
              K = random,
              tol = tol,
              ne = ne,
              ite = iterations,
              burn = burnin,
              thin = thining,
              verbose = verbose)
  end = Sys.time()

  cat(paste0('Start at: ', start,' | Ended at: ', end,'\n'))

  return(list(fitted = fit, VarComps = Vcomp.BGGE(model = fit,env = env,gid = gid, digits=digits)))

}



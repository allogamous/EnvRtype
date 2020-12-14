# **Genomic Prediction using Environmental Covariates and Considering Reaction-Norms**
* author: Germano Costa Neto

*last update: 13th December 2020*

* [Software](#P1)
* [Data Sets](#P2)
* [Environmental Covariables (ECs)](#P3)
* [Environmental Relatedness Kernels](#P4)
* [Preparing the Kernels for Prediction](#P5)
* [Fitting bayesian kernel models ](#P6)
* [Cross-validation to assess predictive ability](#P7)
* [Author's Suggestions](#P8)

<div id="P1" />


# Software

```{r, eval=FALSE}
library(devtools)
install_github('allogamous/EnvRtype')
library(EnvRtype)
```

<div id="P2" />

### Data sets

> * Toy example using 
```{r, eval=FALSE}
library(EnvRtype)
data("maizeG")
data("maizeWTH")
data("maizeYield")
Y  = maizeYield
```

<div id="P3" />

### Environmental Covariables (ECs) for **W** Matrix (W.matrix function)

```{r, eval=FALSE}
## Organizing Environmental Covariables (ECs) in W matrix

> * Data were organized for different development stages in maize. We assume fixed time intervals from the days after planting.

stages    = c('VE','V1_V6','V6_VT','VT_R1','R1_R3','R3_R6',"H")
interval = c(0,7,30,65,70,84,105)
id.vars  = names(maizeWTH)[c(10:15,23,25:30)]

W.matrix = W_matrix(env.data = maizeWTH,env.id = 'env',
                    var.id = id.vars,by.interval = T,time.window = interval,
                    names.window = stages,center = F,scale = F )

```

<div id="P4" />

### Environmental Relatedness Kernels (env_kernel function)

```{r, eval=FALSE}
## Kernel for the involving all development stages
K_F <- env_kernel(env.data = W.matrix,gaussian = T)[[2]]

## Kernels for each development stage
K_S <- env_kernel(env.data = W.matrix,gaussian = T,stages = stages[2:5])[[2]]

# K_G (genotype) and K_E (environment) must be a list of kernels
# So:
K_G = list(G = maizeG)
# And:
K_F <- list(E = K_F)

# for K_S, we dont need to create a list because K_S is already a list of kernels for each development stage
```

<div id="P5" />

### Preparing the Kernels for Prediction (get_kernel function)

> * the get_kernel function creates the modeling structure for predictive purposes

> * In this example, we show the use of the Reaction-Norm Main Effect Model (RNMM), assuming: y = intercept + fixed effects + enviromic + genomic + enviromic x genomic + error.

```{r, eval=FALSE}
## Assembly Genomic and Enviromic Kernel Models
M1 = get_kernel(K_G = K_G, Y = Y, model = "MDs") # baseline model
M2 = get_kernel(K_G = K_G, K_E = K_F, Y = Y, model = "RNMM",dimension_KE = 'q') # reaction-norm 1
M3 = get_kernel(K_G = K_G, K_E = K_S, Y = Y, model = "RNMM",dimension_KE = 'q') # reaction-norm 2
```

<div id="P6" />


### Fitting bayesian kernel models (kernel_model function)

> * kernel_model runs a bayesian kernel model regression.
> * in the example below, we run the function in order to compute the variance components for each genomic and enviromic effect

```{r, eval=FALSE}
model = c('Baseline Genomic','Reaction-Norm','Reaction-Norm for each Dev.Stage')
Models = list(M1,M2,M3) # for running all models in loop

iter = 10E3 # number of iterations
burn = 5E3  # number of burn in
thin = 10   # number for thining

Z_E = model.matrix(~0+env,data=Y) # fixed environmental effects

Vcomp <- c()
for(MODEL in 1:length(Models)){
  fit <- kernel_model(phenotypes = Y$value,env = Y$env,gid = Y$gid,
                      random = Models[[MODEL]],fixed = Z_E,
                      iterations = iter,burnin = burn,thining = thin)
  
  Vcomp <- rbind(Vcomp,data.frame(fit$VarComps,Model=model[MODEL]))
}

```

<div id="P7" />

### Cross-validation to assess predictive ability of GP models (kernel_model function)

> * creating the training sets (prediction scenario: prediction of novel genotypes, CV1)
> * to speed up your analysis, we suggest to use foreach to run the cross-validation

```{r, eval=FALSE}
source('https://raw.githubusercontent.com/gcostaneto/SelectivePhenotyping/master/cvrandom.R')
gid  = Y$gid
env  = Y$env
rep  = 10
seed = 7121
f    = 0.80
iter = 5E3
burn = 1E3
thin = 10
TS = Sampling.CV1(gids = Y$gid,f = f,seed = seed,rep = rep,gidlevel = F)
```

> * running predictive models

```{r, eval=FALSE}
require(foreach)
results <-foreach(REP = 1:rep, .combine = "rbind")%:%
  foreach(MODEL = 1:length(model), .combine = "rbind")%dopar% {
    
    
    yNA      <- Y$value
    tr       <- TS[[REP]]
    yNA[-tr] <- NA
    
    Z_E = model.matrix(~0+env,data=Y) # fixed environmental effects
    
    fit <- kernel_model(phenotypes = yNA,env = Y$env,gid = Y$gid,
                        random = Models[[MODEL]],fixed = Z_E,
                        iterations = iter,burnin = burn,thining = thin)
   

    df<-data.frame(Model = model[MODEL],rep=REP,
                   rTr=cor(Y$value[tr ], fit$fitted$yHat[tr ],use = 'complete.obs'),
                   rTs=cor(Y$value[-tr], fit$fitted$yHat[-tr],use = 'complete.obs'))
    
    write.table(x = df,file = 'PA_models.txt',sep=',',append = T,row.names=T)
    return(df)
  }
  
saveRDS(object = results, file = 'PA_cv1_2' )
require(plyr)

# predictive ability
ddply(results,.(Model),summarise, pa = round(mean(rTs),3),sd = round(sd(rTs),4))

```

<div id="P8" />

### Our Suggestions

> * Although we provide a function to run predictions, the user can also use the kernels created in get_kernel to run analyzes on other packages, such as BGLR or asreml. These packages are interesting because they also allow the modeling of the variance-covariance structure of the residues (errors), in addition to allows the inclusion of factor analytic structures, which can increase the predictive ability of the models;

> * When running genomic prediction (GP) for multiple environments, we enviasage that CV schemes such as CV0 (prediction of novel environments) and CV00 (prediction of novel genotypes at novel environments) are indispensable to better assess the potential predictive ability of yours GP models. The power of models with ECs will be better highlighted when prediction scenarios involve a lack of environmental inforamtion. Schemes such as CV1 or CV2, perhpas, will not be enought to discriminate the predictive ability of the  non-enviromic and enviromic-based GP models.


> * Finally, in addition to the predictive ability assessed at model level, for multi-environment conditons is necessary to assess the predictive ability of the models at each environment. Thus, a better discrimination of the model's potentialities can be achvied.


### CV1 and CV00

> * CV1
```{r, eval=FALSE}
gid  = Y$gid
env  = Y$env
rep  = 30
seed = 1010
f    = 0.20
iter = 2E3
burn = 5E2
thin = 10

source('https://raw.githubusercontent.com/gcostaneto/SelectivePhenotyping/master/cvrandom.R')
TS = Sampling.CV1(gids = Y$gid,f = f,seed = seed,rep = rep,gidlevel = F)


require(foreach)
require(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)


results <-foreach(REP = 1:rep, .combine = "rbind")%:%
  foreach(MODEL = 1:length(model), .combine = "rbind")%dopar% {
    
    
    yNA      <- Y$value
    tr       <- TS[[REP]]
    yNA[-tr] <- NA
    
    Z_E = model.matrix(~0+env,data=Y)
    fit <- kernel_model(phenotypes = yNA,env = Y$env,gid = Y$gid,
                        random = Models[[MODEL]],fixed = Z_E,
                        iterations = iter,burnin = burn,thining = thin)
    
    
    df<-data.frame(Model = model[MODEL],rep=REP,
                   rTr=cor(Y$value[tr ], fit$fitted$yHat[tr ],use = 'complete.obs'),
                   rTs=cor(Y$value[-tr], fit$fitted$yHat[-tr],use = 'complete.obs'))
    
    cat(paste0(model[MODEL],' ',REP," r = ",round(cor(Y$value[-tr], fit$fitted$yHat[-tr]),3),'\n'))
    
    write.table(x = df,file = 'PA_models.txt',sep=',',append = T,row.names=T)
    
    output <- data.frame(obs=Y$value,pred=fit$fitted$yHat,
                         gid=Y$gid, env=Y$env,
                         Model = model[MODEL],rep=REP,pop=NA)
    
    
    output$pop[tr ] <- 'training'
    output$pop[-tr] <- 'testing'
    
    return(output)
  }

stopCluster(cl)

pa = ddply(results,.(rep,pop,Model),summarise,r = cor(obs,pred))
ddply(pa,.(pop,Model),summarise, pa = round(mean(r),3),sd = round(sd(r),3))
```

> * CV00 

```{r, eval=FALSE}
seed = 8172
TS=Sampling.CV0(gids = gid,envs = env,out.env = 2,f = f,seed = seed,rep = rep)

cl <- makeCluster(3)
registerDoParallel(cl)

results <-foreach(REP = 1:rep, .combine = "rbind")%:%
  foreach(MODEL = 1:length(model), .combine = "rbind")%dopar% {
    
    
    yNA      <- Y$value
    tr       <- TS[[REP]]$training
    yNA[-tr] <- NA
    
    Z_E = model.matrix(~0+env,data=Y)
    fit <- kernel_model(phenotypes = yNA,env = Y$env,gid = Y$gid,
                        random = Models[[MODEL]],fixed = Z_E,
                        iterations = iter,burnin = burn,thining = thin)
    
    
    output <- data.frame(obs=Y$value,pred=fit$fitted$yHat,
                         gid=Y$gid, env=Y$env,
                         Model = model[MODEL],rep=REP,pop=NA)
    
    
    output$pop[tr ] <- 'training'
    output$pop[-tr] <- 'testing'
    
    return(output)
  }
  
  stopCluster(cl)

pa = ddply(results,.(rep,pop,Model),summarise,r = cor(obs,pred))
ddply(pa,.(pop,Model),summarise, pa = round(mean(r),3),sd = round(sd(r),3))
```

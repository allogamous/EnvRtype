# **Genomic Prediction using Environmental Covariates and Considering Reaction-Norms**

*last update: 13th December 2020*

* [Software](#P1)
* [Data Sets](#P2)
* [Soil Features](#P3)
* [Environmental Typing](#P4)
* [Environmental Covaraibles](#P5)
* [Environmental Similarity](#P6)
              
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

### Preparing the Kernels for Prediction (get_kernel function)

> * In this example, we show the use of the Reaction-Norm Main Effect Model, assuming: $y = 1 \mu + X \Beta$
```{r, eval=FALSE}
## Assembly Genomic and Enviromic Kernel Models
M1 = get_kernel(K_G = K_G, Y = Y, model = "MDs") # baseline model
M2 = get_kernel(K_G = K_G, K_E = K_F, Y = Y, model = "RNMM",dimension_KE = 'q') # reaction-norm 1
M3 = get_kernel(K_G = K_G, K_E = K_S, Y = Y, model = "RNMM",reaction = T,dimension_KE = 'q') # reaction-norm 2
```

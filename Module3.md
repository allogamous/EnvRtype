# Reaction Norm Module
                        
> * We provide Genomic and Envirotypic kernels for reaction norm prediction. After generate the kernels, the user must use the [BGGE](https://github.com/italo-granato/BGGE) package to run the models
                      
### **Toy Example: genomic prediction for grain yield in tropical maize**

```{r}
require(EnvRtype)
data("maizeYield") # 150 maize hybrids over 5 environments (grain yield data)
data("maizeG")     # GRM for maizeYield
data('maizeWTH')   # weather data for maize Yield

Y <- maizeYield
G <- maizeG
df.clim <- maizeWTH
```

### **Statistical Models**
                      
<p align="center">
<img src="/fig/summary_pkg.png" width="70%" height="70%">
</p>
                        
> * Returns benchmark main effect model: 
                        
<p align="center">
<img width="120" height="18" src="/fig/mod1.png">
 </p>
                        
```{r}
MM <- get_kernel(K_G = list(G=G),Y = Y,reaction = F,model = 'MM')
```
> *  Returns benchmark main GxE deviation model:
                        
<p align="center">
<img width="160" height="18" src="/fig/mod2.png">
</p>

```{r}
MDs <-get_kernel(K_G = list(G=G),Y = Y,reaction = F,model = 'MDs')
```
> * Obtaining environmental variables based on quantiles
                      
```{r}
Env.data<-W.matrix(env.data = df.clim,by.interval = T,statistic = 'quantile',time.window = c(0,14,35,60,90,120))
dim(Env.data)
```
> * Creating Env Kernels from W matrix and Y dataset

```{r}
H <- EnvKernel(env.data = Env.data,Y = Y,merge = T,env.id = 'env')
dim(H)
dim(H$varCov) # variable relationship
dim(H$envCov) # environmental relationship

require(superheat)
superheat(H$envCov,row.dendrogram = T,col.dendrogram = T)
 ```
 
> * Parametrization by 
                      
<p align="center">
<img width="110" height="50" src="/fig/mod3.png">
</p>
                        
```{r}
H <- EnvKernel(env.data = Env.data,Y = Y,merge = T,env.id = 'env',bydiag=FALSE)
dim(H)
dim(H$varCov) # variable relationship
dim(H$envCov) # environmental relationship
superheat(H$envCov,row.dendrogram = T,col.dendrogram = T)
```

> * Parametrization by 
                      
<p align="center">
<img width="130" height="50" src="/fig/mod4.png">
</p>
resulting in diag(K_W) = 1
                      
```{r}
H <- EnvKernel(env.data = Env.data,Y = Y,merge = T,env.id = 'env',bydiag=TRUE)
dim(H)
dim(H$varCov) # variable relationship
dim(H$envCov) # environmental relationship
superheat(H$envCov,row.dendrogram = T,col.dendrogram = T)
 ```
 
> * Gaussian parametrization by 
                      
<p align="center">
<img width="130" height="50" src="/fig/mod5.png">
</p>

which d = dist(W), q = median(d) and h = gaussian parameter (default = 1)
                      
```{r}
H <- EnvKernel(env.data = Env.data,Y = Y,merge = T,env.id = 'env',gaussian=TRUE)
dim(H)
dim(H$varCov) # variable relationship
dim(H$envCov) # environmental relationship
superheat(H$envCov,row.dendrogram = T,col.dendrogram = T)
```
**________________________________________________________________________________________________________**  
                        
**Attention**:\
K_G = list of genomic kernels;\
K_E = list of environmental kernels;\
reaction = TRUE, build the haddamard's product between genomic and envirotype-based kernels;\
reaction = FALSE, but K_E != NULL, only random environmental effects using K_E are incorporated in the model  

**________________________________________________________________________________________________________**  

- Returns benchmark main effect model plus random environmental covariables:

<p align="center">
  <img width="140" height="18" src="/fig/mod6.png">
</p>

```{r}
EMM <-get_kernel(K_G = list(G=G),K_E = list(W=H$envCov), Y = Y,model = 'EMM') 
```
- Returns benchmark main GxE deviation model plus random environmental covariables: 

<p align="center">
  <img width="180" height="18" src="/fig/mod7.png">
</p>

```{r}
EMDs <-get_kernel(K_G = list(G=G),Y = Y,K_E = list(W=H$envCov),model = 'EMDs') # or model = MDs

```
- Returns reaction norm model: 

<p align="center">
  <img width="200" height="18" src="/fig/mod8.png">
</p>

```{r}
RN <-get_kernel(K_G = list(G=G),K_E = list(W=H$envCov), Y = Y,model = 'RNMM')

```

- Returns a full reaction norm model with GE and GW kernels: 

<p align="center">
  <img width="220" height="18" src="/fig/mod9.png">
</p>

```{r}
fullRN <-get_kernel(K_G = list(G=G),K_E = list(W=H$envCov), Y = Y,model = 'RNMDs')

```

- Advanced options: **lets build again the W matrix**

```{r}
Env.data<-W.matrix(env.data = df.clim,by.interval = T,statistic = 'quantile',
                time.window = c(0,14,35,60,90,120))

W <- EnvKernel(env.data = Env.data,Y = Y,merge = T,env.id = 'env',bydiag=TRUE)

# by using size_E = 'environment', get_kernel directly takes a W of q x q environments and builds a n x n matrix as EnvKernel()
EMM <-get_kernel(K_G = list(G=G),K_E = list(W=W$envCov), Y = Y,,model = 'EMM',size_E = 'environment')


# Its possible to integrate more than one environmental kernel
T.cov<- EnvTyping(env.data=df.clim,var.id =  c('T2M','PRECTOT','WS2M'),env.id='env',format = 'wide')
eT <- EnvKernel(env.data =T.cov,Y = Y,merge = T,env.id = 'env',bydiag=TRUE)


EMM <-get_kernel(K_G = list(G=G),K_E = list(W=W$envCov,eT=eT$envCov), Y = Y,model = 'EMM',size_E = 'environment')
EMM$KE_W # kernel from W
EMM$KE_eT # kernel from T (envirotype)

```

- Integration with **BGGE package**

```{r}
require(BGGE)

 ne <- as.vector(table(maizeYield$env))
      fit <- BGGE(y = maizeYield$value,
                  K = EMM,
                  ne = ne,
                  ite = 1000,
                  burn = 100,
                  thin = 2,
                  verbose = TRUE)
```


# *EnvRtype*: a tool for envirotyping analysis and genomic prediction considering reaction norms  <img align="right" src="/fig/pkg_i.png" width="20%" height="20%">



<div id="menu" />
  
  ---------------------------------------------
  
  ## Background
  
  > Environmental typing (envirotyping) has proven useful in identifying the non-genetic drivers of phenotypic adaptation in plant breeding. Combined with phenotyping and genotyping data, the use of envirotyping data may leverage the molecular breeding strategies to cope with environmental changing scenarios. Over the last ten years, this data has been incorporated in genomic-enabled prediction models aiming to better model genotype x environment interaction (GE) as a function of reaction norm. However, there is difficult for most breeders to deal with the interplay between envirotyping, ecophysiology, and genetics. 
> Here we present the EnvRtype R package as a new toolkit developed to facilitate the interplay between envirotyping and genomic prediction. This package offers three modules: (1) collection and processing data set, (2) environmental characterization, (3) build of ecophysiological enriched genomic prediction models accounting for three different structures of reaction norm. Here we focus our efforts to present a practical use of EnvRtype package in supporting the genome-wide prediction of reaction norms. We provide a intuitive framework to integrate different reaction norm models in Bayesian Genomic Genotype x Environment Interaction (BGGE) package.

<div id="menu" />
  
  ---------------------------------------------
  ## Resources
  
  > EnvRtype consists in three modules (sections 2-4), which collectively generate a simple workflow to collect, process and integrates envirotyping data into genomic prediction over multiple environments.

* [1. Install and Required Packages](#Instal)
* [2. Environmental Sensing Module](#P1)
* [3. Environmental Characterization Module](#P2)
* [4. Reaction Norm Module](#P3)
* [5.Authorship](#P4)
* [6. Acknowledgments](#P5)
* [7. Simplified R script (Tutorial)](#P6)
              
              
<div id="Instal" />
                
## Install

### Using devtools in R

```{r}
library(devtools)
install_github('allogamous/EnvRtype')
require(EnvRtype)
  ```
### Manually installing

> If the method above doesn't work, use the next lines by downloading the EnvRtype-master.zip file

```{r}
setwd("~/EnvRtype-master.zip") # ~ is the path from where you saved the file.zip
unzip("EnvRtype-master.zip") 
file.rename("EnvRtype-master", "EnvRtype") 
shell("R CMD build EnvRtype") # or system("R CMD build EnvRtype")
install.packages("EnvRtype_0.1.9.tar.gz", repos = NULL, type="source") # Make sure to use the current verision
```
 
 ### Required packages
 
> * **[EnvRtype](https://github.com/allogamous/EnvRtype)** 
> * **[raster](https://CRAN.R-project.org/package=raster)** 
> * **[nasapower](https://github.com/ropensci/nasapower)** 
> * **[BGGE](https://github.com/italo-granato/BGGE)**
> * **[foreach](https://github.com/cran/foreach)**
> * **[doParalell](https://github.com/cran/doparallel)**
                
```{r}
install.packages("foreach")
install.packages("doParallel")
install.packages("raster")
install.packages("nasapower")
install.packages("BGGE")
              
or
              
source("https://raw.githubusercontent.com/gcostaneto/Funcoes_naive/master/instpackage.R");
              
inst.package(c("BGGE",'foreach','doParalell','raster','nasapower'));

library(EnvRtype)
              
```
<!-- toc -->
[Menu](#menu)
                  
 <div id="P1" />
                    
 ## Environmental Sensing Module
 
 ### Geographic Information Dabases
 
 > The collection, organization and processing of environmental data is a step that requires equipment installed in the field. Such equipment can be expensive or difficult to access for some research groups in certain regions or countries. For this reason, we decided to insert a routine for collecting climatic data through the [NASA POWER base](https://power.larc.nasa.gov/), which can access information on a daily scale anywhere on the globe.
 
 > The [Raster Package](https://cran.r-project.org/web/packages/raster/raster.pdf) also offers a digital platform for downloading files in raster format of climatic data (from the [WorldClim database](https://www.worldclim.org/)) and [SRTM (elevation)](http://srtm.csi.cgiar.org/) using only geographical coordinates (Latitude and Longitude).
 
 > To facilitate the use by researchers, especially in the field of genetics and plant breeding, we have integrated these platforms in the functions below:
 
> * Preparing de informations (latitude, longitude, start day and end date)

```{r}
lat = c(-13.05,-12.32,-18.34,-18.90,-23.03)  # vector of latitude WGS84
lon = c(-56.05,-55.42,-46.31,-49.56,-51.02)  # vector of lontitude WGS84
env = c("NM","SO","PM","IP","SE")            # vector of environment/site ID
plant.date = c("2015-02-15","2015-02-13",    # vector of start period
                                 "2015-02-26","2015-03-01",
                                 "2015-02-19") 
harv.date =rep("2015-06-15",5)               # vector of end period
```
> * So we can use this information to collect weather data from NASAPOWER

```{r}
df.clim <- get_weather(env.id = env,lat = lat,lon = lon, start.day = plant.date,end.day = harv.date, asdataframe = F) # returns a list of dataframes by environments
                  
df.clim <- get_weather(env.id = env,lat = lat,lon = lon,start.day = plant.date,end.day = harv.date,country = 'BRA') # returns a dataframe with all environments by default

head(df.clim)

```
<!-- toc -->
[Menu](#menu)

 ### Additional variables (ecophysiological)
 
 
> * Basic processing of get_weather() 


 ```{r}
 df.clim <-processWTH(env.data = df.clim)
```
> * Basic summary statistics for environmental data

```{r}
summaryWTH(df.clim)
# or
summaryWTH(df.clim,env.id = 'env')
 ```

> * Summary a particular environmental variable

```{r}
summaryWTH(df.clim,env.id = 'env',var.id = 'T2M')
summaryWTH(df.clim,env.id = 'env',var.id = c('T2M','T2M_MAX')) # or more than one
```

> * Summary by time intervals. Dividing the development cycle into time intervals (e.g., phenology), whether phenological or fixed time intervals (e.g. 10-day intervals) helps to understand the temporal variation of environmental factors during the crop growth cycle.

```{r}
summaryWTH(df.clim,env.id = 'env',by.interval = T)
```
> * Summary by time intervals given by *time.window* argument.

 ```{r}
summaryWTH(df.clim,env.id = 'env',by.interval = T,time.window = c(0,14,35,60,90,120))
```

> * Summary by time intervals given by *time.window* and *names.window*.

```{r}
summaryWTH(df.clim,env.id = 'env',by.interval = T,time.window = c(0,14,35,60,90,120), names.window = c('P-E','E-V1','V1-V4','V4-VT','VT-GF','GF-PM'))
 ```
 
> * Returns only mean values

```{r}
summaryWTH(df.clim,env.id = 'env',statistic = 'mean')
```

> * Returns only sum values

```{r}
summaryWTH(df.clim,env.id = 'env',statistic = 'sum')
```

> * Returns quantile values (default = 25%, 50% and 75%)

```{r}
summaryWTH(df.clim,env.id = 'env',statistic = 'quantile')
```

> * For specific quantiles (e.g., 20%, 76% and 90%)

```{r}
summaryWTH(df.clim,env.id = 'env',statistic = 'quantile',probs = c(.20,.76,.90))
```

<!-- toc -->
[Menu](#menu)


### Building Environmental Covariable Matrices
                  
> * Environmental variables can be used as indicators of the quality of an environment (experiment, location). A double entry table (*q* environments x *k* environmental factors) can be built as suggested by Jarquin et al (2014). Hereafther we will refer to this matrix as **W**, and therefore, it will be obtained by the function *W.matrix*:

```{r}
W.matrix(env.data = df.clim,by.interval = F)
```

> * As shown in the *summaryWTH* function, we can create time windows to capture the temporal variability between environmental information.To do this, we use the arguments *by.interval = TRUE* and *time.window* to define time limits. Such limits refer to the beginning of climate information. In the example below, for example, the first interval is for 0 to 14 days after planting, the second for 15 to 35 days after planting, and so on respectively.

```{r}
W.matrix(env.data = df.clim,by.interval = T, time.window = c(0,14,35,60,90,120))
```

> * Different statistics can be used, as in *summaryWTH*. The statistic argument is used to select between *mean* or *quantile*. If the selection is made in quantile, the complementary argument *prob* is used to choose the quantiles to be used.

```{r}
W.matrix(env.data = df.clim,by.interval = T,statistic = 'mean'    ,time.window = c(0,14,35,60,90,120))
W.matrix(env.data = df.clim,by.interval = T,statistic = 'quantile',time.window = c(0,14,35,60,90,120))
```

> *  We can perform a Quality Control (QC) based on the maximum sd tolered.

```{r}
W.matrix(env.data = df.clim,by.interval = F,QC = T)
```
> *  We can perform a Quality Control (QC) based on the maximum sd tolered

```{r}
W.matrix(env.data = df.clim,by.interval = F,QC = T,sd.tol = 3)
W.matrix(env.data = df.clim,by.interval = F,QC = T,sd.tol = 2)
```

> *  Create for specific variables. To do this, insert the name of the variables in the *id.var* argument.

```{r}
id.var = c('T2M_MAX','T2M_MIN','T2M') # maximum temperature, minimum temperature and average temperature
W.matrix(env.data = df.clim,var.id = id.var)
```

> *  Or even combine with summaryWTH by using *is.processed=TRUE*:

```{r}
data<-summaryWTH(df.clim,env.id = 'env',statistic = 'quantile')
W.matrix(env.data = data,is.processed = T)
```

[Menu](#menu)

<div id="P2" />
                      
## Environmental Characterization Module
                      
                      ### Environmental Typologies based on Cardinal Limits
                      
                      ```{r}
                    EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M')
                    ```
                    - Typologies by.intervals (generic time intervals)
                    ```{r}
                    EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T)
                    ```
                    - Typologies by.intervals (specific time intervals)
                    ```{r}
                    EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T,time.window = c(0,15,35,65,90,120))
                    ```
                    - Typologies by.intervals (specific time intervals and with specific names)
                    ```{r}
                    names.window = c('1-intial growing','2-leaf expansion I','3-leaf expansion II',
                                     '4-flowering','5-grain filling','6-maturation')
                    out<-EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T,
                                   time.window = c(0,15,35,65,90,120),
                                   names.window = names.window)
                    ```
                    - OBS: some possible plots with ggplot2....
                    ```{r}
                    
                    
                    ```
                    - For more than one variable, we can use the quantiles for all environments
                    ```{r}
                    EnvTyping(env.data = df.clim,var.id =  c('T2M','PRECTOT','WS2M'),env.id='env',by.interval = T)
                    ```
                    - We can define the cardinals for each variable
                    ```{r}
                    (cardinals= list(T2M=c(0,9,22,32,45),PRECTOT=c(0,5,10),WS2M=c(0,1,5)))
                    
                    EnvTyping(env.data = df.clim,var.id =  c('T2M','PRECTOT','WS2M'),
                              cardinals = cardinals,env.id='env')
                    
                    ```
                    - However, we do not always have ecophysiological information about the best possible cardinals ... so we use quantiles!
                      If quantiles = NULL, 1%, 25%, 50%, 99% is assumed
                    ```{r}
                    (cardinals= list(T2M=c(0,9,22,32,45),PRECTOT=c(0,5,10),WS2M=NULL))
                    EnvTyping(env.data = df.clim,var.id =  c('T2M','PRECTOT','WS2M'),
                              cardinals = cardinals,env.id='env')
                    ```
                    - All analyses can also be run considering centered on the mean and scaled x ~ N (0.1)
                    ```{r}
                    EnvTyping(env.data = df.clim,var.id = 'PRECTOT',env.id='env',scale = T)
                    EnvTyping(env.data = df.clim,var.id =  c('T2M','PRECTOT','WS2M'),env.id='env',scale = T) 
                    ```
                    [Menu](#menu)
                      
                      <div id="P3" />
                        
                        ------------------------------------------------------------
                        
                        # 3. Reaction Norm Module
                        
                        We provide Genomic and Envirotypic kernels for reaction norm prediction. After generate the kernels, the user must use the [BGGE](https://github.com/italo-granato/BGGE) package to run the models
                      
                      - Toy Example: genomic prediction for grain yield in tropical maize
                      ```{r}
                      data("maizeYield") # 150 maize hybrids over 5 environments (grain yield data)
                      data("maizeG")     # GRM for maizeYield
                      data('maizeWTH')   # weather data for maize Yield
                      
                      Y <- maizeYield
                      G <- maizeG
                      df.clim <- maizeWTH
                      
                      ```
                      - Statistical Models
                      
                      <p align="center">
                        <img src="/fig/summary_pkg.png" width="70%" height="70%">
                        </p>
                        
                        - Returns benchmark main effect model: 
                        
                        <p align="center">
                        <img width="120" height="18" src="/fig/mod1.png">
                        </p>
                        
                        ```{r}
                      MM <- get_kernel(K_G = list(G=G),Y = Y,reaction = F,model = 'MM')
                      ```
                      - Returns benchmark main GxE deviation model:
                        
                        <p align="center">
                        <img width="160" height="18" src="/fig/mod2.png">
                        </p>
                        
                        ```{r}
                      MDs <-get_kernel(K_G = list(G=G),Y = Y,reaction = F,model = 'MDs')
                      ```
                      - Obtaining environmental variables based on quantiles
                      
                      ```{r}
                      W.cov<-W.matrix(env.data = df.clim,by.interval = T,statistic = 'quantile',
                                      time.window = c(0,14,35,60,90,120))
                      dim(W.cov)
                      
                      ```
                      - Creating Env Kernels from W matrix and Y dataset
                      
                      ```{r}
                      H <- EnvKernel(env.data = W.cov,Y = Y,merge = T,env.id = 'env')
                      dim(H)
                      dim(H$varCov) # variable relationship
                      dim(H$envCov) # environmental relationship
                      
                      #env.plots(H$envCov,row.dendrogram = T,col.dendrogram = T) # superheat
                      superheat(H$envCov,row.dendrogram = T,col.dendrogram = T)
                      
                      ```
                      - Parametrization by 
                      
                      <p align="center">
                        <img width="110" height="50" src="/fig/mod3.png">
                        </p>
                        
                        ```{r}
                      H <- EnvKernel(env.data = W.cov,Y = Y,merge = T,env.id = 'env',bydiag=FALSE)
                      dim(H)
                      dim(H$varCov) # variable relationship
                      dim(H$envCov) # environmental relationship
                      
                      #env.plots(H$envCov,row.dendrogram = T,col.dendrogram = T) # superheat
                      superheat(H$envCov,row.dendrogram = T,col.dendrogram = T)
                      
                      ```
                      - Parametrization by 
                      
                      <p align="center">
                        <img width="130" height="50" src="/fig/mod4.png">
                        </p>
                        
                        resulting in diag(K_W) = 1
                      
                      ```{r}
                      H <- EnvKernel(env.data = W.cov,Y = Y,merge = T,env.id = 'env',bydiag=TRUE)
                      dim(H)
                      dim(H$varCov) # variable relationship
                      dim(H$envCov) # environmental relationship
                      
                      #env.plots(H$envCov,row.dendrogram = T,col.dendrogram = T) # superheat
                      superheat(H$envCov,row.dendrogram = T,col.dendrogram = T)
                      
                      ```
                      - Gaussian parametrization by 
                      
                      <p align="center">
                        <img width="130" height="50" src="/fig/mod5.png">
                        </p>
                        
                        which d = dist(W), q = median(d) and h = gaussian parameter (default = 1)
                      
                      ```{r}
                      H <- EnvKernel(env.data = W.cov,Y = Y,merge = T,env.id = 'env',gaussian=TRUE)
                      dim(H)
                      dim(H$varCov) # variable relationship
                      dim(H$envCov) # environmental relationship
                      
                      #env.plots(H$envCov,row.dendrogram = T,col.dendrogram = T) # superheat
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
W.cov<-W.matrix(env.data = df.clim,by.interval = T,statistic = 'quantile',
                time.window = c(0,14,35,60,90,120))

W <- EnvKernel(env.data = W.cov,Y = Y,merge = T,env.id = 'env',bydiag=TRUE)

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


<div id="P4" />

------------------------------------------------------------

## 5. Authorship

This package is a initiative from the [Allogamous Plant Breeding Lab (University of São Paulo, ESALQ/USP, Brazil)](http://www.genetica.esalq.usp.br/en/lab/allogamous-plant-breeding-laboratory).

**Developer**

> * [Germano Costa Neto](https://github.com/gcostaneto), PhD Candidate in Genetics and Plant Breeding


<div id="P5" />

------------------------------------------------------------

## 6. Acknowledgments

> * [Giovanni Galli](https://github.com/giovannigalli), PhD in Genetics and Plant Breeding

> * [Humberto Fanelli](https://github.com/humbertofanelli), PhD in Genetics and Plant Breeding

> * [Roberto Fritsche-Neto](roberto.neto@usp.br), PhD Candidate in Genetics and Plant Breeding

> * [University of São Paulo (ESALQ/USP)](https://www.esalq.usp.br/)

> * [Conselho Nacional de Desenvolvimento Científico e Tecnológico](http://www.cnpq.br/) for the PhD scholarship granted to the authors of the package

> * [Pedro L. Longhin](https://github.com/pedro-longhin) for additional support in Git Hub

<div id="P6" />

------------------------------------------------------------

## 7. Simplified R script (Tutorial)

[Simplified Tutorial (R script)](https://raw.githubusercontent.com/allogamous/EnvRtype/master/tutorial_script_R.R)



<img align="right" width="110" height="100" src="/fig/logo_alogamas.png">


<div id="menu" />

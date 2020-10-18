## Remote Environmental Sensing
 
 ### Geographic Information Databases
 
 > The collection, organization and processing of environmental data is a step that requires equipment installed in the field. Such equipment can be expensive or difficult to access for some research groups in certain regions or countries. For this reason, we decided to insert a routine for collecting climatic data through the [NASA POWER base](https://power.larc.nasa.gov/), which can access information on a daily scale anywhere on the globe.
 
 > The [Raster Package](https://cran.r-project.org/web/packages/raster/raster.pdf) also offers a digital platform for downloading files in raster format of climatic data (from the [WorldClim database](https://www.worldclim.org/)) and [SRTM (elevation)](http://srtm.csi.cgiar.org/) using only geographical coordinates (Latitude and Longitude).
 
  ### Remote Data Collection
 
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

> The argument *country* were setted to collect elevation data from Brazil (BRA). For other countries please consult the ISO in the following table:

<img align="center" src="/fig/ISO3.png" width="90%" height="90%">

<!-- toc -->
[Menu](#menu)


 ### Data Processing
 
**Additional variables (ecophysiological)**
 
 
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


<div id="P2" />
                      
## Environmental Characterization Module

> *  Environmental characterization is a fundamental step to understand how the environment regulates the phenotypic expression and adaptation of the genotypes under different growing conditions. For this reason, based on envirotyping (environmental + typing) studies alredy published (see references), we provide mechanisms that enable the typing of environmental factors in terms of frequency of occurrence. We have also developed functions for collecting environmental factors and organizing them as covariates to be used in reaction norm studies.


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
> *  We can perform a Quality Control (QC) based on the maximum sd tolered.

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


### Environmental Typologies

> *  The *EnvTyping* function provides mechanisms for organizing environmental types through the selection of specific environmental factors and their analysis over time and space
                      
```{r}
EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M')
```

> * Typologies can be defined across different time intervals by setting the argument *by.interval = TRUE* (generic time intervals)

```{r}
EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T)
```

> * Typologies can be defined across **specific** time intervals by setting the argument *by.interval = TRUE* and and defining the time windows (in days after begining of the data) using the function *time.window*

```{r}
EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T,time.window = c(0,15,35,65,90,120))
```
> * Typologies by.intervals (specific time intervals and with specific names)

```{r}
names.window <- c('1-intial growing','2-leaf expansion I','3-leaf expansion II','4-flowering','5-grain filling','6-maturation')
time.window  <- c(0,15,35,65,90,120)
EnvTyping(env.data = df.clim,env.id = 'env',var.id='T2M',by.interval = T,time.window = time.window, names.window = names.window)
```
> * The two-way table of typologies can be plotted based on this [code](https://raw.githubusercontent.com/allogamous/EnvRtype/master/plot.R)

### **Example 1 : Air Temperature (**T2M**)**


**Option 1: facet by developmental stages**


 <p align="center">
 <img width="100%" height="100%" src="/fig/t2m_full.png">
  </p>


**Option 2: envirotypes as a combination of environmental factor x cardinal class x developmental stage**


 <p align="center">
 <img width=50%" height="50%" src="/fig/t2m_freq.png">
  </p>
  
  
  **Option 3: envirotypes per environnment**


 <p align="center">
 <img width=100%" height="100%" src="/fig/per_env_t2m.png">
  </p>
  
  
  

### **Example 2 : Solar Radiation (**SRAD**)**


**Option 1: facet by developmental stages**


 <p align="center">
 <img width="100%" height="100%" src="/fig/srad_full.png">
  </p>


**Option 2: envirotypes as a combination of environmental factor x cardinal class x developmental stage**


 <p align="center">
 <img width=50%" height="50%" src="/fig/srad_freq.png">
  </p>
  

  **Option 3: envirotypes per environnment**


 <p align="center">
 <img width=100%" height="100%" src="/fig/per_env_srad.png">
  </p>
                                                        
### **Example 3 : Valour Pressure Defict (**VPD**)**


**Option 1: facet by developmental stages**


 <p align="center">
 <img width="100%" height="100%" src="/fig/vpd_full.png">
  </p>


**Option 2: envirotypes as a combination of environmental factor x cardinal class x developmental stage**


 <p align="center">
 <img width=50%" height="50%" src="/fig/vpd_freq.png">
  </p>
  
  
   **Option 3: envirotypes per environnment**


 <p align="center">
 <img width=100%" height="100%" src="/fig/per_env_vpd.png">
  </p>
  
> * in *var.id* you can put the names of the desirable variables:

```{r}
EnvTyping(env.data = df.clim,var.id =  c('T2M','PRECTOT','WS2M'),env.id='env',by.interval = T)
```
> * Then, we can define the different cardinals for each variable

```{r}
# Create a list of cardinals
cardinals <- list(T2M=c(0,9,22,32,45),PRECTOT=c(0,5,10),WS2M=c(0,1,5))
EnvTyping(env.data = df.clim,var.id =  c('T2M','PRECTOT','WS2M'),cardinals = cardinals,env.id='env')
```

> * These cardinals must respect ecophysiological limits for each crop, germplasm or region. For that, we recommend looking for ecophysiology literature and crop growth modeling, such as Soltani and Sinclar (2012) [**Modeling physiology of crop development, growth and yield**](https://www.amazon.com.br/Modeling-Physiology-Development-Growth-Yield/dp/1845939700); However, we do not always have ecophysiological information about the best possible cardinals ... so we use quantiles!

> * If quantiles = NULL, 1%, 25%, 50%, 99% is assumed
```{r}
cardinals= list(T2M=c(0,9,22,32,45),PRECTOT=c(0,5,10),WS2M=NULL)
EnvTyping(env.data = df.clim,var.id =  c('T2M','PRECTOT','WS2M'),cardinals = cardinals,env.id='env')
```

> * All analyses can also be run considering centered on the mean and scaled x ~ N (0.1)

```{r}
EnvTyping(env.data = df.clim,var.id = 'PRECTOT',env.id='env',scale = T)
EnvTyping(env.data = df.clim,var.id =  c('T2M','PRECTOT','WS2M'),env.id='env',scale = T) 
```

[Menu](#menu)

<div id="P3" />
                        
------------------------------------------------------------


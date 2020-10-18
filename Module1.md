# **Remote Environmental Sensing**

* [Geographic Information Databases](#P1)
* [Remote Data Collection](#P2)
* [Raw Data Processing](#P3)
* [Raw Data Summary](#P4)

              
<div id="P1" />

 
 ### Geographic Information Databases
 
 > The collection, organization and processing of environmental data is a step that requires equipment installed in the field. Such equipment can be expensive or difficult to access for some research groups in certain regions or countries. For this reason, we decided to insert a routine for collecting climatic data through the [NASA POWER base](https://power.larc.nasa.gov/), which can access information on a daily scale anywhere on the globe.
 
 > The [Raster Package](https://cran.r-project.org/web/packages/raster/raster.pdf) also offers a digital platform for downloading files in raster format of climatic data (from the [WorldClim database](https://www.worldclim.org/)) and [SRTM (elevation)](http://srtm.csi.cgiar.org/) using only geographical coordinates (Latitude and Longitude).
 
 
 <div id="P2" />
 
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


<div id="P3" />

 ### Raw Data Processing
 
**Additional variables (ecophysiological)**
 
 
> * Basic processing of get_weather() 


 ```{r}
 df.clim <-processWTH(env.data = df.clim)
```
<div id="P4" />

### Raw Data Summary
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
[Menu](https://github.com/allogamous/EnvRtype)



------------------------------------------------------------


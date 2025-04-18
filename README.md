<img align="left" src="/fig/logo3_300.png" width="20%" height="20%">

#  [*EnvRtype*: Envirotyping Tools in R](https://github.com/allogamous/EnvRtype/blob/master/README.md)

### A R Interplay between Quantitative Genetics and Ecophysiology for GxE analysis


#### **Current Version**: 1.1.2 (Mar 2025 | Last Version:  1.1.1  (Oct 2024)

Please consider migrating to [envirotypeR repo](https://github.com/gcostaneto/envirotypeR)
<div id="menu" />
  
  [![DOI](https://img.shields.io/badge/DOI-doi.org%2F10.1093%2Fg3journal%2Fjkab040-orange)](https://doi.org/10.1093/g3journal/jkab040)
  [![SUPPORT](https://img.shields.io/badge/SUPPORT-R-yellowgreen)](https://github.com/gcostaneto/EnvRtype_course/blob/main/README.md)
 





##





Envirotyping has proven useful in identifying the non-genetic drivers of phenotypic adaptation in plants cultivaded in diverse growing conditions. Combined with phenotyping and genotyping data, the use of envirotyping data may leverage the molecular breeding strategies to cope with environmental changing scenarios. Over the last ten years, this data has been incorporated in genomic-enabled prediction models aiming to better model genotype x environment interaction (GE) as a function of reaction-norm. However, there is difficult for most breeders to deal with the interplay between envirotyping, ecophysiology, and genetics. 
  
It also can be useful for several fields of agricultural, livestook and ecology research, by delivering high-quality environmental information and environmental grouping appraoches.

Here we present the EnvRtype R package as a new toolkit developed to facilitate the interplay between envirotyping and fields of plant research such as genomic prediction. This package offers three modules: (1) collection and processing data set, (2) environmental characterization, (3) build of ecophysiological enriched predictive models accounting for three different structures of reaction-norm over different sources of genomic relatedness. Thus, EnvRtype is useful for exploratory purposes and predctive breeding for multiple growing conditions.

<div id="menu" />
  

  ## Resources
  
 The envirotyping pipeline provided by EnvRtype consists in three modules (1 - Environmental Sensing, 2- Macro-Environmental Characterization and 3 - Enviromic Similarity and Phenotype Prediction). Collectively, the EnvRtyping functions generate a simple workflow to collect, process and integrates envirotyping data into several fields of agricultural research, specially for predictive breeding that may include the use of genomic x enviromic relatedness information.
  
  <img align="center" src="/fig/workflow_2.png" width="90%" height="90%">
  
 

  ## Updates and Maintence

* EnvRtype 1.1.2 adding get_soil() function to download data from SoilGrids


* EnvRtype 1.1.2 is online (Mar 2025) - correcting bug on get_weather()

* EnvRtype 1.1.1 is online (17 Oct 2024)
  
* raster bug fixed (when downloading altitude data)
  
* rgdal is retired! We fixed it by Oct 21 2023. More info about rgdal's retirement [here](https://r-spatial.org/r/2022/04/12/evolution.html)
  
* EnvRtype 1.1.0 is online (1st June 2022)

* New version of the get_weather() function -- now running in parallel (Jun 2022)

* Coming Soon: get_soil() function and EPA() [environmental-phenotype associations]
* [SOLVED] Date 2021-10-14: NASA POWER server off "Error: Something went wrong with the query, no data were returned. Please see <https://power.larc.nasa.gov> for potential server issues."
* Join our [DISCUSSION FORUM](https://groups.google.com/u/1/g/envrtype)
* [SOLVED] PRECTOT variable (rainfall precipitation) is currently off from NASA POWER
* Coming soon: tutorial for using external sources of environmental data (from field micro-stations) 
* [SOLVED] Coming soon (Dec 2021): tutorial for colecting soil data from SoilGrids data base
* The current version of the package is 0.0.2 (May 22th 2021)
* [SOLVED] From December 15th 2020 to January 10th 2021 this page will be under maintence. This means that we are now working in several updates and some changes will be made in some functions.
* The current version of the package is 0.0.1 (Nov 20th 2020)


 ## Tutorials

* [Envirotyping pipeline](https://github.com/allogamous/EnvRtype/blob/master/Enviromic_pipeline.md)
* [Genomic Prediction using Environmental Covariates](https://github.com/allogamous/EnvRtype/blob/master/Genomic%20Prediction.md)
* [Full examples and R Codes for the G3 Paper](https://github.com/allogamous/EnvRtype/blob/master/EnvRtype_full_tutorial.R)

**Information**
* [Authorship](#P4)
* [Acknowledgments](#P5)
* [Publications](#P6)
* [Getting help](https://groups.google.com/u/1/g/envrtype)

              
<div id="Instal" />
                
## Install

### Using devtools in R

```{r}
library(devtools)
devtools::install_github('allogamous/EnvRtype',force=TRUE) # current version:  1.1.0 (June 2022)
# Enter one or more numbers, or an empty line to skip updates: 3
require(EnvRtype)
  ```
### Manually installing

> If the method above doesn't work, use the next lines by downloading the EnvRtype-master.zip file

```{r}
setwd("~/EnvRtype-master.zip") # ~ is the path from where you saved the file.zip
unzip("EnvRtype-master.zip") 
file.rename("EnvRtype-master", "EnvRtype") 
shell("R CMD build EnvRtype") # or system("R CMD build EnvRtype")
install.packages("EnvRtype_1.1.0.tar.gz", repos = NULL, type="source") # Make sure to use the current verision
```
 
 ### Required packages
 
 For some users, it seems that the packages below must downloaded...(I am not a IT guy, I am really don't know why, sorry).
 
 * **[EnvRtype](https://github.com/allogamous/EnvRtype)** 
 * **[raster](https://CRAN.R-project.org/package=raster)** 
 * **[nasapower](https://github.com/ropensci/nasapower)** 
 * **[BGGE](https://github.com/italo-granato/BGGE)**
 * **[foreach](https://github.com/cran/foreach)**
 * **[doParalell](https://github.com/cran/doparallel)**
                
```{r}
install.packages("foreach")
install.packages("doParallel")
install.packages("raster")
install.packages("nasapower")
install.packages("rgdal")
install.packages("BGGE")
              
or
              
#source("https://raw.githubusercontent.com/gcostaneto/Funcoes_naive/master/instpackage.R");
#inst.package(c("BGGE",'foreach','doParallel','raster','rgdal','nasapower'));

install.packages(c("BGGE",'foreach','doParallel','raster','rgdal','nasapower'))

library(EnvRtype)
              
```
<!-- toc -->
[Menu](#menu)
                  
 <div id="P1" />
  


## Authorship

This package is a initiative from the [Allogamous Plant Breeding Lab (University of São Paulo, ESALQ/USP, Brazil)](http://www.genetica.esalq.usp.br/en/lab/allogamous-plant-breeding-laboratory).

**Developer**

 * [Germano Costa Neto](https://github.com/gcostaneto), University of Sao Paulo/ Cornell University


**Maintence**

 * [Germano Costa Neto](https://github.com/gcostaneto)



<div id="P5" />

## Publications

*Last update: 2021-10-10*

* Costa-Neto, G., Crossa, J., and Fritsche-Neto, R. (2021). Enviromic Assembly Increases Accuracy and Reduces Costs of the Genomic Prediction for Yield Plasticity in Maize. **Frontiers in Plant Science** 12. doi:10.3389/fpls.2021.717552.

* Costa-Neto, G., Galli, G., Carvalho, H. F., Crossa, J., and Fritsche-Neto, R. (2021). EnvRtype: a software to interplay enviromics and quantitative genomics in agriculture. **G3 Genes|Genomes|Genetics**. doi:10.1093/g3journal/jkab040.

* Galli G, Horne DW, Collins SD, Jung J, Chang A, Fritsche‐Neto R, et al. (2020). Optimization of UAS‐based high‐throughput phenotyping to estimate plant health and grain yield in sorghum. **Plant Phenome** J 3: 1–14.

* Costa-Neto G, Fritsche-Neto R, Crossa J (2020). Nonlinear kernels, dominance, and envirotyping data increase the accuracy of genome-based prediction in multi-environment trials. **Heredity** (Edinb).

  
  <div id="P7" />

## Acknowledgments

 * [Giovanni Galli](https://github.com/giovannigalli)

 * [Humberto Fanelli](https://github.com/humbertofanelli)

 * Jose Crossa, Biometrics and Statistic Unit at CIMMYT.

 * [Roberto Fritsche-Neto](roberto.neto@usp.br)

 * [University of São Paulo (ESALQ/USP)](https://www.esalq.usp.br/)

 * [Conselho Nacional de Desenvolvimento Científico e Tecnológico](http://www.cnpq.br/) for the PhD scholarship granted to the authors of the package

 * [Pedro L. Longhin](https://github.com/pedro-longhin) for additional support in Git Hub

<div id="P6" />
  


<img align="right" width="110" height="100" src="/fig/logo_alogamas.png">


<div id="menu" />


<div align='center'>

<a href='https://www.free-website-hit-counter.com'><img src='https://www.free-website-hit-counter.com/c.php?d=9&id=159093&s=1' border='0' alt='Free Website Hit Counter'></a><br / ><small><a href='https://www.free-website-hit-counter.com' title="Free Website Hit Counter">Free website hit counter</a></small>

</div>


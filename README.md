<img align="left" src="/fig/logo3_300.png" width="20%" height="20%">

#  [*EnvRtype*: Envirotyping Tools in R](https://github.com/allogamous/EnvRtype/blob/master/README.md)

### A R Interplay between Quantitative Genetics and Ecophysiology for GxE analysis


#### **Current Version of this repo**: 1.1.2 (Mar 2025 | Last Version:  1.1.1  (Oct 2024)

Please consider migrating to new repo and updated package (Latest version, Oct 2026, and now in CRAN!) [EnvRtype](https://github.com/gcostaneto/EnvRtype)
<div id="menu" />
  
  [![DOI](https://img.shields.io/badge/DOI-doi.org%2F10.1093%2Fg3journal%2Fjkab040-orange)](https://doi.org/10.1093/g3journal/jkab040)
  [![SUPPORT](https://img.shields.io/badge/SUPPORT-R-yellowgreen)](https://github.com/gcostaneto/EnvRtype_course/blob/main/README.md)
 





##



Please migrate to the new repo [EnvRtype 0.1.2](https://github.com/gcostaneto/EnvRtype)

Envirotyping has proven useful in identifying the non-genetic drivers of phenotypic adaptation in plants cultivaded in diverse growing conditions. Combined with phenotyping and genotyping data, the use of envirotyping data may leverage the molecular breeding strategies to cope with environmental changing scenarios. Over the last ten years, this data has been incorporated in genomic-enabled prediction models aiming to better model genotype x environment interaction (GE) as a function of reaction-norm. However, there is difficult for most breeders to deal with the interplay between envirotyping, ecophysiology, and genetics. 
  
It also can be useful for several fields of agricultural, livestook and ecology research, by delivering high-quality environmental information and environmental grouping appraoches.

Here we present the EnvRtype R package as a new toolkit developed to facilitate the interplay between envirotyping and fields of plant research such as genomic prediction. This package offers three modules: (1) collection and processing data set, (2) environmental characterization, (3) build of ecophysiological enriched predictive models accounting for three different structures of reaction-norm over different sources of genomic relatedness. Thus, EnvRtype is useful for exploratory purposes and predctive breeding for multiple growing conditions.

<div id="menu" />
  

  ## Resources
  
 The envirotyping pipeline provided by EnvRtype consists in three modules (1 - Environmental Sensing, 2- Macro-Environmental Characterization and 3 - Enviromic Similarity and Phenotype Prediction). Collectively, the EnvRtyping functions generate a simple workflow to collect, process and integrates envirotyping data into several fields of agricultural research, specially for predictive breeding that may include the use of genomic x enviromic relatedness information.
  
  <img align="center" src="/fig/workflow_2.png" width="90%" height="90%">
  
 

 ## Tutorials

Please migrate to the new repo [EnvRtype 0.1.2](https://github.com/gcostaneto/EnvRtype)


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


```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("gcostaneto/EnvRtype")
```

```r
library(EnvRtype)
```
### Manually installing

Please migrate to the new repo [EnvRtype 0.1.2](https://github.com/gcostaneto/EnvRtype)
 
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


cat("it is ugly, but works! :) \n")

if (!requireNamespace('nasapower', quietly = TRUE)) {utils::install.packages("nasapower")} 
if (!requireNamespace('plyr', quietly = TRUE)) {utils::install.packages("plyr")}
if (!requireNamespace('BGGE', quietly = TRUE)) {utils::install.packages("BGGE")}
if (!requireNamespace('foreach', quietly = TRUE)) {utils::install.packages("foreach")}
if (!requireNamespace('reshape2', quietly = TRUE)) {utils::install.packages("reshape2")}
if (!requireNamespace('raster', quietly = TRUE)) {utils::install.packages("raster")}

source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/AtmosphericPAram.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/SradPARAM.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/SupportFUnction.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/EnvTyping.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/Wmatrix.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/covfromraster.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/envKernel.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/gdd.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/getGEenriched.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/get_weather_gis.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/processWTH.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/met_kernel_model.R')
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/R/summary_weather.R')

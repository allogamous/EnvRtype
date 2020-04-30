#'##############################################################################
#' Basic tutorial for EnvRtype package
#' gcneto and rfn, april 08 2020
#'##############################################################################

############# Installing the package  ##################

library(devtools)
install_github('allogamous/EnvRtype')

require(EnvRtype)

############# Collecting daily weather from GIS data sources  ##################
lat = c(-13.05,-12.32,-18.34,-18.90,-23.03) # vector of latitude WGS84
lon = c(-56.05,-55.42,-46.31,-49.56,-51.02) # vector of lontitude WGS84
env = c("NM","SO","PM","IP","SE")           # vector of environment/site ID
plant.date = c("2015-02-15",'2015-02-13', # vector with the beggining of growings seasons
               "2015-02-26","2015-03-01",
               '2015-02-19')
harv.date = rep("2015-06-15", 5) # vector of end period

df.clim <- get_weather(env.id = env, lat = lat, lon = lon,
                       start.day = plant.date, end.day = harv.date, country = 'BRA')

head(df.clim)

#if this does not work try:
#df.clim = EnvRtype::maizeWTH[,-(21:28)]

############# Basic processing of get_weather() output ###############################

df.clim2 <- processWTH(weather.data = df.clim)

############# Basic summary statistics of environmental data ##################

summaryWTH(df.clim)
summaryWTH(df.clim, env.id = c('env')) #atenção o argumento env nao esta afetando o resultado

### Summary a particular environmental variable
summaryWTH(df.clim,env.id = 'env', var.id = 'T2M')
summaryWTH(df.clim,env.id = 'env', var.id = c('T2M', 'T2M_MAX'))

### Summary by time intervals (e.g., phenology)
summaryWTH(df.clim, env.id = 'env', by.interval = T)

### Summary by time intervals given by time.window
summaryWTH(df.clim, env.id = 'env', by.interval = T, time.window = c(0, 14, 35, 60, 90, 120))

### Summary by time intervals given by time.window and names.window
summaryWTH(df.clim, env.id = 'env', by.interval = T,
           time.window = c(0, 14, 35, 60, 90, 120),
           names.window = c('P-E', 'E-V1', 'V1-V4', 'V4-VT', 'VT-GF', 'GF-PM'))


### returning only mean values
summaryWTH(df.clim, env.id = 'env', statistic = 'mean')

### or only sum values
summaryWTH(df.clim, env.id = 'env', statistic = 'sum')

### or quantile values (default = 25%, 50% and 75%)
summaryWTH(df.clim, env.id = 'env', statistic = 'quantile')

### or specific quantiles (e.g., 20%, 76% and 90%)
summaryWTH(df.clim, env.id = 'env', statistic = 'quantile', probs = c(.20, .76, .90))


############### Building the Environmental Covariable Matrix  ##################

### Mean-centered and scaled matrix
W <- W.matrix(weather.data = df.clim, by.interval = F)
dim(W)
W[,1:5]

## same as SummaryWTH, we can add time.windows
W <- W.matrix(weather.data = df.clim, by.interval = T,
         time.window = c(0, 14, 35, 60, 90, 120))
dim(W)
W[,1:5]

## and select the statistic to be used
W <- W.matrix(weather.data = df.clim, by.interval = T, statistic = 'mean',
         time.window = c(0, 14, 35, 60, 90, 120))
dim(W)
W[,1:5]

W <- W.matrix(weather.data = df.clim, by.interval = T, statistic = 'quantile',
         time.window = c(0, 14, 35, 60, 90, 120))
dim(W)
W[,1:5]

## We can perform a Quality Control (QC) based on the maximum sd tolerated
W <- W.matrix(weather.data = df.clim, by.interval = F, QC = T)
dim(W)
W[,1:5]

## We can perform a Quality Control (QC) based on the maximum sd tolered
W <- W.matrix(weather.data = df.clim, by.interval = F, QC = T, sd.tol = 3)
dim(W)
W[,1:5]

## or create W for specific variables
id.var = c('T2M_MAX','T2M_MIN','T2M')
W <- W.matrix(weather.data = df.clim, var.id = id.var)
dim(W)
W[,1:3]

## or combine with summaryWTH by using is.processed = T
data <- summaryWTH(df.clim, env.id = 'env', statistic = 'quantile')
W <- W.matrix(weather.data = data, is.processed = T)
dim(W)
W[,1:3]


################ Environmental Typologies based on Cardinal Limits #############

EnvTyping(weather.data = df.clim, env.id = 'env', var.id = 'T2M')

# or by.intervals (generic time intervals)
EnvTyping(weather.data = df.clim, env.id = 'env', var.id = 'T2M', by.interval = T)

# or by.intervals (specific time intervals)
EnvTyping(weather.data = df.clim,
          env.id = 'env',
          var.id = 'T2M',
          by.interval = T,
          time.window = c(0, 15, 35, 65, 90, 120))

# or by.intervals (specific time intervals and with specific names
names.window = c('1-intial growing',
                 '2-leaf expansion I',
                 '3-leaf expansion II',
                 '4-flowering',
                 '5-grain filling',
                 '6-maturation')

out <- EnvTyping(weather.data = df.clim,
                 env.id = 'env',
                 var.id = 'T2M',
                 by.interval = T,
                 time.window = c(0, 15, 35, 65, 90, 120),
                 names.window = names.window)


######## Visualing the output #######

# plot 1: enviromental variables panel
require(ggplot2)
ggplot() +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradientn(colours= rainbow(15))+
  geom_tile(data=out,aes(x=reorder(env,Freq), y=reorder(env.variable,Freq),fill=Freq))+
  ylab('Envirotype ID\n')+
  xlab("\nEnvironment ID")+
  labs(fill='Frequency')+
  theme(axis.title = element_text(size=19),
        legend.text = element_text(size=9),
        strip.text.y  = element_text(size=13,angle=360),
        legend.title = element_text(size=17),
        strip.background = element_rect(fill="gray95",size=1),
        legend.position = 'bottom')


# plot 2: distribution of envirotypes
ggplot() +
  #theme_void()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # theme_pubclean()+
  #scale_fill_manual(values=c('red2','orange2','green3','blue3','violet'))+
  facet_grid(var~interval)+
  geom_bar(data=out, aes(y=Freq, x=env,fill=env.variable),
           position = "fill",stat = "identity",width = 1,size=.2)+
  # scale_y_continuous(labels = scales::percent,expand = c(0,0))+ #coord_flip()+
  ylab('Absolute Frequency\n of Occurence\n')+
  xlab("\nEnvironment ID")+
  labs(fill='Envirotype')+
  theme(axis.title = element_text(size=19),
        #   axis.text = element_blank(),
        legend.text = element_text(size=9),
        strip.text = element_text(size=17),
        legend.title = element_text(size=17),
        strip.background = element_rect(fill="gray95",size=1),
        legend.position = 'bottom')


# For more than one variable, we can use the quantiles for all environments
EnvTyping(weather.data = df.clim,
          var.id = c('T2M','PRECTOT','WS2M'),
          env.id = 'env',
          by.interval = T)

# Also, we can define the cardinals for each variable
(cardinals = list(T2M = c(0, 9, 22, 32, 45),
                  PRECTOT = c(0, 5, 10),
                  WS2M = c(0, 1, 5)))

EnvTyping(weather.data = df.clim,
          var.id = c('T2M','PRECTOT','WS2M'),
          cardinals = cardinals,
          env.id = 'env')

# However, in most cases we know litle about cardinals, Thus, for variables without a cardinal, NULL is assumed, therefore, the quantile delimited in quantiles () is run.
# If quantiles = NULL, 1%, 25%, 50%, 99% is assumed

(cardinals = list(T2M = c(0, 9, 22, 32, 45),
                  PRECTOT = c(0, 5, 10),
                  WS2M = NULL))

EnvTyping(weather.data = df.clim,
          var.id = c('T2M','PRECTOT','WS2M'),
          cardinals = cardinals,
          env.id = 'env')

# all analyses can also be run considering centered on the mean and scaled x ~ N (0.1)
EnvTyping(weather.data = df.clim, var.id = 'PRECTOT', env.id = 'env', scale = T)

EnvTyping(weather.data = df.clim, var.id = c('T2M','PRECTOT','WS2M'),env.id = 'env', scale = T)


################  Genomic Prediction-based reaction norm models   ###################
#' We provide Genomic and Envirotypic kernels for reaction norm prediction
#' After generate the kernels, the user must use the BGGE package to run the predictions

# loading the G matrix of hybrids
data('maizeG')
dim(maizeG)
maizeG[1:5, 1:5]

# loading the phenotypes of hybrids
data('maizeYield')
dim(maizeYield)
head(maizeYield)

# loading weather data
data("maizeWTH")
dim(maizeWTH)
maizeWTH[1:6, 1:12]


# obtaining the benchmark main effect model: Y = fixed + G using BGGE
MM <- get_kernel(K_G = list(G = as.matrix(maizeG)),Y = maizeYield, model = 'MM')

# or the benchmark main GxE deviation model: Y = fixed + G +GE using BGGE
MDs <- get_kernel(K_G = list(G = as.matrix(maizeG)),Y = maizeYield, model = 'MDs')

# Creating the envirotyping covariance matrix based on phenology time intervals
W.cov <- W.matrix(weather.data = maizeWTH,by.interval = T, statistic = 'quantile',
                  time.window = c(0, 14, 35, 60, 90, 120))
dim(W.cov)
W.cov[,1:4]


# Creating Env Kernels from W matrix and Y dataset
H <- EnvKernel(weather.data = W.cov, Y = maizeYield, merge = T, env.id = 'env')
class(H)

dim(H$varCov) # variable relationship
H$varCov[1:4,1:4]

dim(H$envCov) # environmental relationship, Parametrization by K_W = WW'/ncol(W)
H$envCov[1:4,1:4]

# Visualising the relationssip between environments and phenological stages
require(superheat)
superheat(H$envCov, row.dendrogram = T, col.dendrogram = T)

# OBS: others parametrizations for envCov (W)

## benchmark: Parametrization by K_W = WW'/ncol(W)
EK <- EnvKernel(weather.data = W.cov,
                Y = maizeYield,
                merge = F,
                env.id = 'env',
                gaussian = F)$envCov

dim(EK)
EK


## Parametrization by K_W = WW'/diag( WW'), resulting in diag(K_W) = 1
EK <- EnvKernel(weather.data = W.cov,
                Y = maizeYield,
                merge = F,
                env.id = 'env',
                bydiag = T)$envCov
dim(EK)
EK

## Gaussian parametrization by K_W = exp(-hd/q), which d = dist(W), q = median(d) and h = gaussian parameter (default = 1)
EK <- EnvKernel(weather.data = W.cov,
                Y = maizeYield,
                merge = F,
                env.id = 'env',
                gaussian = T,
                h.gaussian = )$envCov
dim(EK)
EK

# Finally, creating the Genomic-based prediction enriched kernnels
# K_G = list of genomic kernels
# K_E = list of environmental kernels
# reaction = TRUE, build the haddamard's product between genomic and envirotype-based kernels
# reaction = FALSE, but K_E != NULL, only random environmental effects using K_E are incorporated in the model

# returning the benchmark main effect model: Y = fixed + W + G
EMM <- get_kernel(K_G = list(G = as.matrix(maizeG)),
                  Y = maizeYield,K_E = list(W = H$envCov),
                  model = 'EMM') # or model = MM

names(EMM)

# or the benchmark main GxE deviation model: Y = fixed + G + W + GE
EMDs <- get_kernel(K_G = list(G = as.matrix(maizeG)),
                   Y = maizeYield,
                   K_E = list(W = H$envCov),
                   model = 'MDs') # or model = MDs
names(EMDs)


# or the benchmark main effect model: Y = fixed + W + G + GW
RN <- get_kernel(K_G = list(G = as.matrix(maizeG)),
                 Y = maizeYield,
                 K_E = list(W = H$envCov),
                 model = 'RNMM')
names(RN)


# or the benchmark main effect model: Y = fixed + W + G + GW + GE
fullRN <- get_kernel(K_G = list(G = as.matrix(maizeG)),
                     Y = maizeYield,
                     K_E = list(W = H$envCov),
                     model = 'RNMDs')

names(fullRN)

# Advanced options:
# Let's build again the W matrix
W.cov <- W.matrix(weather.data = maizeWTH,
                  by.interval = T,
                  statistic = 'quantile',
                  time.window = c(0, 14, 35, 60, 90, 120))

dim(W.cov)

# by using size_E = 'environment', get_kernel directly takes a W of q x q environments and builds a n x n matrix as EnvKernel()
EMM <- get_kernel(K_G = list(G = as.matrix(maizeG)),
                  Y = maizeYield,
                  K_E = list(W = H$envCov),
                  reaction = F,
                  model = 'EMM',
                  size_E = 'environment')

names(EMM)

# Its possible to integrate more than one environmental kernel
T.cov<- EnvTyping(weather.data=df.clim,var.id =  c('T2M','PRECTOT','WS2M'),env.id='env',format = 'wide')
eT <- EnvKernel(weather.data =T.cov,Y = maizeYield,merge = T,env.id = 'env',bydiag=TRUE)
dim(T.cov)
T.cov[, 1:4]



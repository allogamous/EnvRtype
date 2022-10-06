#'##############################################################################
#'
#' EnvRtype: tutorial code for the G3 manuscript
#' author: GCN, December 17th 2020
#' if you have any troubles, please contact <germano.cneto@gmail.com>
#' 
#'##############################################################################


################ BOX 1: Install EnvRtype ################

install.packages("devtools")
devtools::install_github("allogamous/EnvRtype")
require("EnvRtype")


################ BOX 2: Data sets ################
data("maizeYield") # toy set of phenotype data (grain yield per environment)
data("maizeG") # toy set of genomic relationship for additive effects 
data("maizeWTH") # toy set of environmental data
maizeYield <- droplevels(maizeYield) # make sure the nlevels of GID matchs with the G matrix dimension

################ BOX 3: Practical use of get_weather ################

env.data = get_weather(env.id = 'NAIROBI',country = 'KEN',
                       lat = -1.367,lon = 36.834,
                       start.day = '2015-03-01',end.day = '2015-04-01')

head(env.data)

################  BOX 4: Practical use of extract_GIS ################ 
data("clay_5_15")
env.data = extract_GIS(covraster = clay_5_15,name.out = 'clay_5_15',env.data = env.data)
head(env.data)

################ BOX 5: Practical use of SummaryWTH ################ 
summaryWTH(env.data = env.data, env.id = 'env', days.id = 'daysFromStart',statistic = 'mean')
summaryWTH(env.data = env.data) # by default

################ BOX 6: Practical use of param_temperature for Dry Beans in Nairobi, Kenya ################ 
TempData = param_temperature(env.data = env.data,Tbase1 = 8,Tbase2 = 45,Topt1 = 30,Topt2 = 35)
head(TempData)
env.data = param_temperature(env.data = env.data,Tbase1 = 8,Tbase2 = 45,Topt1 = 30,Topt2 = 35,merge = TRUE)
head(env.data) # merging TempData automatically

################ BOX 7: Practical use of param_atmospheric for Dry Bean Crop in Nairobi, Kenya ################ 
RadData = param_radiation(env.data = env.data) # first need to compute radiation parameters
head(RadData)
env.data = param_radiation(env.data = env.data,merge = TRUE)
AtmData  = param_atmospheric(env.data = env.data, Alt = 1628) 
head(AtmData)
env.data = param_atmospheric(env.data = env.data, Alt = 1628,merge = TRUE)
head(env.data)

################ BOX 8: Basic use of env_typing for typing temperature in Los Baños, Philipines, from 2000 to 2020 ################ 
env.data = get_weather(env.id = 'LOSBANOS',country = 'PHL',
                       lat = 14.170,lon = 121.241,variables.names = 'T2M',
                       start.day = '2000-03-01',end.day = '2020-03-01')

card = list(T2M=c(0,8,15,28,40,45,Inf)) # a list of vectors containing empirical and cardinal thresholds
env_typing(env.data = env.data,env.id = 'env', var.id = 'T2M', cardinals = card)

################ BOX 9: Basic use of env_typing for more than one variable ################
var = c("PRECTOT", "T2MDEW") # variables
env.data = get_weather(env.id = 'LOSBANOS',country = 'PHL',
                       lat = 14.170,lon = 121.241,variables.names = var,
                       start.day = '2000-03-01',end.day = '2020-03-01')
card = list(PRECTOT = c(0,5,10,25,40,100), T2MDEW = NULL) # cardinals and data-driven limits
env_typing(env.data = env.data,env.id = 'env', var.id = var, cardinals = card)
################ BOX 10: Basic use of env_typing for more than one variable ################
data("maizeWTH") # toy set of environmental data
var = c("PRECTOT", "T2MDEW", "T2M_MAX", "T2M_MIN") # variables
W = W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100,],
             var.id=var, statistic="mean", by.interval=TRUE)
dim(W)

################ BOX 11: Basic use of env_kernel ################
env_kernel(env.data = W, gaussian = FALSE)
env_kernel(env.data = W, gaussian = TRUE)

################ BOX 12: Basic usage of get_kernel function ################
data("maizeYield") # toy set of phenotype data (grain yield per environment)
data("maizeG"    ) # toy set of genomic relationship for additive effects 
data("maizeWTH")   # toy set of environmental data
y   = "value"      # name of the vector of phenotypes
gid = "gid"        # name of the vector of genotypes
env = "env"        # name of the vector of environments
maizeYield <- droplevels(maizeYield) # make sure the nlevels of GID matchs with the G matrix dimension

ECs  = W_matrix(env.data = maizeWTH, var.id = c("FRUE",'PETP',"SRAD","T2M_MAX"),statistic = 'mean')
## KG and KE might be a list of kernels
KE = list(W = env_kernel(env.data = ECs)[[2]])
KG = list(G=maizeG);
## Creating kernel models with get_kernel
MM    = get_kernel(K_G = KG, y = y,gid = gid,env = env, data = maizeYield,model = "MM") 
MDs   = get_kernel(K_G = KG, y = y,gid = gid,env = env,  data = maizeYield, model = "MDs") 
EMM   = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "EMM") 
EMDs  = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "EMDs") 
RMMM  = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMM") 
RNMDs = get_kernel(K_G = KG, K_E = KE, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs") 

################ BOX 13: Basic usage of kernel_model ################
fixed = model.matrix(~0+env, maizeYield)
MDs   = get_kernel(K_G = KG, y = y,gid = gid,env = env,  data = maizeYield, model = "MDs") 
fit   = kernel_model(y = y,env = env,gid = gid, data = maizeYield,random = MDs,fixed = fixed)

fit$yHat    # predicted phenotype values
fit$VarComp # variance components and confidence intervals
fit$BGGE    # full output of Hierarchical Bayesian Modeling

################ BOX 14: Remote Sensing for Several Places ################
env = c('GOI','TEX','BRI','MON','LOS','PON','CAL','PAL','DAV')
lat = c(-16.67,19.25,-27.47,43.61,14.170,6.294,3.261,-10.168,38.321)
lon = c(-49.25,-99.50,153.02,3.87,121.241,2.361,-76.312,-48.331,-121.442)
start = c('2020-03-15','2019-05-15','2018-09-15',
            '2017-06-18','2017-05-18','2016-07-18',
            '2017-11-18','2017-12-18','2018-07-18')
end = c('2020-04-15','2019-06-15','2018-10-15',
          '2017-07-18','2017-06-18','2016-08-18',
          '2017-12-18','2018-01-18','2018-08-18')
env.data = get_weather(env.id = env, lat = lat, lon = lon, start.day = start, end.day = end)


################ BOX 15: Discovering ETs and similarity among locations ################

ET = env_typing(env.data = env.data,env.id = 'env',var.id = 'T2M',format = 'wide')
EC = W_matrix(env.data = env.data,var.id = 'T2M')
distances = env_kernel(env.data = ET,gaussian = T)[[2]] 
kinship   = env_kernel(env.data = EC,gaussian = F, sd.tol = 3)[[2]]

ET = env_typing(env.data = env.data,env.id = 'env',var.id = 'T2M',format = 'wide')

# plot
require(ggplot2)
ggplot() + 
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(labels = scales::percent,expand = c(0,0))+
  geom_bar(data=ET, aes(y=Freq, x=env,fill=env.variable), 
           position = "fill",stat = "identity",colour='white',width = 1,size=.2)+
  # scale_y_continuous(labels = scales::percent,expand = c(0,0))+ #coord_flip()+
  ylab('Frequency of Occurence\n')+ 
  xlab("\nEnvironment")+
  labs(fill='Envirotype')+
  scale_fill_manual(values = c('blue4','green4',"orange",'red'))+
  theme(axis.title = element_text(size=19),
           axis.text = element_text(size=15),
        legend.text = element_text(size=12),
        strip.text = element_text(size=20),
        legend.title = element_text(size=17),
        strip.background = element_rect(fill="gray95",size=1),
        legend.position = 'right')

require(superheat)
require(viridis)
## other plot
superheat(distances,
        #  pretty.order.rows = TRUE,
       #   pretty.order.cols = TRUE,
       row.dendrogram = F,
       col.dendrogram = T,
          grid.vline.col = "white",
          grid.hline.col = "white",
          #row.dendrogram = T,
          legend.width = 4,
          left.label.size = 0.1,
          bottom.label.text.size = 5,
          bottom.label.size = 0.2,
          bottom.label.text.angle = 90,
          heat.pal = viridis::inferno(100),
          #heat.pal = viridis::magma(100),
          legend.text.size = 17,
          #   X.text = round(as.matrix(a),1),X.text.col="white",
          legend.height=0.2)

################  BOX 16: Envirotyping levels and model structures for GP with ECs ################ 
data("maizeYield") # toy set of phenotype data (grain yield per environment)
data("maizeG"    ) # toy set of genomic relationship for additive effects 
data("maizeWTH")   # toy set of environmental data
y   = "value"      # name of the vector of phenotypes
gid = "gid"        # name of the vector of genotypes
env = "env"        # name of the vector of environments
maizeYield <- droplevels(maizeYield) # make sure the nlevels of GID matchs with the G matrix dimension

### 1- Environmental Covariables (ECs)
stages    = c('VE','V1_V6','V6_VT','VT_R1','R1_R3','R3_R6',"H")
interval  = c(0,7,30,65,70,84,105)
EC1  = W_matrix(env.data = maizeWTH, var.id = 'FRUE')
EC2  = W_matrix(env.data = maizeWTH, var.id = 'PETP')
EC3  = W_matrix(env.data = maizeWTH, var.id = c('FRUE','PETP'))
EC4  = W_matrix(env.data = maizeWTH, var.id = 'FRUE',
                by.interval = T,time.window = interval,names.window = stages)
EC5  = W_matrix(env.data = maizeWTH, var.id = 'PETP',
                by.interval = T,time.window = interval,names.window = stages)
EC6  = W_matrix(env.data = maizeWTH, var.id = c('FRUE','PETP'),
                by.interval = T,time.window = interval,names.window = stages)

### 2- Kernels
K1 = list(FRUE = env_kernel(env.data = EC1)[[2]])
K2 = list(PETP = env_kernel(env.data = EC2)[[2]])
K3 = list(FRUE_PETP = env_kernel(env.data = EC3)[[2]])
K4 = list(FRUE = env_kernel(env.data = EC4)[[2]])
K5 = list(PETP = env_kernel(env.data = EC5)[[2]])
K6 = list(FRUE_PETP = env_kernel(env.data = EC6)[[2]])
### 3- Obtain Kernel Models
M0  = get_kernel(K_G = KG,  y = y,gid = gid,env = env,  data = maizeYield, model = "MDs") 
M1  = get_kernel(K_G = KG, K_E = K1, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs") 
M2  = get_kernel(K_G = KG, K_E = K2, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs") 
M3  = get_kernel(K_G = KG, K_E = K3, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs") 
M4  = get_kernel(K_G = KG, K_E = K4, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs") 
M5  = get_kernel(K_G = KG, K_E = K5, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs") 
M6  = get_kernel(K_G = KG, K_E = K6, y = y,gid = gid,env = env,  data = maizeYield, model = "RNMDs") 

## plot for each EC (e.g., EC1, EC2, EC3...)
superheat(env_kernel(env.data = EC1)[[2]],
            pretty.order.rows = TRUE,
             pretty.order.cols = TRUE,
          row.dendrogram = F,
          col.dendrogram = F,
          grid.vline.col = "white",
          grid.hline.col = "white",
          #row.dendrogram = T,
          legend.width = 4,
          left.label.size = 0.1,
          bottom.label.text.size = 5,
          bottom.label.size = 0.2,
          bottom.label.text.angle = 90,
          heat.pal = viridis::inferno(100),
        #  heat.pal = viridis::magma(100),
          legend.text.size = 17,
          #   X.text = round(as.matrix(a),1),X.text.col="white",
          legend.height=0.2)

################  BOX 17: Fitting Genomic-enabled models with enviromic data ################ 
fixed = model.matrix(~0+env,maizeWTH)

iter = 1000
burn = 500
seed = 78172
thin = 10
model = paste0('M',0:6)
Models = list(M0,M1,M2,M3,M4,M5,M6)
Vcomp <- c()
for(MODEL in 1:length(Models)){
  set.seed(seed)
  fit <- kernel_model(data = maizeYield,y = y,env = env,gid = gid,
                      random = Models[[MODEL]],fixed = Z_E,
                      iterations = iter,burnin = burn,thining = thin)
  
  Vcomp <- rbind(Vcomp,data.frame(fit$VarComp,Model=model[MODEL]))
}

Vcomp

reshape2::dcast(Vcomp,Model~Type,sum,value.var = 'Var')
reshape2::dcast(Vcomp,Model~Type,sum,value.var = 'CI_lower')
reshape2::dcast(Vcomp,Model~Type,sum,value.var = 'CI_upper')

################  BOX 18: Bulding enviromic kernels for each development stage ################ 
data("maizeYield") # toy set of phenotype data (grain yield per environment)
data("maizeG"    ) # toy set of genomic relationship for additive effects 
data("maizeWTH")   # toy set of environmental data
y   = "value"      # name of the vector of phenotypes
gid = "gid"        # name of the vector of genotypes
env = "env"        # name of the vector of environments

## Organizing Environmental Covariables (ECs) in W matrix
stages    = c('VE','V1_V6','V6_VT','VT_R1','R1_R3','R3_R6',"H")
interval = c(0,7,30,65,70,84,105)
id.vars  = names(maizeWTH)[c(10:15,23,25:30)]

W.matrix = W_matrix(env.data = maizeWTH,env.id = 'env',
                    var.id = id.vars,by.interval = T,time.window = interval,
                    names.window = stages,center = F,scale = F )


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
## Assembly Genomic and Enviromic Kernel Models
M1 = get_kernel(K_G = K_G,data = maizeYield,env = env,gid = gid,y = y, model = "MDs") # baseline model
M2 = get_kernel(K_G = K_G, K_E = K_F, data = maizeYield,env = env,gid = gid,
                y = y, model = "RNMM",dimension_KE = 'q') # reaction-norm 1
M3 = get_kernel(K_G = K_G, K_E = K_S, data = maizeYield,env = env,gid = gid,
                y = y,model = "RNMM",dimension_KE = 'q') # reaction-norm 2

model = c('Baseline Genomic','Reaction-Norm','Reaction-Norm for each Dev.Stage')
Models = list(M1,M2,M3) # for running all models in loop

iter = 10E3 # number of iterations
burn = 5E3  # number of burn in
thin = 10   # number for thining

Z_E = model.matrix(~0+env,data=maizeYield) # fixed environmental effects

Vcomp <- c()
for(MODEL in 1:length(Models)){
  set.seed(seed)
  fit <- kernel_model(data = maizeYield,y = y,env = env,gid = gid,
                      random = Models[[MODEL]],fixed = Z_E,
                      iterations = iter,burnin = burn,thining = thin)
  
  Vcomp <- rbind(Vcomp,data.frame(fit$VarComp,Model=model[MODEL]))
}

Vcomp



## plots for each kernel per stage (list K_S object)
superheat(K_S$R1_R3,
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          row.dendrogram = F,
          col.dendrogram = F,
          grid.vline.col = "white",
          grid.hline.col = "white",
          #row.dendrogram = T,
          legend.width = 4,
          left.label.size = 0.1,
          bottom.label.text.size = 5,
          bottom.label.size = 0.2,
          bottom.label.text.angle = 90,
          heat.pal = viridis::inferno(100),
          #  heat.pal = viridis::magma(100),
          legend.text.size = 17,
          #   X.text = round(as.matrix(a),1),X.text.col="white",
          legend.height=0.2)


################  BOX 19: Genomic Prediction using kernel_model ################ 
## Cross-validation to assess predictive ability of GP models (kernel_model function)
source('https://raw.githubusercontent.com/gcostaneto/SelectivePhenotyping/master/cvrandom.R')

rep  = 3
seed = 1010
f    = 0.20
iter = 5E3
burn = 1E3
thin = 10

## CV1
TS = Sampling.CV1(gids = maizeYield$gid,f = f,seed = seed,rep = rep,gidlevel = F)

require(foreach)
require(EnvRtype)

Y = maizeYield 
results <-foreach(REP = 1:rep, .combine = "rbind")%:%
  foreach(MODEL = 1:length(model), .combine = "rbind")%dopar% {
    
    
    yNA      <- Y
    tr       <- TS[[REP]]
    yNA$value[-tr] <- NA
    
    Z_E = model.matrix(~0+env,data=yNA) # fixed environmental effects
    
    fit <- kernel_model(data = yNA,y = y,env = env,gid = gid,
                        random = Models[[MODEL]],fixed = Z_E,
                        iterations = iter,burnin = burn,thining = thin)
    
    
    df<-data.frame(Model = model[MODEL],rep=REP,
                   rTr=cor(Y$value[tr ], fit$yHat[tr ],use = 'complete.obs'),
                   rTs=cor(Y$value[-tr], fit$yHat[-tr],use = 'complete.obs'))
    
    write.table(x = df,file = 'point-estimate_r.txt',sep=',',append = T,row.names=T)
    
    output <- data.frame(obs=Y$value,pred=fit$yHat,
                         gid=Y$gid, env=Y$env,
                         Model = model[MODEL],rep=REP,pop=NA)
    
    
    output$pop[tr ] <- 'training'
    output$pop[-tr] <- 'testing'
    
    return(output)
  }
#stopCluster(cl)

results
require(plyr)
pa = ddply(results,.(rep,pop,Model),summarise,r = cor(obs,pred))
ddply(pa,.(pop,Model),summarise, pa = round(mean(r),3),sd = round(sd(r),3))


## CV00 (sampling genotypes and environments)

seed = 8172
# out environments (as testing set) = 2
TS=Sampling.CV0(gids = Y$gid,envs = Y$env,out.env = 2,f = f,seed = seed,rep = rep)

cl <- makeCluster(3)
registerDoParallel(cl)

results <-foreach(REP = 1:rep, .combine = "rbind")%:%
  foreach(MODEL = 1:length(model), .combine = "rbind")%dopar% {
    
    
    yNA      <- Y
    tr       <- TS[[REP]]$training
    yNA$value[-tr] <- NA
    
    Z_E = model.matrix(~0+env,data=Y)
    fit <- kernel_model(data = yNA,y = y,env = env,gid = gid,
                        random = Models[[MODEL]],fixed = Z_E,
                        iterations = iter,burnin = burn,thining = thin)
    
    
    output <- data.frame(obs=Y$value,pred=fit$yHat,
                         gid=Y$gid, env=Y$env,
                         Model = model[MODEL],rep=REP,pop=NA)
    
      
    output$pop[tr ] <- 'training'
    output$pop[-tr] <- 'testing'
    
    return(output)
  }

stopCluster(cl)

pa = ddply(results,.(rep,pop,Model),summarise,r = cor(obs,pred))
ddply(pa,.(pop,Model),summarise, pa = round(mean(r),3),sd = round(sd(r),3))

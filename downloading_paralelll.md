# **Download weather and elevation data into a faster way **

*last update: June 21th 2021 by Germano Costa Neto (germano.cneto@gmail.com)*

- This tutorial aims to provide tools for speeding up the download of a wide number of environments using get_wether function.

# Contents

* [Case 1: More than 50 environments](#P1)
* [Case 2: Diverse locations across diverse countries arround the world](#P2)
* [Case 3 : A wide number of environments for a same given location](#P3)
* [Case 4: Using WorldClim data base](#P4)

<div id="P1" />

# Case 1 : More than 50 environments

```{r}
# coordinates
id_info <- read.csv(file = 'https://raw.githubusercontent.com/allogamous/EnvRtype/master/Supplementary%20Source%20and%20Data/Brazil_city.csv')
# step 1: check coordinates
# step 2: check identifications of city
str(id_info)


# step 3: check dates
id_info$start
levels(id_info$start)
id_info    =   droplevels(id_info[!id_info$start %in% levels(id_info$start)[1],])

# step 4: translate dates using as.Date
id_info$start = as.Date(id_info$start)
id_info$end = id_info$start+120

# step 4: create an environment identification
id_info$environment = paste(id_info$city,id_info$season,id_info$month,sep = '_')

#id_info <- id_info[1:10,]

## download covariables
require("EnvRtype")
require("foreach") # install foreach
require('doParallel') # install doParallel

# running in parallel is faster!
nclust = 3
cl <- makeCluster(nclust,outfile="environments.txt")
registerDoParallel(cl)

climate = foreach(E = 1:.Ne,.combine="rbind") %dopar% 
  {
    cat(id_info$environment[E])
    climate = EnvRtype::get_weather(env.id = id_info$environment[E],lat = id_info$Latitude[E],
                                    lon = id_info$Longitude[E],start.day = id_info$start[E],
                                    end.day =  id_info$end[E],country = 'BRA')
    return(climate)
  }

stopImplicitCluster()
stopCluster(cl)
```

<div id="P2" />

# Case 2 : different locations across diverse countries arround the world


<div id="P3" />

# Case 3 : A wide number of environments for a same given location


<div id="P4" />

# Case 4: Using WorldClim data base




              



# Software

```{r, eval=FALSE}
library(devtools)
install_github('allogamous/EnvRtype')
library(EnvRtype)
```

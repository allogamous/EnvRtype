#' Maize Phenotypic Data Set
#'
#'This data set was included in Souza et al. (2017) and Cuevas et al (2019) is from
#'the Helix Seeds Company (HEL). It consists of grain yied from 452 maize hybrids
#'obtained by crossing 111 pure lines (inbreds); the hybrids were evaluated in 2015 at five Brazilian sites (E1-E5). The experimental design
#'used in each site was a randomized block with two replicates per hybrid. However, to facilitate the demonstration of functions,
#'only 150 hybrids per environment are being considered, thus counting 750 genotype x environment observations.
#'Grain yield data are mean-centered and scaled.
#'
#'@docType data
#'
#'@usage data(maizeYield)
#'
#'@examples
#' data(maizeYield); head(maizeYield)
#'
#' @format A data frame with 750 rows and 3 variables:
#' \describe{
#'   \item{env}{environmental id (factor)}
#'   \item{gid}{genotypic id (factor)}
#'   \item{value}{ grain yield value of genotype plus genotype by location interaction effects in kg ha-1 (numeric)}
#'   ...
#' }
"maizeYield"

#save(maizeWTH,file = 'maizeWTH.RData')
#save(maizeYield,file = 'maizeYield.RData')
#save(maizeG,file = 'maizeG.RData')

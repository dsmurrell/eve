library(devtools)
setwd("~/Dropbox/projects/eve/roxygen")
document('../eve')
setwd("~/Dropbox/projects/eve/examples/QSPR/LogS/")

source("http://bioconductor.org/biocLite.R")
biocLite("multicore")

library(eve)

ls("package:eve")
?BuildCaretEVEstimator
?BuildEVEstimator
?CEC
?PredictSigmas

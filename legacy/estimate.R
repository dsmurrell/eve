# creates f distance matricies where f is the numer of folds
# each distance matrix is composed of the Euclidean distance of each instance in a fold to all the instances in the other folds
# this is the non-multicore version of this function and is no longer used
CreateDistanceMatrices <- function(x.train, folds) {
  dms <- list()
  for(fold in 1:length(folds)) {
    out <- as.matrix(x.train[folds[[fold]], ])
    ref <- as.matrix(x.train[-folds[[fold]], ])
    dm <- matrix(NA, nrow=nrow(ref), ncol=nrow(out))
    for(i in 1:nrow(out)) {
      y <- as.vector(out[i,], mode="numeric")
      y[is.na(y)] <- 0
      dm[,i] <- colSums((t(ref)-y)^2)
    }
    dms[[fold]] <- dm
  }
  dms
}
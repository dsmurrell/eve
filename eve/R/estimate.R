# creates a distance matrix from one fold of the data to all the other folds
CreateDistanceFold <- function(fold, x.train, folds) {
  out <- as.matrix(x.train[folds[[fold]], ])
  ref <- as.matrix(x.train[-folds[[fold]], ])
  dm <- matrix(NA, nrow=nrow(ref), ncol=nrow(out))
  for(i in 1:nrow(out)) {
    y <- as.vector(out[i,], mode="numeric")
    y[is.na(y)] <- 0
    dm[,i] <- colSums((t(ref)-y)^2)
  }
  dm
}

# this is not used at the moment (was tested once and was found to work better for kernel based models such as the SVM, not tested on the RVM or GP yet)
# creates a distance matrix from one fold of the data to all the other folds (this version used the the distance in a kernel space)
# library(kernlab)
CreateKernelDistanceFold <- function(fold, x.train, folds) {
  k <- rbfdot(sigma=0.00391)
  out <- as.matrix(x.train[folds[[fold]], ])
  ref <- as.matrix(x.train[-folds[[fold]], ])
  dm <- matrix(NA, nrow=nrow(ref), ncol=nrow(out))
  for(i in 1:nrow(out)) {
    y <- as.vector(out[i,], mode="numeric")
    y[is.na(y)] <- 0 
    dm[,i] <- 1-apply(ref, 1, k, y)
  }
  dm
}
  
# uses CreateDistanceFold to create a list of distance matricies from each fold to the rest of the data points
# this can be carried out in parralel
CreateDistanceMatricesMC <- function(x.train, folds, cores = 1) {
  mclapply(1:length(folds), CreateDistanceFold, x.train, folds, mc.cores=max(length(folds), cores))
}

# confine error-variance score 
confine <- function(col, dm, errors, obs, preds, N) {
  x <- dm[,col]
  sum((errors[order(x)[1:N]])^2)/N
}

# convive error-variance score
convive <- function(col, dm, errors, obs, preds, N) {
  if(N==1) return(0)
  x <- dm[,col]
  var(obs[order(x)[1:N]])
}

# confuse error-variance score
confuse <- function(col, dm, errors, obs, preds, N) {
  if(N==1) return(0)
  x <- dm[,col]
  var(errors[order(x)[1:N]])
}

# distance error-variance score
distance <- function(col, dm, errors, obs, preds, N) {
  x <- dm[,col]
  mean(sqrt(x[order(x)][1:N]))
}

# diffNN error-variance score
diffNN <- function(col, dm, errors, obs, preds, N) {
  x <- dm[,col]
  abs(mean(obs[order(x)[1:N]]) - preds[col])
}

# applies one of the estimators, specified by func, to the data points in a particular fold
DoFold <- function(fold, folds, dms, errors, obs, preds, N, func) {
  sapply(1:ncol(dms[[fold]]), func, dms[[fold]], abs(errors)[-folds[[fold]]], obs[-folds[[fold]]], preds[folds[[fold]]], N)
}

# applies one of the estimators to all the data points by iterating over the all the folds
CalculateErrorEstimates <- function(folds, dms, errors, obs, preds, N, func) {
  num.folds <- length(folds)
  estimates <- lapply(1:num.folds, DoFold, folds, dms, errors, obs, preds, N, func)
  estimates.vec <- do.call(c, estimates)
  fold.vec <- do.call(c, folds)
  match <- match(1:length(fold.vec), fold.vec) # reordering is performed here to get the result into the original folds indexes order
  estimates.vec[match]
}

# for a particular number of neighbours create a matrix which whose columns relate to an implemented estimator
CombineErrorEstimates <- function(N, folds, dms, errors, obs, preds) {
  m <- CalculateErrorEstimates(folds, dms, errors, obs, preds, N, confine)
  m <- cbind(m, CalculateErrorEstimates(folds, dms, errors, obs, preds, N, confuse))
  m <- cbind(m, CalculateErrorEstimates(folds, dms, errors, obs, preds, N, convive))
  m <- cbind(m, CalculateErrorEstimates(folds, dms, errors, obs, preds, N, distance))
  m <- cbind(m, CalculateErrorEstimates(folds, dms, errors, obs, preds, N, diffNN))
  m
}

# combine all the matrices above (cbind) from 1 neighbour to the maximum number of neighbours considered (Nmax)
# not used at the moment because not multicore
CreateEstimatorMatrix <- function(Nmax, folds, dms, errors, obs, preds) {
  m <- do.call(cbind, lapply(1:Nmax, CombineErrorEstimates, folds, dms, errors, obs, preds))
  # remove the columns where all the results are 0 (this is columns 2 and 3 for the moment)
  # m <- m[,-which(apply(m, 2, sum)==0)]
  m <- m[,-c(2,3)]
}

# combine all the matrices above (cbind) from 1 neighbour to the maximum number of neighbours considered (Nmax)
# this is a multicore implementation of the above function
CreateEstimatorMatrixMC <- function(Nmax, folds, dms, errors, obs, preds, cores = 1) {
  m <- do.call(cbind, mclapply(1:Nmax, CombineErrorEstimates, folds, dms, errors, obs, preds, mc.cores=cores))
  # remove the columns where all the results are 0 (this is columns 2 and 3 for the moment)
  # m <- m[,-which(apply(m, 2, sum)==0)]
  m <- m[,-c(2,3)]
}

# shortcut for finding maximum likelihood of the variance assuming 0 mean
FitSigma <- function(x) {
  sum(x^2)/length(x)
}

# for each estimator function at a specific number of considered neighbours (1 column in the estimator matrix)
# find two splines
# each data point has an estimator value and an absolute error
# order the data points by estimator value
# define a window width at 1/5 the number of data points as are in the training set
# half that with (hwidth) is 1/10 the number of data point as are in the training set
# for each data point inside the whole window use the absolute errors to estimate the variance for that index
# make a mapping from an estimator value to an index into the sorted estimators
# make a mapping from the index to the window sigma
# now, when we obtain an estimator value for a new data point, the index can be found 
# using the first mapping and and consequently the variance estimate of the window using the second mapping

# indexes are used instead of using a fixed size window over the estimator space because
# the distribution of the space could vary and outliers might affect the window construction algorithm
# this method was proven to work better than using the range of the dataset and dividing up that space into windows

ConvertToSigmaSplines <- function(x, errors, spar = 0.8, p = FALSE) {
  splines <- list()
  order <- order(x, decreasing=TRUE)
  error.sort <- errors[order]
  abs.error.sort <- abs(error.sort)
  x.sort <- x[order]
  hwidth <- floor(length(x)/10)
  # initialise a 'interval sigma estimate'
  ise <- rep(NA, length(x))
  inside <- hwidth:(length(x)-hwidth)
  for(i in inside) {
    set <- abs.error.sort[(i-hwidth):(i+hwidth)]
    ise[i] <- FitSigma(abs.error.sort[(i-hwidth):(i+hwidth)])
  }
  
  if(p) {
    print(plot(x.sort, 1:length(x), xlab="Sorted estimator values", ylab="Index"))
    print(plot(inside, ise[inside], xlab="Index", ylab="Interval Sigma Estimate"))
  }
  
  # index.spline creates a mapping from the estimator value to the index into the sorted estimators
  splines$index.spline <- smooth.spline(x.sort, 1:length(x), spar = spar)
  # ise.spline creates a mapping from the index into the sorted estimators to the interval sigma
  splines$ise.spline <- smooth.spline(inside, ise[inside], spar = spar)
  # using these two splines, one can obtain the interval sigma estimate from the estimator value for a new data point
  return(splines)
}

# for each column of the estimator matrix get the two sigma splines which create a mapping from the estimator of a new datapoint to the interval sigma
# return a list of these spline pairs (the pairs are a sublist)
GetSigmaSplineList <- function(m, errors) {
  sl <- list()
  for(i in 1:ncol(m)) {
    sl[[i]] <- ConvertToSigmaSplines(m[,i], errors, i==4)
  }
  sl # return the spline list
} 

# find the sigma estimate from the two sigma splines given the original estimator values 
GetSigmasFromSplines <- function(i, m, sl) {
  indexes <- predict(sl[[i]]$index.spline, m[,i])$y
  predict(sl[[i]]$ise.spline, indexes)$y
}

#' Confidence–error correlation is a measure of the performance of the error estimates against the absolute errors
#' 
#' CEC uses the Pearson's correlation between the sigma estimates of the assumed normal error distributions 
#' and the absolute errors. This is normalised by the perfect confidence estimator
#' defined as the Pearson's correlation between the sorted sigmas and the sorted absolute errors.
#' CEC = cor(sigmas, |errors|) / cor(sort(sigmas), sort(|errors|))
#' 
#' @param x First variable
#' @param y Second variable
#' @export
#' @return CEC value
#' @references No Longer Confidential: Estimating the Confidence of Individual Regression Predictions:  
#' \url{http://dx.doi.org/10.1371/journal.pone.0048723}
#' @author Daniel Murrell <dsmurrell@@gmail.com>
CEC <- function(x, y) {
  if(sum(x)==0) return(0)
  cor(x,y)/cor(sort(x), sort(y))
}

# Greedy ensemble algorithm
GreedOpt <- function(X, Y, optFunc = CEC, iter = 100L) {
  bests       <- c()
  N           <- ncol(X)
  weights     <- rep(0L, N)
  pred        <- 0 * X
  sum.weights <- 0
  absY <- abs(Y)
  
  while(sum.weights < iter) {
    sum.weights   <- sum.weights + 1L
    pred          <- (pred + X) * (1L / sum.weights)
    CECs          <- apply(pred, 2, optFunc, absY)
    best          <- which.max(CECs)
    bests         <- c(bests, max(CECs))
    weights[best] <- weights[best] + 1L
    pred          <- pred[, best] * sum.weights
  } 
  l <- list()
  l$weights <- weights
  l$bests <- bests
  return(l)
}

# Application of the greedy ensemble to select sigma matrix weights that maximises the CEC function 
GetWeights <- function(m, errors, sl, cores = 1, optFunc = CEC) {
  sigma.matrix <- do.call(cbind, mclapply(1:ncol(m), GetSigmasFromSplines, m, sl, mc.cores=cores))
  l <- GreedOpt(sigma.matrix, errors, optFunc, iter=50)
  weights <- l$weights
  #plot(l$bests)
  weights/sum(weights)
}  

#' Builds an error variance estimator using training set, fold indexes, observed values and CV predictions
#' 
#' An error variance estimator is built using the training set, the fold indexes used during training, 
#' the observed values of the training set and the cross-validation predictions during training. 
#' If the cross-validation predictions are kept during training, this should allow \code{eve} to
#' be backwards compatible to any trained model. The estimator can be applied to make error estimates for
#' new data points using the PredictSigmas function.
#' 
#' @param x Training set
#' @param folds Folds used in cross validation
#' @param obs Observed values
#' @param preds Crossvalidated predicitions
#' @param Nmax Maximum number of neighbours to consider
#' @param cores Number of cores to utilise during estimator building
#' @param optFunc Optimisation function used to score the correlation
#' @export
#' @return Trained error variance estimator
#' @author Daniel Murrell <dsmurrell@@gmail.com>
BuildEVEstimator <- function(x, folds, obs, preds, Nmax = 20, cores = 1, optFunc = CEC) {
  errors <- preds - obs
  dms <- CreateDistanceMatricesMC(x, folds, cores)
  m <- CreateEstimatorMatrixMC(Nmax, folds, dms, errors, obs, preds, cores)
  sl <- GetSigmaSplineList(m, errors)
  weights <- GetWeights(m, errors, sl, cores, optFunc)
  list(x = x, sl = sl, weights = weights, obs = obs, preds = preds, errors = errors, Nmax = Nmax)
}

#' Builds an error variance estimator using a caret model as input (caret model must have saved CV predictions)
#' 
#' An error variance estimator is built using a caret model that was trained using cross-validation.
#' The \code{savePredictions} argument to the \code{trainControl} function needs to be set to TRUE
#' so that caret saves the cross-validation predictions during model training.
#' The estimator can be applied to make error estimates for
#' new data points using the PredictSigmas function.
#' 
#' @param x Training set
#' @param model Trained caret model (must have saved CV predictions)
#' @param Nmax Maximum number of neighbours to consider
#' @param cores Number of cores to utilise during estimator building
#' @param optFunc Optimisation function used to score the correlation
#' @export
#' @return Trained error variance estimator
#' @author Daniel Murrell <dsmurrell@@gmail.com>
BuildCaretEVEstimator <- function(x, model, Nmax = 20, cores=1, optFunc = CEC) {
  indexes <- which(apply(as.data.frame(model$pred[,names(model$bestTune)]), 1, function(x,y) {identical(as.numeric(x), as.numeric(y))}, model$bestTune))
  best <- model$pred[indexes,]
  order <- order(best$rowIndex, decreasing=FALSE)
  best <- best[order,]
  obs <- best$obs
  preds <- best$pred
  fold.strings <- unique(best$Resample)
  folds <- list()
  for(fold in 1:length(fold.strings)) {
    folds[[fold]] <- which(best$Resample == fold.strings[fold])
  }
  BuildEVEstimator(x, folds, obs, preds, Nmax, cores, optFunc)
}
 
#-----------------------------
# Estimates on new predictions
#-----------------------------

# distances matrix from the new set of datapoints to the rest of the training set
CreateNewDistanceMatrix <- function(x, x.new) {
  out <- x.new
  ref <- as.matrix(x)
  dm <- matrix(NA, nrow=nrow(ref), ncol=nrow(out))
  for(i in 1:nrow(out)) {
    y <- as.vector(out[i,], mode="numeric")
    y[is.na(y)] <- 0
    dm[,i] <- colSums((t(ref)-y)^2)
  }
  #max <- max(dm)  # not sure if this is necessary
  #for(i in 1:nrow(dm)) {
  #  dm[i,i] <- max
  #}
  dm
}

# kernel version of CreateNewDistanceMatrix, not used yet
CreateNewKernelDistanceMatrix <- function(x, x.new) {
  k <- rbfdot(sigma=0.00391)
  out <- x.new
  ref <- as.matrix(x)
  dm <- matrix(NA, nrow=nrow(ref), ncol=nrow(out))
  for(i in 1:nrow(out)) {
    y <- as.vector(out[i,], mode="numeric")
    y[is.na(y)] <- 0
    dm[,i] <- 1-apply(t(ref), 2, k, y)
  }
  #max <- max(dm) # not sure if this is necessary
  #for(i in 1:nrow(dm)) {
  #  dm[i,i] <- max
  #}
  dm
}

DoNew <- function(dm, errors, obs, preds, N, func) {
  sapply(1:ncol(dm), func, dm, abs(errors), obs, preds, N)
}

CombineNewErrorEstimates <- function(N, dm, errors, obs, preds) {
  m <- DoNew(dm, errors, obs, preds, N, confine)
  m <- cbind(m, DoNew(dm, errors, obs, preds, N, confuse))
  m <- cbind(m, DoNew(dm, errors, obs, preds, N, convive))
  m <- cbind(m, DoNew(dm, errors, obs, preds, N, distance))
  m <- cbind(m, DoNew(dm, errors, obs, preds, N, diffNN))
  m
}

CreateNewEstimatorMatrix <- function(Nmax, dm, errors, obs, preds) {
  m <- do.call(cbind, lapply(1:Nmax, CombineNewErrorEstimates, dm, errors, obs, preds))
  m[,-c(2,3)]
}

GetNewSigmaMatrix <- function(nm, sl, cores=1) {
  do.call(cbind, mclapply(1:ncol(nm), GetSigmasFromSplines, nm, sl, mc.cores=cores))
}

GetNewSigmas <- function(sigma.matrix, weights) {
  apply(sweep(sigma.matrix,2,weights,`*`), 1, sum)
}

#' Predict the error variance of new data points given the dataset and a trained estimator
#' 
#' Given new datapoints, \code{PredictSigmas} uses a trained error variance estimator to predict the
#' sigmas of the error distribution for those datapoints. It is important that the descriptors are 
#' transformed in exactly the same way as the descriptors were transformed before model training.
#' 
#' @param x Descriptors of the new datapoints transformed in exactly the same way as they were during training
#' @param estimator Trained estimator returned by either \code{\link{BuildEVEstimator}} or \code{\link{BuildCaretEVEstimator}}
#' @export
#' @return Predicted error variances for the new data points
#' @author Daniel Murrell <dsmurrell@@gmail.com>
PredictSigmas <- function(x, estimator) {
  dm <- CreateNewDistanceMatrix(estimator$x, x)
  nm <- CreateNewEstimatorMatrix(estimator$Nmax, dm, estimator$errors, estimator$obs, estimator$preds)
  sigma.matrix <- GetNewSigmaMatrix(nm, estimator$sl)
  sigmas <- GetNewSigmas(sigma.matrix, estimator$weights)
  prediction <- list()
  prediction$dm <- dm
  prediction$nm <- nm
  prediction$sigma.matrix <- sigma.matrix
  prediction$sigmas <- sigmas
  prediction
}



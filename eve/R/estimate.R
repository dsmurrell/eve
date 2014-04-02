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
# sort the estimator by the 
# find two splines
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
  sl
} 

# find the sigma estimate from the two sigma splines given the original estimator values 
GetSigmasFromSplines <- function(i, m, sl) {
  indexes <- predict(sl[[i]]$index.spline, m[,i])$y
  predict(sl[[i]]$ise.spline, indexes)$y
}

# some measure of the correlation between two variables
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

# Application of the greed ensemble to select sigma matrix weights that moximise the CEC function 
GetWeights <- function(m, errors, sl, cores = 1, optFunc = CEC) {
  sigma.matrix <- do.call(cbind, mclapply(1:ncol(m), GetSigmasFromSplines, m, sl, mc.cores=cores))
  l <- GreedOpt(sigma.matrix, errors, optFunc, iter=50)
  weights <- l$weights
  plot(l$bests)
  weights/sum(weights)
}  

# builds an error variance estimator using training set, fold indexes, observed values and CV predictions
BuildEVEstimator <- function(x, folds, obs, preds, Nmax = 20, cores = 1, optFunc = CEC) {
  errors <- preds - obs
  dms <- CreateDistanceMatricesMC(x, folds, cores)
  m <- CreateEstimatorMatrixMC(Nmax, folds, dms, errors, obs, preds, cores)
  sl <- GetSigmaSplineList(m, errors)
  weights <- GetWeights(m, errors, sl, cores, optFunc)
  list(x = x, sl = sl, weights = weights, obs = obs, preds = preds, errors = errors, Nmax = Nmax)
}

# builds an error variance estimator using a caret model as input (caret model must have saved CV predictions)
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

GetNewSigmas <- function(nm, sl, weights, cores=1) {
  sigma.matrix <- do.call(cbind, mclapply(1:ncol(nm), GetSigmasFromSplines, nm, sl, mc.cores=cores))
  apply(sweep(sigma.matrix,2,weights,`*`), 1, sum)
}

GetNewSigmaMatrix <- function(nm, sl, cores=1) {
  do.call(cbind, mclapply(1:ncol(nm), GetSigmasFromSplines, nm, sl, mc.cores=cores))
}

GetNewSigmas <- function(sigma.matrix, weights) {
  apply(sweep(sigma.matrix,2,weights,`*`), 1, sum)
}

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



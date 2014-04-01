# creates f distance matricies where f is the numer of folds
# each distance matrix is composed of the Euclidean distance of each instance in a fold to all the instances in the other folds
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

#library(kernlab)

# creates a distance matrix from one fold of the data to all the other folds
# this is so multiple cores can be used in the function below
CreateDistanceFold <- function(fold, x.train, folds) {
  k <- rbfdot(sigma=0.00391)
  out <- as.matrix(x.train[folds[[fold]], ])
  ref <- as.matrix(x.train[-folds[[fold]], ])
  dm <- matrix(NA, nrow=nrow(ref), ncol=nrow(out))
  for(i in 1:nrow(out)) {
    y <- as.vector(out[i,], mode="numeric")
    y[is.na(y)] <- 0
    dm[,i] <- 1-apply(ref, 1, k, y)
    #dm[,i] <- colSums((t(ref)-y)^2)
  }
  dm
}

#m <- matrix(rnorm(20), 4, 5)
#v <- rnorm(5)
#k <- rbfdot(sigma=0.00391)
#1 - apply(m, 1, k, v)
  
# uses the above function to create a list of distance matricies
CreateDistanceMatricesMC <- function(x.train, folds, cores = 1) {
  mclapply(1:length(folds), CreateDistanceFold, x.train, folds, mc.cores=max(length(folds), cores))
}

# calculate the confine error-variance score for the instance in the col column of dm
# given the 
confine <- function(col, dm, errors, obs, preds, N) {
  x <- dm[,col]
  sum((errors[order(x)[1:N]])^2)/N
}

confuse <- function(col, dm, errors, obs, preds, N) {
  if(N==1) return(0)
  x <- dm[,col]
  var(errors[order(x)[1:N]])
}

convive <- function(col, dm, errors, obs, preds, N) {
  if(N==1) return(0)
  x <- dm[,col]
  var(obs[order(x)[1:N]])
}

distance <- function(col, dm, errors, obs, preds, N) {
  x <- dm[,col]
  mean(sqrt(x[order(x)][1:N]))
}

diffNN <- function(col, dm, errors, obs, preds, N) {
  x <- dm[,col]
  abs(mean(obs[order(x)[1:N]]) - preds[col])
}

DoFold <- function(fold, folds, dms, errors, obs, preds, N, func) {
  sapply(1:ncol(dms[[fold]]), func, dms[[fold]], abs(errors)[-folds[[fold]]], obs[-folds[[fold]]], preds[folds[[fold]]], N)
}

CalculateErrorEstimates <- function(folds, dms, errors, obs, preds, N, func) {
  num.folds <- length(folds)
  estimates <- lapply(1:num.folds, DoFold, folds, dms, errors, obs, preds, N, func)
  estimates.vec <- do.call(c, estimates)
  fold.vec <- do.call(c, folds)
  match <- match(1:length(fold.vec), fold.vec)
  estimates.vec[match]
}

CombineErrorEstimates <- function(N, folds, dms, errors, obs, preds) {
  m <- CalculateErrorEstimates(folds, dms, errors, obs, preds, N, confine)
  m <- cbind(m, CalculateErrorEstimates(folds, dms, errors, obs, preds, N, confuse))
  m <- cbind(m, CalculateErrorEstimates(folds, dms, errors, obs, preds, N, convive))
  m <- cbind(m, CalculateErrorEstimates(folds, dms, errors, obs, preds, N, distance))
  m <- cbind(m, CalculateErrorEstimates(folds, dms, errors, obs, preds, N, diffNN))
  m
}

CreateEstimatorMatrix <- function(Nmax, folds, dms, errors, obs, preds) {
  m <- do.call(cbind, lapply(1:Nmax, CombineErrorEstimates, folds, dms, errors, obs, preds))
  # remove the columns where all the results are 0 (this is columns 2 and 3 for the moment)
  # m <- m[,-which(apply(m, 2, sum)==0)]
  m <- m[,-c(2,3)]
}

CreateEstimatorMatrixMC <- function(Nmax, folds, dms, errors, obs, preds, cores = 1) {
  m <- do.call(cbind, mclapply(1:Nmax, CombineErrorEstimates, folds, dms, errors, obs, preds, mc.cores=cores))
  # remove the columns where all the results are 0 (this is columns 2 and 3 for the moment)
  # m <- m[,-which(apply(m, 2, sum)==0)]
  m <- m[,-c(2,3)]
}

# shortcut for finding maximum likelihood assuming 0 mean
FitSigma <- function(x) {
  sum(x^2)/length(x)
}

ConvertToSigmaSplines <- function(x, errors, spar = 0.8) {
  splines <- list()
  order <- order(x, decreasing=TRUE)
  error.sort <- errors[order]
  abs.error.sort <- abs(error.sort)
  x.sort <- x[order]
  hwidth <- floor(length(x)/10)
  cis <- rep(NA, length(x))
  inside <- hwidth:(length(x)-hwidth)
  for(i in inside) {
    set <- abs.error.sort[(i-hwidth):(i+hwidth)]
    cis[i] <- FitSigma(abs.error.sort[(i-hwidth):(i+hwidth)])
  }
  
  splines$index.spline <- smooth.spline(x.sort, 1:length(x), spar = spar)
  splines$ci.spline <- smooth.spline(inside, cis[inside], spar = spar)
  return(splines)
}

# tbd, convert to multiple core implementation
GetSigmaSplineList <- function(m, errors) {
  sl <- list()
  for(i in 1:ncol(m)) {
    sl[[i]] <- ConvertToSigmaSplines(m[,i], errors)
  }
  sl
} 

GetSigmasFromSplines <- function(i, m, sl) {
  indexes <- predict(sl[[i]]$index.spline, m[,i])$y
  predict(sl[[i]]$ci.spline, indexes)$y
}

CEC <- function(x, y) {
  if(sum(x)==0) return(0)
  cor(x,y)/cor(sort(x), sort(y))
}

GreedOptCEC <- function(X, Y, optFunc = CEC, iter = 100L){
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

GetWeights <- function(m, errors, sl, cores = 1, optFunc = CEC) {
  sigma.matrix <- do.call(cbind, mclapply(1:ncol(m), GetSigmasFromSplines, m, sl, mc.cores=cores))
  l <- GreedOptCEC(sigma.matrix, errors, optFunc, iter=50)
  weights <- l$weights
  plot(l$bests)
  weights/sum(weights)
}  

BuildEstimator <- function(x, folds, obs, preds, Nmax = 20, cores = 1, optFunc = CEC) {
  errors <- preds - obs
  dms <- CreateDistanceMatricesMC(x, folds, cores)
  m <- CreateEstimatorMatrixMC(Nmax, folds, dms, errors, obs, preds, cores)
  sl <- GetSigmaSplineList(m, errors)
  weights <- GetWeights(m, errors, sl, cores, optFunc)
  list(x = x, sl = sl, weights = weights, obs = obs, preds = preds, errors = errors, Nmax = Nmax)
}

# builds an error estimator using a caret model as input
BuildCaretErrorEstimator <- function(x, model, Nmax = 20, cores=1, optFunc = CEC) {
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
  BuildEstimator(x, folds, obs, preds, Nmax, cores, optFunc)
}
 
# estimates on new predictions
CreateNewDistanceMatrix <- function(x, x.new) {
  k <- rbfdot(sigma=0.00391)
  out <- x.new
  ref <- as.matrix(x)
  dm <- matrix(NA, nrow=nrow(ref), ncol=nrow(out))
  for(i in 1:nrow(out)) {
    y <- as.vector(out[i,], mode="numeric")
    y[is.na(y)] <- 0
    dm[,i] <- 1-apply(t(ref), 2, k, y)
    #dm[,i] <- colSums((t(ref)-y)^2)
  }
  #max <- max(dm)
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



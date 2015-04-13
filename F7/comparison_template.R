#!/usr/bin/Rscript

args<-commandArgs(TRUE)
name_dataset <- args[1]
print(name_dataset)
direct <- paste("/Users/icortes/Desktop/eve/",name_dataset,sep="")
#####!/c5/shared/R/3.0.3/bin/Rscript

# load required packages
#library(eve)
library(camb)
library(doMC)
filename <- paste(name_dataset,".smi",sep="")

StandardiseMolecules(structures.file=filename,
standardised.file="standardised.sdf",
removed.file="removed.sdf",
#output="standardisation_info.csv",
remove.inorganic=TRUE,
fluorine.limit=-1,
chlorine.limit=-1,
bromine.limit=-1,
iodine.limit=-1,
min.mass.limit=-1, #suggested value 20 
max.mass.limit=-1) #suggested 

properties <- read.table("properties.csv", header=TRUE, sep="\t") 
properties <- properties[properties$Kept==1, ]
head(properties)
targets <- read.table(paste(name_dataset,".bio",sep=""),header=T)$pIC50
print(length(targets))

descriptor.types <- c("2D")
descriptors <- GeneratePadelDescriptors(standardised.file = "standardised.sdf",
    types = descriptor.types, threads = 1)
descriptors <- RemoveStandardisedPrefix(descriptors)

descriptors <- ReplaceInfinitesWithNA(descriptors)
descriptors <-  ImputeFeatures(descriptors)
saveRDS(descriptors, file = "descriptors.rds")

descriptors <- readRDS("descriptors.rds")
print(names(descriptors[1]))
tot

# the number of folds used during the error estimation cross validation
NF <- 30

dataset <- SplitSet(ids=descriptors$Name, x=descriptors[,2:ncol(descriptors)], y=targets, percentage=0.1)
dataset <- GetCVTrainControl(dataset, folds=NF, savePredictions=TRUE)
method <- "svmRadial"
tune.grid <- expand.grid(.sigma = expGrid(-10, -6, 1, 2), .C = c(0.5, 1, 2, 4, 8, 16, 32))
registerDoMC(cores=15)
model <- train(dataset$x.train, dataset$y.train, method, tuneGrid=tune.grid, trControl=dataset$trControl, savePredictions=T)
saveRDS(model, file=paste(method,".rds",sep=""))
#model <- readRDS("svmRadial.rds")

# get the data out of the model and the dataset in preparation for error estimation cross validation
indexes <- which(apply(as.data.frame(model$pred[,names(model$bestTune)]), 1, function(x,y) {identical(as.numeric(x), as.numeric(y))}, model$bestTune))
best <- model$pred[indexes,]
order <- order(best$rowIndex, decreasing=FALSE)
best <- best[order,]
obs <- best$obs
preds <- best$pred
errors <- preds - obs
fold.strings <- unique(best$Resample)
folds <- list()
for(fold in 1:length(fold.strings)) {
  folds[[fold]] <- which(best$Resample == fold.strings[fold])
}
x <- as.matrix(dataset$x.train)

# functions used currently for conformal error estimation
get_error_model <- function(model,x,trControl,method="svm",tune.grid=tune.grid){
  # this works for SVM, To make it general, one has to check the following line:
  predObsCV <- model$pred[which(model$pred$sigma == model$bestTune$sigma & model$pred$C == model$bestTune$C), 1:3]
  predObsCV <- predObsCV[order(predObsCV$rowIndex),]
  error_model <- train(x, abs(predObsCV$pred - predObsCV$obs), method, 
                       tuneGrid=tune.grid, trControl=trControl,
                       savePredictions=T)
  return(error_model)
}
h = function(modell){
  model$pred[which(modell$pred$sigma == modell$bestTune$sigma & modell$pred$C == modell$bestTune$C), 1:3]
}
Conformity <- function(model, error_model,confidence=0.80,data.new){
  if (confidence > 1 || confidence < 0)
    stop("confidence must be between 0 and 1")
  if (is.null(data.new))
    stop("new datapoints are required as input")
  print("Calculating alphas..")
  cat('\n')
  predObsCV <- model$pred[which(model$pred$sigma == model$bestTune$sigma & model$pred$C == model$bestTune$C), 1:3]
  predObsCV <- predObsCV[order(predObsCV$rowIndex),]
  predObsCV_error <- error_model$pred[which(error_model$pred$sigma == error_model$bestTune$sigma & error_model$pred$C == error_model$bestTune$C), 1:3]
  predObsCV_error <- predObsCV_error[order(predObsCV_error$rowIndex),]
  
  conformityA <- function(obs,pred,error){
    alpha <- abs(obs-pred) / error
    return(sort(alpha))
  }
  
  alphas <- conformityA(predObsCV$obs,predObsCV$pred,predObsCV_error$pred)
  print("Predicting for the new data..")
  cat('\n')
  pred <- as.vector(predict(model$finalModel, newdata = data.new))
  pred_error <- as.vector(predict(error_model$finalModel, newdata = data.new)) 
  
  intervals <- ((alphas[length(alphas)*confidence]) * pred_error) 
  return(intervals)
}

# train the various methods over 
# for each left out fold, build an estimator on the other folds
for(lo in 1:NF) {
  nfolds <- list()
  nx <- c()
  nobs <- c()
  nerrors <- c()
  inset <- c(1:NF)[-lo]
  f <- 1
  for(fold in inset) {
    indexes <- folds[[fold]]
    nx <- rbind(nx, x[indexes,])
    start <- length(nobs) + 1
    nobs <- c(nobs, obs[indexes])
    end <- length(nobs)
    nerrors <- c(nerrors, errors[indexes])
    nfolds[[f]] <- start:end
    f <- f + 1
  }
  npreds <- nerrors + nobs
  
  # use the ensemble method
  print("Ensemble")
  dms <- CreateDistanceMatricesMC(nx, nfolds, cores=15)
  m <- CreateEstimatorMatrixMC(20, nfolds, dms, nerrors, nobs, npreds, cores=15)
  sl <- GetSigmaSplineList(m, nerrors)
  weights <- GetWeights(m, nerrors, sl, cores=15)
  estimator <- list(x = nx, sl = sl, weights = weights, obs = nobs, preds = npreds, errors = nerrors, Nmax = 20)
  #saveRDS(estimator, file=paste("estimator", lo, ".rds", sep=""))
  
  # make predictions for the left out fold
  indexes <- folds[[lo]]
  tx <- x[indexes,]
  prediction <- PredictSigmas(tx, estimator)
  saveRDS(prediction, file=paste("prediction", lo, ".rds", sep=""))
  
  
  # use 10 fold cross validation for each holdout set
  trControl <- trainControl(method = "cv", number = 10, 
                            repeats = 1, returnResamp = "none", returnData = TRUE, 
                            savePredictions = TRUE, verboseIter = TRUE, allowParallel = TRUE, 
                            index = createMultiFolds(nobs, k = 10, 
                                                     times = 1))
  
  # do the conformity predictions
  print("SVM")
  tune.grid <- expand.grid(.sigma = expGrid(-10, -6, 1, 2), .C = c(0.5, 1, 2, 4, 8, 16, 32))
  registerDoMC(cores=15)
  model <- train(nx, nobs, method="svmRadial", tuneGrid=tune.grid, trControl=trControl, savePredictions=T)
  saveRDS(model, file=paste("svm_", lo, ".rds",sep=""))
  error_model <- get_error_model(model, nx, trControl, method="svmRadial", tune.grid=tune.grid)
  saveRDS(error_model, file=paste("error_model_", lo, ".rds",sep=""))
  conformity_intervals <- Conformity(model, error_model, confidence=0.80, data.new=tx)
  saveRDS(conformity_intervals, file=paste("conformity_interval_", lo, ".rds",sep=""))
  
  # train the GP models
  print("GP")
  registerDoMC(cores=15)
  tune.grid <- expand.grid(.sigma = expGrid(-10, -6, 1, 2))
  gp <- train(nx, nobs, method = "gaussprRadial", trControl = trControl, type = "regression", variance.model=TRUE, tuneGrid = tune.grid)
  saveRDS(gp, file=paste("gp_", lo, ".rds",sep=""))
  
  # train the RF models
  print("RF")
  registerDoMC(cores=15)
  tune.grid <- expand.grid(.mtry = floor(seq(10, 140, length.out = 15)))
  rf <- train(nx, nobs, method = "rf", trControl = trControl, importance=TRUE, tuneGrid = tune.grid)
  saveRDS(rf, file=paste("rf_", lo, ".rds",sep=""))
}


# either load or calculate the variance measures here
CECs <- matrix(nrow=98, ncol=NF)
ensemble <- rep(NA, NF)
conformity <- rep(NA, NF)
gpCEC <- rep(NA, NF)
rfCEC <- rep(NA, NF)
for(lo in 1:NF) {
  indexes <- folds[[lo]]
  terrors <- errors[indexes]
  prediction <- readRDS(file=paste("prediction", lo, ".rds", sep=""))
  CECs[,lo] <- apply(prediction$sigma.matrix, 2, CEC, abs(terrors))
  ensemble[lo] <- CEC(prediction$sigmas, abs(terrors))
  conformity_intervals <- readRDS(file=paste("conformity_interval_", lo, ".rds",sep=""))
  conformity[lo] <- CEC(conformity_intervals, abs(terrors))
  tx <- x[indexes,]
  gp <- readRDS(file=paste("gp_", lo, ".rds",sep=""))
  gp_var <- as.vector(predict(gp$finalModel, newdata = tx, type="variance"))
  gpCEC[lo] <- CEC(gp_var, abs(terrors))
  rf <- readRDS(file=paste("rf_", lo, ".rds",sep=""))
  rf_preds <- predict(rf$finalModel, newdata = tx, predict.all=TRUE)
  rf_var <- apply(X=rf_preds$individual, MARGIN=1, FUN=var)
  rfCEC[lo] <- CEC(rf_var, abs(terrors))
  plot(rf_var, abs(terrors))
}

# put the variance measure vectors into a matrix to use for plotting
aves <- apply(CECs, 1, mean)
aves <- c(aves, rep(mean(ensemble), 20))
aves <- c(aves, rep(mean(conformity), 20))
aves <- c(aves, rep(mean(gpCEC), 20))
aves <- c(aves, rep(mean(rfCEC), 20))
sds <- apply(CECs, 1, sd)
sds <- c(sds, rep(sd(ensemble), 20))
sds <- c(sds, rep(sd(conformity), 20))
sds <- c(sds, rep(sd(gpCEC), 20))
sds <- c(sds, rep(sd(rfCEC), 20))
estimators <- c("CONFINE", "VarError", "CONVIVE", "AvgDist", "DiffNN")
method <- rep(estimators, 20)[-c(2,3)]
method <- c(method, rep("Ensemble", 20))
method <- c(method, rep("Conformity", 20))
method <- c(method, rep("GP Var", 20))
method <- c(method, rep("RF Var", 20))
N <- floor(((1:100)+4)/5)[-c(2,3)]
N <- c(N, 1:20)
N <- c(N, 1:20)
N <- c(N, 1:20)
N <- c(N, 1:20)

# checking lengths
length(N)
length(aves)
length(sds)
length(method)
length(N)

pd <- data.frame(aves = aves, sds = sds, method = method, N = N)

f1 = ggplot(pd, aes(x = N, y = aves, group = method, col=method) ) + # method becomes a classifying factor 
  geom_errorbar(aes(ymin = aves - sds, ymax = aves + sds), width=0.3) + # add error bars (do so before geom_point so the points are on top of the error bars)
  geom_line() + # join points with lines (specify this before geom_point, or the lines will be drawn over the shapes)
  geom_point(aes(shape=method, fill=method), size=5) + # add a scatterplot; constant size, shape/fill depends on lesion
  scale_x_continuous("Number of Neighbours", breaks=1:20) + # have tick marks for each session
  scale_y_continuous("Average CEC +-1 sd") + # rescale Y axis slightly
  theme_bw() + # make the theme black-and-white rather than grey (do this before font changes, or it overrides them)
  theme(plot.title = element_text(face="bold", size=14), # use theme_get() to see available options
        axis.title.x = element_text(face="bold", size=12),
        axis.title.y = element_text(face="bold", size=12, angle=90),
        panel.grid.major = element_blank(), # switch off major gridlines
        panel.grid.minor = element_blank(), # switch off minor gridlines
        #legend.position = c(0.2,0.8), # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
        legend.title = element_blank(), # switch off the legend title
        legend.text = element_text(size=12),
        legend.key.size = unit(1.5, "lines"),
        legend.key = element_blank() # switch off the rectangle around symbols in the legend
  )
print(f1)




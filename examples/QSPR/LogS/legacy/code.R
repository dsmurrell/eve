training.sigmas <- PredictSigmas(x=dataset$x.train, estimator=eve)$sigmas
training.errors <- model$pred  dataset$y.train
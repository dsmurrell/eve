ConvertToCISplines <- function(x, percent, spar = 0.8) {
  splines <- list()
  order <- order(x, decreasing=TRUE)
  error.sort <- errors[order]
  abs.error.sort <- abs(error.sort)
  x.sort <- x[order]
  hwidth <- floor(length(x)/10)
  cis <- rep(NA, length(x))
  inside <- hwidth:(length(x)-hwidth)
  offset <- (100 - percent)/200
  for(i in inside) {
    set <- abs.error.sort[(i-hwidth):(i+hwidth)]
    qs <- quantile(set, c(offset, 1-offset))
    cis[i] <- qs[2] - qs[1]
  }
  
  splines$index.spline <- smooth.spline(x.sort, 1:length(x), spar = spar)
  splines$ci.spline <- smooth.spline(inside, cis[inside], spar = spar)
  return(splines)
}

############## CRAP

# mlist <- lapply(apply(m, 2, FUN=list), unlist)
# length(mlist)
# mlist[[1]]
# 
# sapply((1:3), FUN=list)
# sapply((1:ncol(m))
# sigma.matrix <- do.call(cbind, mclapply(1:ncol(m), GetCIsFromSplines, m, sl, mc.cores=1))
# dim(sigma.matrix)
# GetCIsFromSplineList(i=1, 
# 
# CEC(sigmas, abs(errors))
# sigmas <- apply(sweep(sigma.matrix,2,weights,`*`), 1, sum)
# dim(sigma.matrix)
# plot(errors, sigmas)
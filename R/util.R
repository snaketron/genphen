



# Description:
# Computes HDI given a vector, taken "Doing Bayesian Analysis"
getHdi <- function(vec, hdi.level) {
  sortedPts <- sort(vec)
  ciIdxInc <- floor(hdi.level * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  return(HDIlim)
}




# Description:
# Given a confusion matrix table(predicted, real), compute the Cohen's
# kappa statistics. Cohen makes the following distinction between the
# different kappa ranges.
getKappa <- function(predicted, real, aas) {
  
  # should not occur, just in case check this
  if(length(unique(aas)) != length(unique(c(predicted, real)))) {
    stop("Error while building confusion matrix (getKappa)")
  }
  
  # Build 2x2 confusion matrix
  buildConfusionMatrix <- function(predicted, real) {
    cm <- matrix(data = 0, nrow = 2, ncol = 2)
    cm[1, 1] <- length(intersect(which(real %in% aas[1]),
                                 which(predicted %in% aas[1])))
    cm[2, 2] <- length(intersect(which(real %in% aas[2]),
                                 which(predicted %in% aas[2])))
    cm[2, 1] <- length(intersect(which(real %in% aas[1]),
                                 which(!predicted %in% aas[1])))
    cm[1, 2] <- length(intersect(which(real %in% aas[2]),
                                 which(!predicted %in% aas[2])))
    return (cm)
  }
  
  cm <- buildConfusionMatrix(predicted = predicted, real = real)
  
  ca.exp <- (sum(cm[1, ])*sum(cm[, 1])+sum(cm[2, ])*sum(cm[, 2]))/sum(cm)^2
  ca <- (cm[1, 1]+cm[2, 2])/sum(cm)
  kappa <- (ca-ca.exp)/(1-ca.exp)
  
  # if NaN, CA_exp = 1
  if(is.nan(x = kappa) == TRUE) {
    kappa <- 0
  }
  
  return (kappa)
}




# Description:
# Bhattacharyya Coefficient of two distribution
# Taken from: source("http://tguillerme.github.io/R/bhatt.coef.R")
getBhattacharyya <- function(x, y, bw = bw.nrd0, ...) {
  #SANITIZING
  #x
  if(class(x) != 'numeric') {
    stop("'x' must be numeric.")
  }
  if(length(x) < 2) {
    stop("'x' need at least two data points.")
  }
  
  #y
  if(class(y) != 'numeric') {
    stop("'y' must be numeric.")
  }
  if(length(y) < 2) {
    stop("'y' need at least two data points.")
  }
  
  #bw
  if(length(bw) != 1) {
    stop("'bw' must be either a single numeric value or a single function.")
  }
  if(class(bw) != 'function') {
    if(class(bw) != 'numeric') {
      stop("'bw' must be either a single numeric value or a single function.")
    }
  }
  #Avoiding non-entire numbers
  if(class(bw) == 'numeric') {
    bw<-round(bw)
  }
  
  #BHATTACHARYYA COEFFICIENT
  #sum(sqrt(x relative counts in bin_i * y relative counts in bin_i))
  
  #Setting the right number of bins (i)
  if(class(bw) == 'function') {
    #Bin width
    band.width<-bw(c(x,y), ...)
    #Bin breaks
    #adding an extra bandwith to the max to be sure to include all the data
    bin.breaks<-seq(from=min(c(x,y)), to=max(c(x,y)+band.width), by=band.width)
    #Number of bins
    bin.n<-length(bin.breaks)-1
  } else {
    #Bin breaks
    bin.breaks<-hist(c(x,y), breaks=bw, plot=FALSE)$breaks
    #Bin width
    band.width<-diff(bin.breaks)[1]
    #Number of bins
    bin.n<-bw
  }
  
  #Counting the number of elements per bin
  histx<-hist(x, breaks=bin.breaks, plot=FALSE)[[2]]
  histy<-hist(y, breaks=bin.breaks, plot=FALSE)[[2]]
  #Relative counts
  rel.histx<-histx/sum(histx)
  rel.histy<-histy/sum(histy)
  
  #Calculating the Bhattacharyya Coefficient (sum of the square root of
  # the multiple of the relative counts of both distributions)
  bc <- sum(sqrt(rel.histx*rel.histy))
  return(list(bc = bc))
}




# Description:
# Given phenotype.type and model.type, the procedure compiles the STAN model.
compileModel <- function(phenotype.type, model.type) {
  if(phenotype.type == "continuous") {
    f.local <- paste("inst/extdata/continuous", model.type, "stan", sep = '.')
    f.pkg <- system.file("extdata", package = "genphen",
                         paste("continuous", model.type, "stan", sep = '.'))
  }
  else if(phenotype.type == "dichotomous") {
    f.local <- paste("inst/extdata/dichotomous", model.type, "stan", sep = '.')
    f.pkg <- system.file("extdata", package = "genphen",
                         paste("dichotomous", model.type, "stan", sep = '.'))
  }
  
  if(file.exists(f.pkg)) {
    model.stan <- rstan::stan_model(file = f.pkg, model_name = "model")
  }
  if(file.exists(f.local)) {
    model.stan <- rstan::stan_model(file = f.local, model_name = "model")
  }
  
  return(model.stan)
}




# Description:
# If an object of type DNAMultipleAlignment
convertMsaToGenotype <- function(genotype) {
  if(is.null(attr(genotype, "class")) == FALSE) {
    genotype <- as.matrix(genotype)
  }
  return (genotype)
}




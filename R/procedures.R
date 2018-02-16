


# Description:
# Get genotype-phenotype data in format for stan
getGenphenData <- function(genotype, phenotype, phenotype.type,
                           min.observations = 2) {

  if(phenotype.type == "continuous") {
    out <- vector(mode = "list")
    out.counter <- 1
    max.sigma <- 0 # max sigma throughout dataset
    for(i in 1:ncol(genotype)) {
      js <- unique(genotype[, i])
      if(length(js) != 1) {
        X <- as.numeric(as.factor(genotype[, i]))
        k <- which(table(X) >= min.observations)
        if(length(k) != 1) {
          Y <- phenotype[X %in% as.numeric(names(k))]
          G <- genotype[X %in% as.numeric(names(k)), i]
          X <- as.numeric(as.factor(X[X %in% as.numeric(names(k))]))

          Ng <- numeric(length = length(unique(X))) # nr. of samples per genotype
          E_mu <- numeric(length = length(unique(X)))
          E_sigma <- numeric(length = length(unique(X)))
          for(j in 1:max(X)) {
            temp.Y <- Y[X == j]
            Ng[j] <- length(temp.Y)
            E_mu[j] <- mean(temp.Y)
            E_sigma[j] <- stats::sd(temp.Y)

            # max observed SD in the dataset, use to set as empirical SD if for
            # a given group only a single observation is avaliable or SD = 0.
            max.sigma <- max(max.sigma, max(c(E_sigma, 0), na.rm = TRUE))
          }

          l <- list(site = i, G = G, Y = Y, X = X, Ng = Ng, Nx = length(Ng),
                    Ny = length(Y), E_mu = E_mu, E_sigma = E_sigma)
          out[[out.counter]] <- l
          out.counter <- out.counter + 1
        }
      }
    }

    # if empty return NULL
    if(length(out) == 0) {
      return(NULL)
    }

    # correct empirical SD
    for(i in 1:length(out)) {
      for(j in 1:length(out[[i]]$E_sigma)) {
        if(is.na(out[[i]]$E_sigma[j]) | out[[i]]$E_sigma[j] == 0) {
          out[[i]]$E_sigma[j] <- max.sigma
        }
      }
    }
  }
  else if(phenotype.type == "dichotomous") {
    out <- vector(mode = "list")
    out.counter <- 1
    for(i in 1:ncol(genotype)) {
      js <- unique(genotype[, i])
      if(length(js) != 1) {
        X <- as.numeric(as.factor(genotype[, i]))
        k <- which(table(X) >= min.observations)
        if(length(k) != 1) {
          Y <- phenotype[X %in% as.numeric(names(k))]
          G <- genotype[X %in% as.numeric(names(k)), i]
          X <- as.numeric(as.factor(X[X %in% as.numeric(names(k))]))

          Ng <- numeric(length = length(unique(X))) # nr. of samples per genotype
          for(j in 1:max(X)) {
            temp.Y <- Y[X == j]
            Ng[j] <- length(temp.Y)
          }

          l <- list(site = i, G = G, Y = Y, X = X, Ng = Ng,
                    Nx = length(Ng), Ny = length(Y))
          out[[out.counter]] <- l
          out.counter <- out.counter + 1
        }
      }
    }
  }

  # if empty return NULL
  if(length(out) == 0) {
    return(NULL)
  }


  return (out)
}




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
# Given a vector, this procedure computes its mean and CI with bootstrapping.
getMeanHdi <- function(vec, hdi.level) {
  boot.mean <- function(data, index) {
    mean(data[index])
  }
  b <- boot::boot(data = vec, statistic = boot.mean, R = 1000)
  return (list(mu = b$t[, 1]))
}




# Description:
# If an object of type DNAMultipleAlignment
convertMsaToGenotype <- function(genotype) {
  if(is.null(attr(genotype, "class")) == FALSE) {
    genotype <- as.matrix(genotype)
  }
  return (genotype)
}




# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInput <- function(genotype, phenotype, phenotype.type, mcmc.chains,
                       mcmc.iterations, mcmc.warmup, mcmc.cores, hdi.level,
                       stat.learn.method, cv.iterations, diagnostics.points,
                       diagnostics.samples, diagnostics.rf.trees) {

  checkGenotypePhenotype <- function(genotype, phenotype) {
    # CHECK: genotype
    if(is.null(attr(genotype, "class")) == FALSE) {
      if(!attr(genotype, "class") %in% c("AAMultipleAlignment",
                                         "DNAMultipleAlignment")) {
        stop("genotype can be one of the following structures: matrix or
             data.frame or DNAMultipleAlignment/AAMultipleAlignment")
      }
      else {
        temp <- as.matrix(genotype)
        if(nrow(temp) < 2 | ncol(temp) == 0) {
          stop("the genotypes cannot have less than two observations the or
               number of genotypes cannot be 0.")
        }
      }
    }
    else {
      if(is.vector(genotype)) {
        genotype <- matrix(data = genotype, ncol = 1)
      }

      if(nrow(genotype) < 2 | ncol(genotype) == 0) {
        stop("the genotypes cannot have less than two observations or the
             number of genotypes cannot be 0.")
      }

      if(!is.matrix(genotype) & !is.data.frame(genotype)) {
        stop("genotype can be one of the following structures: matrix or
             data.frame or DNAMultipleAlignment/AAMultipleAlignment")
      }

      if(typeof(genotype) != "character") {
        stop("if it is structured in matrix/data.frame the genotype must
           be of character type.")
      }
    }

    # CHECK: phenotype
    if(!is.vector(phenotype)) {
      stop("the phenotype must be a vector.")
    }

    if(length(phenotype) < 2) {
      stop("phenotype cannot contain fewer than 2 elements.")
    }

    if(!is.numeric(phenotype)) {
      stop("the phenotype must be of numeric type.")
    }

    # CHECK: genotype & phenotype
    if(nrow(genotype) != length(phenotype)) {
      stop("length(genotype)!=length(phenotype), they must be equal in length.")
    }
  }

  checkPhenotypeValidity <- function(phenotype, phenotype.type) {
    if(phenotype.type == "dichotomous") {
      if(length(unique(phenotype)) != 2) {
        stop("The dichotomous phenotype must be a vector with exactly two
             categories (classes) \n")
      }
    }

    if(phenotype.type == "continuous") {
      if(length(unique(phenotype)) <= 2) {
        warning("The continuous phenotype has less then 3 unique elements,
                are you sure this is a continuous vector? \n")
      }
    }
  }

  checkPhenotypeType <- function(phenotype.type) {
    # CHECK: phenotype.type
    if(length(phenotype.type) != 1) {
      stop("phenotype.type must be a string (default = 'continuous')")
    }

    if(!is.character(phenotype.type)) {
      stop("phenotype.type must be a string: 'continuous' or 'dichotomous'")
    }

    if(!phenotype.type %in% c("continuous", "dichotomous")) {
      stop("phenotype.type must be a string: 'continuous' or 'dichotomous'")
    }
  }

  checkMcmcIterations <- function(mcmc.iterations) {
    # CHECK: mcmc.iterations
    if(length(mcmc.iterations) != 1) {
      stop("the mcmc.iterations must be a number > 0 (default = 10000).")
    }

    if(!is.numeric(mcmc.iterations)) {
      stop("mcmc.iterations must be a numeric argument (default = 10000).")
    }

    if(mcmc.iterations <= 0) {
      stop("mcmc.iterations must be larger than 0 (default = 10000).")
    }
  }

  checkMcmcWarmup <- function(mcmc.warmup) {
    # CHECK: mcmc.warmup
    if(length(mcmc.warmup) != 1) {
      stop("the mcmc.warmup must be a number > 0 (default = 5000).")
    }

    if(!is.numeric(mcmc.warmup)) {
      stop("mcmc.warmup must be a numeric argument (default = 5000).")
    }

    if(mcmc.warmup <= 0) {
      stop("mcmc.warmup must be larger than 0 (default = 5000).")
    }
  }

  checkMcmcChains <- function(mcmc.chains) {
    # CHECK: mcmc.chains
    if(length(mcmc.chains) != 1) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }

    if(!is.numeric(mcmc.chains)) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }

    if(mcmc.chains <= 0) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }
  }

  checkMcmcCores <- function(mcmc.cores) {
    # CHECK: mcmc.cores
    if(length(mcmc.cores) != 1) {
      stop("mcmc.cores is numeric parameter.")
    }

    if(is.numeric(mcmc.cores) == FALSE) {
      stop("mcmc.cores is numeric parameter.")
    }

    if(mcmc.cores <= 0) {
      stop("mcmc.cores is numeric parameter >=1.")
    }
  }

  checkHdi <- function(hdi.level) {
    if(length(hdi.level) != 1) {
      stop("The HDI level must be in range (0, 1).")
    }

    if(is.numeric(hdi.level) == FALSE) {
      stop("The HDI level must be in range (0, 1).")
    }

    if(hdi.level >= 1 | hdi.level <= 0) {
      stop("The HDI level must be in range (0, 1).")
    }
  }

  checkMlMethod <- function(stat.learn.method) {
    # CHECK: phenotype.type
    if(length(stat.learn.method) != 1) {
      stop("stat.learn.method must be a string (default = 'rf')")
    }

    if(!is.character(stat.learn.method)) {
      stop("stat.learn.method must be a string (default = 'rf')")
    }

    if(!stat.learn.method %in% c("rf", "svm")) {
      stop("stat.learn.method must be a string (either 'rf' or 'svm')")
    }
  }

  checkCv <- function(cv.iterations) {
    if(length(cv.iterations) != 1) {
      stop("cv.iterations must be a number (default = 1,000).")
    }

    if(is.numeric(cv.iterations) == FALSE) {
      stop("cv.iterations must be a number (default = 1,000).")
    }

    if(cv.iterations < 500) {
      stop("cv.iterations >= 500 recomended (default = 1,000).")
    }
  }

  checkDiagnostics <- function(diagnostics.points,
                               diagnostics.samples,
                               diagnostics.rf.trees) {
    if(length(diagnostics.points) != 1) {
      stop("diagnostics.points must be a number (default = 10).")
    }
    if(length(diagnostics.samples) != 1) {
      stop("diagnostics.samples must be a number (default = 10).")
    }
    if(length(diagnostics.rf.trees) != 1) {
      stop("diagnostics.samples must be a number (default = 50,000).")
    }


    if(is.numeric(diagnostics.points) == FALSE) {
      stop("diagnostics.points must be a number (default = 10).")
    }
    if(is.numeric(diagnostics.samples) == FALSE) {
      stop("diagnostics.samples must be a number (default = 10).")
    }
    if(is.numeric(diagnostics.rf.trees) == FALSE) {
      stop("diagnostics.rf.trees must be a number (default = 50,000).")
    }


    if(diagnostics.points <= 1) {
      stop("diagnostics.points >= 2 accepted (default = 10).")
    }
    if(diagnostics.samples <= 0) {
      stop("diagnostics.samples >= 1 accepted (default = 10).")
    }
    if(diagnostics.rf.trees <= 10000) {
      stop("diagnostics.rf.trees >= 10,000 accepted (default = 50,000).")
    }
  }

  if(is.null(genotype) | missing(genotype) |
     is.null(phenotype) | missing(phenotype) |
     is.null(phenotype.type) | missing(phenotype.type) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.iterations) | missing(mcmc.iterations) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(mcmc.cores) | missing(mcmc.cores) |
     is.null(hdi.level) | missing(hdi.level) |
     is.null(stat.learn.method) | missing(stat.learn.method) |
     is.null(cv.iterations) | missing(cv.iterations) |
     is.null(diagnostics.points) | missing(diagnostics.points) |
     is.null(diagnostics.samples) | missing(diagnostics.samples) |
     is.null(diagnostics.rf.trees) | missing(diagnostics.rf.trees)) {
    stop("arguments must be non-NULL/specified")
  }

  checkGenotypePhenotype(genotype = genotype, phenotype = phenotype)
  checkPhenotypeType(phenotype.type = phenotype.type)
  checkPhenotypeValidity(phenotype = phenotype, phenotype.type = phenotype.type)
  checkMcmcIterations(mcmc.iterations = mcmc.iterations)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkMcmcCores(mcmc.cores = mcmc.cores)
  checkHdi(hdi.level = hdi.level)
  checkMlMethod(stat.learn.method = stat.learn.method)
  checkDiagnostics(diagnostics.points = diagnostics.points,
                   diagnostics.samples = diagnostics.samples,
                   diagnostics.rf.trees = diagnostics.rf.trees)
}





# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using random forests.
getRfCa <- function(data.list, cv.fold, cv.steps, hdi.level, ntree) {


  # Description:
  # Perform N number of classifications and compute N number of:
  # - classification accuracy
  # - number of successful classifications (ideally == N)
  booter <- function(X, Y, cv.fold, cv.steps, ntree) {
    # number of total data entries
    rows <- length(Y)

    # initialize result vectors
    accuracies <- c()
    successful.boots <- 0

    # get the size of the smaller sample
    unique.Y <- unique(Y)
    min.sample <- min(table(Y))
    cv.sample <- floor(x = min.sample * cv.fold)
    ys.1 <- which(Y == unique.Y[1])
    ys.2 <- which(Y == unique.Y[2])

    for(b in 1:cv.steps) {
      # sample at random
      s.1 <- sample(x = ys.1, size = min.sample, replace = FALSE)
      s.2 <- sample(x = ys.2, size = min.sample, replace = FALSE)

      train <- data.frame(Y = c(Y[s.1[1:cv.sample]], Y[s.2[1:cv.sample]]),
                          X = c(X[s.1[1:cv.sample]], X[s.2[1:cv.sample]]),
                          stringsAsFactors = FALSE)
      test <- data.frame(Y = c(Y[s.1[(cv.sample + 1):min.sample]],
                               Y[s.2[(cv.sample + 1):min.sample]]),
                         X = c(X[s.1[(cv.sample + 1):min.sample]],
                               X[s.2[(cv.sample + 1):min.sample]]),
                         stringsAsFactors = FALSE)

      # only one type of predictor (no continous variable)
      if(length(unique(train$X)) <= 1) {
        accuracies <- c(accuracies, NA)
      }
      else {
        # train classification model (try condition to avoid errors in case
        # only one-level predictor is train data)
        rf.out <- try(randomForest::randomForest(as.factor(Y) ~ X,
                                                 data = train,
                                                 ntree = ntree),
                      silent = TRUE)
        if(attr(rf.out, "class")[1] == "try-error") {
          accuracies <- c(accuracies, NA)
        }
        else {
          # test classification model
          prediction <- stats::predict(object = rf.out, newdata = test)

          # compute classification accuracy (1 - classification error)
          ac <- length(
            which(as.character(test$Y)==as.character(prediction)))/nrow(test)
          accuracies <- c(accuracies, ac)

          # record bootstrap
          successful.boots <- successful.boots + 1
        }
      }
    }

    # build result list and return
    result <- list(accuracies = accuracies,
                   successful.boots = successful.boots)
    return (result)
  }


  ca.out <- c()
  for(i in 1:(max(data.list$X) - 1)) {
    for(j in (i + 1):max(data.list$X)) {


      # general data
      site <- data.list$site
      n.i <- data.list$Ng[i]
      n.j <- data.list$Ng[j]
      g.i <- data.list$G[data.list$X == i][1]
      g.j <- data.list$G[data.list$X == j][1]
      mutation <- paste(g.i, "->", g.j, sep = '')
      general <- paste(g.i, ":", n.i, ", ", g.j, ":", n.j, sep = '')


      # subset predictor/response
      X <- data.list$Y[data.list$X %in% c(i, j)]
      Y <- as.factor(data.list$X[data.list$X %in% c(i, j)])


      # perform classification
      class.obj <- booter(Y = Y, X = X, cv.fold = cv.fold,
                          cv.steps = cv.steps, ntree = 1000)


      # get 95% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = class.obj$accuracies, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      ca.hdi <- paste("(", round(x = ca.L, digits = 2), ", ",
                      round(x = ca.H, digits = 2), ")", sep = '')


      # pack outputs
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          ca = mean(class.obj$accuracies, na.rm = TRUE),
                          ca.L = ca.L,
                          ca.H = ca.H,
                          ca.hdi = ca.hdi,
                          ca.boots = class.obj$successful.boots)
      ca.out <- rbind(ca.out, stats)
    }
  }

  return(ca.out)
}




# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using support vector machines.
getSvmCa <- function(data.list, cv.fold, cv.steps, hdi.level) {


  # Description:
  # Perform N number of classifications and compute N number of:
  # - classification accuracy
  # - number of successful classifications (ideally == N)
  booter <- function(X, Y, cv.fold, cv.steps, ntree) {
    # number of total data entries
    rows <- length(Y)

    # initialize result vectors
    accuracies <- c()
    successful.boots <- 0

    # get the size of the smaller sample
    unique.Y <- unique(Y)
    min.sample <- min(table(Y))
    cv.sample <- floor(x = min.sample * cv.fold)
    ys.1 <- which(Y == unique.Y[1])
    ys.2 <- which(Y == unique.Y[2])

    for(b in 1:cv.steps) {
      # sample at random
      s.1 <- sample(x = ys.1, size = min.sample, replace = FALSE)
      s.2 <- sample(x = ys.2, size = min.sample, replace = FALSE)

      train <- data.frame(Y = c(Y[s.1[1:cv.sample]], Y[s.2[1:cv.sample]]),
                          X = c(X[s.1[1:cv.sample]], X[s.2[1:cv.sample]]),
                          stringsAsFactors = FALSE)
      test <- data.frame(Y = c(Y[s.1[(cv.sample + 1):min.sample]],
                               Y[s.2[(cv.sample + 1):min.sample]]),
                         X = c(X[s.1[(cv.sample + 1):min.sample]],
                               X[s.2[(cv.sample + 1):min.sample]]),
                         stringsAsFactors = FALSE)

      # only one type of predictor (no continous variable)
      if(length(unique(train$X)) <= 1) {
        accuracies <- c(accuracies, NA)
      }
      else {
        # train classification model (try condition to avoid errors in case
        # only one-level predictor is train data)
        svm.out <- try(e1071::svm(as.factor(Y) ~ X,
                                  data = train,
                                  type = "C-classification"),
                       silent = TRUE)
        if(attr(svm.out, "class")[1] == "try-error") {
          accuracies <- c(accuracies, NA)
        }
        else {
          # test classification model
          prediction <- stats::predict(object = svm.out, newdata = test)

          # compute classification accuracy (1 - classification error)
          ac <- length(
            which(as.character(test$Y)==as.character(prediction)))/nrow(test)

          # record bootstrap
          successful.boots <- successful.boots + 1
        }
      }
    }

    # build result list and return
    result <- list(accuracies = accuracies, successful.boots = successful.boots)
    return (result)
  }


  ca.out <- c()
  for(i in 1:(max(data.list$X) - 1)) {
    for(j in (i + 1):max(data.list$X)) {


      # general data
      site <- data.list$site
      n.i <- data.list$Ng[i]
      n.j <- data.list$Ng[j]
      g.i <- data.list$G[data.list$X == i][1]
      g.j <- data.list$G[data.list$X == j][1]
      mutation <- paste(g.i, "->", g.j, sep = '')
      general <- paste(g.i, ":", n.i, ", ", g.j, ":", n.j, sep = '')


      # subset predictor/response
      X <- data.list$Y[data.list$X %in% c(i, j)]
      Y <- as.factor(data.list$X[data.list$X %in% c(i, j)])


      # perform classification
      class.obj <- booter(Y = Y, X = X, cv.fold = cv.fold,
                          cv.steps = cv.steps, ntree = 1000)


      # get HDI for the classification accuracy
      ca.hdi <- getHdi(vec = class.obj$accuracies, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      ca.hdi <- paste("(", round(x = ca.L, digits = 2), ", ",
                      round(x = ca.H, digits = 2), ")", sep = '')

      # pack outputs
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          ca = mean(class.obj$accuracies, na.rm = TRUE),
                          ca.L = ca.L,
                          ca.H = ca.H,
                          ca.hdi = ca.hdi,
                          ca.boots = class.obj$successful.boots)
      ca.out <- rbind(ca.out, stats)
    }
  }

  return(ca.out)
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

  #Calculating the Bhattacharyya Coefficient (sum of the square root of the multiple of the relative counts of both distributions)
  bhatt.coeff<-sum(sqrt(rel.histx*rel.histy))
  return(bhatt.coeff)
  #End
}




# Description:
# Computes a Bayesian t-test
runContinuous <- function(data.list, mcmc.chains, mcmc.iterations, mcmc.warmup,
                          mcmc.cores, hdi.level, model.stan) {

  # get initial parameter values
  posterior <- sampling(object = model.stan,
                        data = data.list,
                        pars = c("mu", "sigma", "nu"),
                        iter = mcmc.iterations,
                        warmup = mcmc.warmup,
                        chains = mcmc.chains,
                        cores = mcmc.cores,
                        control = list(adapt_delta = 0.99, max_treedepth = 10),
                        verbose = FALSE,
                        refresh = -1)


  # check for divergence
  divergence.params <- get_sampler_params(posterior, inc_warmup = FALSE)
  divergence <- FALSE
  treedepth <- FALSE
  for(chain in 1:mcmc.chains) {
    if(any(divergence.params[[chain]][, "divergent__"] != 0)) {
      divergence <- TRUE
    }

    if(any(divergence.params[[chain]][, "treedepth__"] > 10)) {
      treedepth <- TRUE
    }
  }


  # convergence data
  convergence <- summary(posterior)$summary[, c("Rhat", "n_eff")]
  convergence.out <- c()
  for(i in 1:max(data.list$X)) {
    mu.rhat <- convergence[paste("mu[", i, "]", sep = ''), "Rhat"]
    mu.ess <- convergence[paste("mu[", i, "]", sep = ''), "n_eff"]
    sigma.rhat <- convergence[paste("sigma[", i, "]", sep = ''), "Rhat"]
    sigma.ess <- convergence[paste("sigma[", i, "]", sep = ''), "n_eff"]
    convergence.row <- data.frame(s = data.list$site,
                                  g = data.list$G[data.list$X == i][1],
                                  n = sum(data.list$X == i),
                                  mu.rhat = mu.rhat, sigma.rhat = sigma.rhat,
                                  mu.ess = mu.ess, sigma.ess = sigma.ess,
                                  divergence = divergence,
                                  treedepth = treedepth,
                                  stringsAsFactors = FALSE)
    convergence.out <- rbind(convergence.out, convergence.row)
  }


  # posterior data
  posterior <- data.frame(extract(posterior))
  statistics.out <- c()
  for(i in 1:(max(data.list$X) - 1)) {
    for(j in (i + 1):max(data.list$X)) {
      # general data
      site <- data.list$site
      n.i <- data.list$Ng[i]
      n.j <- data.list$Ng[j]
      g.i <- data.list$G[data.list$X == i][1]
      g.j <- data.list$G[data.list$X == j][1]
      mutation <- paste(g.i, "->", g.j, sep = '')
      general <- paste(g.i, ":", n.i, ", ", g.j, ":", n.j, sep = '')


      # extract posterior
      mu.i <- posterior[, paste("mu.", i, sep = '')]
      sigma.i <- posterior[, paste("sigma.", i, sep = '')]
      mu.j <- posterior[, paste("mu.", j, sep = '')]
      sigma.j <- posterior[, paste("sigma.", j, sep = '')]
      nu <- posterior[, "nu"]


      # alternative Cohen's d based on Hedges 1981
      # pool.sd <- sqrt(((sigma.i^2)*(n.i-1)+(sigma.j^2)*(n.j-1))/(n.i+n.j-2))
      # cohens.d <- (mu.i - mu.j)/pool.sd
      cohens.d <- (mu.i - mu.j)/sqrt((sigma.i^2 + sigma.j^2)/2)
      cohens.d.mean <- mean(cohens.d)
      cohens.d.hdi <- getHdi(vec = cohens.d, hdi.level = hdi.level)
      cohens.d.L = cohens.d.hdi[1]
      cohens.d.H = cohens.d.hdi[2]
      cohens.d.hdi <- paste("(", round(x = cohens.d.L, digits = 2), ", ",
                            round(x = cohens.d.H, digits = 2), ")", sep = '')


      # Bhat coeff
      ppc.i <- mean(mu.i)+mean(sigma.i)*stats::rt(n = 10^4, df = mean(nu))
      ppc.j <- mean(mu.j)+mean(sigma.j)*stats::rt(n = 10^4, df = mean(nu))
      b.coef <- getBhattacharyya(x = ppc.i, y = ppc.j)


      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          cohens.d = cohens.d.mean,
                          cohens.d.L = cohens.d.L,
                          cohens.d.H = cohens.d.H,
                          cohens.d.hdi = cohens.d.hdi,
                          b.coef = b.coef)
      statistics.out <- rbind(statistics.out, stats)
    }
  }

  return (list(statistics.out = statistics.out,
               convergence.out = convergence.out))
}




# Description:
# Computes a Bayesian odds-ratio test
runDichotomous <- function(data.list, mcmc.chains, mcmc.iterations, mcmc.warmup,
                           mcmc.cores, hdi.level, model.stan) {


  # get initial parameter values
  posterior <- sampling(object = model.stan,
                        data = data.list,
                        pars = c("mu"),
                        iter = mcmc.iterations,
                        warmup = mcmc.warmup,
                        chains = mcmc.chains,
                        cores = mcmc.cores,
                        control = list(adapt_delta = 0.99, max_treedepth = 10),
                        verbose = FALSE,
                        refresh = -1)


  # check for divergence
  divergence.params <- get_sampler_params(posterior, inc_warmup = FALSE)
  divergence <- FALSE
  treedepth <- FALSE
  for(chain in 1:mcmc.chains) {
    if(any(divergence.params[[chain]][, "divergent__"] != 0)) {
      divergence <- TRUE
    }

    if(any(divergence.params[[chain]][, "treedepth__"] > 10)) {
      treedepth <- TRUE
    }
  }


  # convergence data
  convergence <- summary(posterior)$summary[, c("Rhat", "n_eff")]
  convergence.out <- c()
  for(i in 1:max(data.list$X)) {
    mu.rhat <- convergence[paste("mu[", i, "]", sep = ''), "Rhat"]
    mu.ess <- convergence[paste("mu[", i, "]", sep = ''), "n_eff"]
    convergence.row <- data.frame(s = data.list$site,
                                  g = data.list$G[data.list$X == i][1],
                                  n = sum(data.list$X == i),
                                  mu.rhat = mu.rhat, mu.ess = mu.ess,
                                  divergence = divergence,
                                  treedepth = treedepth,
                                  stringsAsFactors = FALSE)
    convergence.out <- rbind(convergence.out, convergence.row)
  }


  # posterior data
  posterior <- data.frame(extract(posterior))
  statistics.out <- c()
  for(i in 1:(max(data.list$X) - 1)) {
    for(j in (i + 1):max(data.list$X)) {
      # general data
      site <- data.list$site
      n.i <- data.list$Ng[i]
      n.j <- data.list$Ng[j]
      g.i <- data.list$G[data.list$X == i][1]
      g.j <- data.list$G[data.list$X == j][1]
      mutation <- paste(g.i, "->", g.j, sep = '')
      general <- paste(g.i, ":", n.i, ", ", g.j, ":", n.j, sep = '')


      # extract posterior
      mu.i <- posterior[, paste("mu.", i, sep = '')]
      mu.j <- posterior[, paste("mu.", j, sep = '')]


      # compute Cohen's d and HDI's
      absolute.d <- mu.i - mu.j
      absolute.d.mean <- mean(absolute.d)
      absolute.d.hdi <- getHdi(vec = absolute.d, hdi.level = hdi.level)
      absolute.d.L <- absolute.d.hdi[1]
      absolute.d.H <- absolute.d.hdi[2]
      absolute.d.hdi <- paste("(", round(x = absolute.d.L, digits = 2), ", ",
                              round(x = absolute.d.H, digits = 2), ")",
                              sep = '')


      # Bhat coeff
      ppc.i <- numeric(length = nrow(posterior))
      ppc.j <- numeric(length = nrow(posterior))
      for(p in 1:nrow(posterior)) {
        ppc.i[p] <- mean(stats::rbinom(n = 100, size = 1, prob = mu.i[p]))
        ppc.j[p] <- mean(stats::rbinom(n = 100, size = 1, prob = mu.j[p]))
      }
      b.coef <- getBhattacharyya(x = ppc.i, y = ppc.j)


      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          absolute.d = absolute.d.mean,
                          absolute.d.L = absolute.d.L,
                          absolute.d.H = absolute.d.H,
                          absolute.d.hdi = absolute.d.hdi,
                          b.coef = b.coef)
      statistics.out <- rbind(statistics.out, stats)
    }
  }

  return (list(statistics.out = statistics.out,
               convergence.out = convergence.out))
}




# Description:
# Given a phenotype.type, the procedure compiles the appropriate STAN model.
compileModel <- function(phenotype.type) {
  cat("============================================================= \n")
  cat("===================== Compiling Model ======================= \n")
  cat("============================================================= \n")
  if(phenotype.type == "continuous") {
    # f <- "inst/extdata/continuous.stan"
    f <- system.file("extdata", "continuous.stan", package = "genphen")
    model.stan <- stan_model(file = f, model_name = "continuous")
  }
  else if(phenotype.type == "dichotomous") {
    # f <- "inst/extdata/dichotomous.stan"
    f <- system.file("extdata", "dichotomous.stan", package = "genphen")
    model.stan <- stan_model(file = f, model_name = "dichotomous")
  }

  return(model.stan)
}



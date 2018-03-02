


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
                       stat.learn.method, cv.iterations) {

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
      stop("stat.learn.method must be a string: 'rf', 'svm' or 'none'")
    }

    if(!is.character(stat.learn.method)) {
      stop("stat.learn.method must be a string: 'rf', 'svm' or 'none'")
    }

    if(!stat.learn.method %in% c("rf", "svm", "none")) {
      stop("stat.learn.method must be a string: 'rf', 'svm' or 'none'")
    }
  }

  checkCv <- function(stat.learn.method, cv.iterations) {
    if(stat.learn.method %in% c("rf", "svm")) {
      if(length(cv.iterations) != 1) {
        stop("cv.iterations must be a number (default = 1,000).")
      }

      if(is.numeric(cv.iterations) == FALSE) {
        stop("cv.iterations must be a number (default = 1,000).")
      }

      if(cv.iterations < 100) {
        stop("cv.iterations >= 100 recomended (default = 1,000).")
      }
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
     is.null(cv.iterations) | missing(cv.iterations)) {
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
  checkCv(stat.learn.method = stat.learn.method, cv.iterations = cv.iterations)
}





# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInputDiagnostics <- function(genotype, anchor.points,
                                  with.anchor.points,
                                  rf.importance.trees){



  checkAnchors <- function(anchor.points, rf.importance.trees) {
    if(length(anchor.points) <= 0) {
      stop("anchor.points must be a numeric vector in range [1:ncol(genotype)]")
    }
    if(length(rf.importance.trees) != 1) {
      stop("diagnostics.samples must be a number (default = 50,000).")
    }


    if(is.numeric(anchor.points) == FALSE) {
      stop("anchor.points must be a numeric vector.")
    }
    if(is.numeric(rf.importance.trees) == FALSE) {
      stop("rf.importance.trees must be a number (default = 50,000).")
    }


    if(any(anchor.points <= 0)) {
      stop("anchor.points must be a numeric vector in range [1:ncol(genotype)]")
    }
    if(rf.importance.trees <= 10000) {
      stop("rf.importance.trees >= 10,000 accepted (default = 50,000).")
    }
  }

  checkWithAnchorPoints <- function(with.anchor.points) {
    if(length(with.anchor.points) != 1) {
      stop("with.anchor.points must be a logical parameter")
    }

    if(is.logical(with.anchor.points) == FALSE) {
      stop("with.anchor.points must be a logical parameter")
    }
  }


  if(is.null(anchor.points) | missing(anchor.points) |
     is.null(with.anchor.points) | missing(with.anchor.points) |
     is.null(rf.importance.trees) | missing(rf.importance.trees)) {
    stop("arguments must be non-NULL/specified")
  }

  checkAnchors(anchor.points = anchor.points,
               rf.importance.trees = rf.importance.trees)
  checkWithAnchorPoints(with.anchor.points = with.anchor.points)

  if(all(anchor.points %in% 1:ncol(genotype)) == FALSE) {
    stop("The anchor points must lie in the genotype space:",
         "[", 1, "-", ncol(genotype), "] \n", sep = '')
  }
}







# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using random forests.
getRfCa <- function(data.list, cv.fold, cv.steps,
                    hdi.level, ntree, mcmc.cores) {


  # Description:
  # Perform N number of classifications and compute N number of:
  # - classification accuracy
  # - kappa statistics
  # - number of successful classifications (ideally == N)
  booter <- function(X, Y, cv.fold, cv.steps, ntree) {
    # number of total data entries
    rows <- length(Y)

    # sample at random
    s <- sample(x = 1:rows, size = round(x = cv.fold * rows, digits = 0),
                replace = TRUE)
    train <- data.frame(Y = Y[s], X = X[s], stringsAsFactors = FALSE)
    test <- data.frame(Y = Y[-s], X = X[-s], stringsAsFactors = FALSE)

    # only one type of predictor (no continous variable)
    if(length(unique(train$X)) <= 1) {
      return (list(ca = NA, kappa = NA))
    }
    else {
      # train classification model (try condition to avoid errors in case
      # only one-level predictor is train data)
      rf.out <- try(ranger::ranger(as.factor(Y) ~ X, data = train,
                                   num.trees = ntree), silent = TRUE)
      if(attr(rf.out, "class")[1] == "try-error") {
        return (list(ca = NA, kappa = NA))
      }
      else {
        # test classification model
        pr <- stats::predict(object = rf.out, data = test)

        # compute classification accuracy (1 - classification error)
        ca <- sum(as.character(test$Y)==as.character(pr$predictions))/nrow(test)

        # compute kappa statistics
        kappa <- getKappa(real = as.character(test$Y),
                          predicted = as.character(pr$predictions),
                          aas = unique(Y))
        return (list(ca = ca, kappa = kappa))
      }
    }
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
      registerDoMC(cores = mcmc.cores)
      class.obj <- (foreach(f = 1:cv.steps) %dopar% booter(Y = Y, X = X,
                                                           cv.fold = cv.fold,
                                                           cv.steps = cv.steps,
                                                           ntree = 1000))

      # get cas and kappas
      class.obj <- unlist(class.obj)
      ca <- as.numeric(class.obj[names(class.obj) == "ca"])
      kappa <- as.numeric(class.obj[names(class.obj) == "kappa"])


      # get 95% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = ca, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      ca.hdi <- paste("(", round(x = ca.L, digits = 2), ", ",
                      round(x = ca.H, digits = 2), ")", sep = '')


      # build 95% HDI for the kappas
      kappa.hdi <- getHdi(vec = kappa, hdi.level = hdi.level)
      kappa.L <- as.numeric(kappa.hdi[1])
      kappa.H <- as.numeric(kappa.hdi[2])
      kappa.hdi <- paste("(", round(x = kappa.L, digits = 2), ", ",
                         round(x = kappa.H, digits = 2), ")", sep = '')

      # pack outputs
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          ca = mean(ca, na.rm = TRUE),
                          ca.L = ca.L,
                          ca.H = ca.H,
                          ca.hdi = ca.hdi,
                          kappa = mean(kappa, na.rm = TRUE),
                          kappa.L = kappa.L,
                          kappa.H = kappa.H,
                          kappa.hdi = kappa.hdi)
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
getSvmCa <- function(data.list, cv.fold, cv.steps,
                     hdi.level, mcmc.cores) {


  # Description:
  # Perform N number of classifications and compute N number of:
  # - classification accuracy
  # - kappa statistics
  # - number of successful classifications (ideally == N)
  booter <- function(X, Y, cv.fold) {
    # number of total data entries
    rows <- length(Y)

    # sample at random
    s <- sample(x = 1:rows, size = round(x = cv.fold * rows, digits = 0),
                replace = TRUE)
    train <- data.frame(Y = Y[s], X = X[s], stringsAsFactors = FALSE)
    test <- data.frame(Y = Y[-s], X = X[-s], stringsAsFactors = FALSE)

    # only one type of predictor (no continous variable)
    if(length(unique(train$X)) <= 1) {
      return (list(ca = NA, kappa = NA))
    }
    else {
      # train classification model (try condition to avoid errors in case
      # only one-level predictor is train data)
      svm.out <- try(e1071::svm(as.factor(Y) ~ X, data = train,
                                type = "C-classification"),
                     silent = TRUE)
      if(attr(svm.out, "class")[1] == "try-error") {
        return (list(ca = NA, kappa = NA))
      }
      else {
        # test classification model
        prediction <- stats::predict(object = svm.out, newdata = test)

        # compute classification accuracy (1 - classification error)
        ca <- sum(as.character(test$Y)==as.character(prediction))/nrow(test)

        # compute kappa statistics
        kappa <- getKappa(real = as.character(test$Y),
                          predicted = as.character(prediction),
                          aas = unique(Y))
        return (list(ca = ca, kappa = kappa))
      }
    }
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
      registerDoMC(cores = mcmc.cores)
      class.obj <- (foreach(f = 1:cv.steps) %dopar% booter(Y = Y, X = X,
                                                           cv.fold = cv.fold))

      # get cas and kappas
      class.obj <- unlist(class.obj)
      ca <- as.numeric(class.obj[names(class.obj) == "ca"])
      kappa <- as.numeric(class.obj[names(class.obj) == "kappa"])



      # get HDI for the classification accuracy
      ca.hdi <- getHdi(vec = ca, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      ca.hdi <- paste("(", round(x = ca.L, digits = 2), ", ",
                      round(x = ca.H, digits = 2), ")", sep = '')


      # build 95% HDI for the kappas
      kappa.hdi <- getHdi(vec = kappa, hdi.level = hdi.level)
      kappa.L <- as.numeric(kappa.hdi[1])
      kappa.H <- as.numeric(kappa.hdi[2])
      kappa.hdi <- paste("(", round(x = kappa.L, digits = 2), ", ",
                         round(x = kappa.H, digits = 2), ")", sep = '')

      # pack outputs
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          ca = mean(ca, na.rm = TRUE),
                          ca.L = ca.L,
                          ca.H = ca.H,
                          ca.hdi = ca.hdi,
                          kappa = mean(kappa, na.rm = TRUE),
                          kappa.L = kappa.L,
                          kappa.H = kappa.H,
                          kappa.hdi = kappa.hdi)
      ca.out <- rbind(ca.out, stats)
    }
  }

  return(ca.out)
}




# Description:
# Given a confusion matrix table(predicted, real), compute the Cohen's kappa
# statistics. Cohen makes the following distinction between the different
# kappa ranges:
# if kappa<0 => "no agreement"
# if 0.0-0.2 => "slignt agreement"
# if 0.2-0.4 => "fair agreement"
# if 0.4-0.6 => "moderate agreement"
# if 0.6-0.8 => "substantial agreement"
# if 0.8-1.0 => "almost perfect agreement"
getKappa <- function(predicted, real, aas) {

  buildConfusionMatrix <- function(predicted, real) {
    conf <- matrix(data = 0, nrow = 2, ncol = 2)
    conf[1, 1] <- length(intersect(which(real %in% aas[1]),
                                   which(predicted %in% aas[1])))
    conf[2, 2] <- length(intersect(which(real %in% aas[2]),
                                   which(predicted %in% aas[2])))
    conf[2, 1] <- length(intersect(which(real %in% aas[1]),
                                   which(!predicted %in% aas[1])))
    conf[1, 2] <- length(intersect(which(real %in% aas[2]),
                                   which(!predicted %in% aas[2])))
    return (conf)
  }

  conf <- buildConfusionMatrix(predicted = predicted, real = real)
  expected.accuracy <- (sum(conf[1, ])/sum(conf) * sum(conf[, 1])/sum(conf))+
    (sum(conf[2, ])/sum(conf) * sum(conf[, 2])/sum(conf))
  real.accuracy <- (conf[1, 1] + conf[2, 2])/sum(conf)
  kappa <- (real.accuracy - expected.accuracy)/(1 - expected.accuracy)

  return (kappa)
}





# Description:
# Dummy CA output
getNoneCa <- function(data.list) {


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


      # pack outputs
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          ca = NA,
                          ca.L = NA,
                          ca.H = NA,
                          ca.hdi = NA,
                          ca.boots = NA)
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

  #Calculating the Bhattacharyya Coefficient (sum of the square root of
  # the multiple of the relative counts of both distributions)
  bc <- sum(sqrt(rel.histx*rel.histy))
  b.coef.x <- sum(sqrt((histx[histx != 0]/sum(histx[histx != 0]))*
                         (histy[histx != 0]/sum(histy[histx != 0]))))
  b.coef.y <- sum(sqrt((histx[histy != 0]/sum(histx[histy != 0]))*
                         (histy[histy != 0]/sum(histy[histy != 0]))))
  b.coef.max <- max(b.coef.x, b.coef.y)
  return(list(bc = bc, b.coef.max = b.coef.max))
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
      ppc.i <- mean(mu.i)+mean(sigma.i)*stats::rt(n = 10^6, df = mean(nu))
      ppc.j <- mean(mu.j)+mean(sigma.j)*stats::rt(n = 10^6, df = mean(nu))
      bhat <- getBhattacharyya(x = ppc.i, y = ppc.j)
      bc <- bhat$bc


      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          cohens.d = cohens.d.mean,
                          cohens.d.L = cohens.d.L,
                          cohens.d.H = cohens.d.H,
                          cohens.d.hdi = cohens.d.hdi,
                          bc = bc)
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
      # ppc.i <- numeric(length = nrow(posterior))
      # ppc.j <- numeric(length = nrow(posterior))
      # for(p in 1:nrow(posterior)) {
      #   ppc.i[p] <- mean(stats::rbinom(n = 100, size = 1, prob = mu.i[p]))
      #   ppc.j[p] <- mean(stats::rbinom(n = 100, size = 1, prob = mu.j[p]))
      # }

      # Bhat coeff
      ppc.i <- stats::rbinom(n = 10^6, prob = mean(mu.i), size = n.i)/n.i
      ppc.j <- stats::rbinom(n = 10^6, prob = mean(mu.j), size = n.j)/n.j
      bhat <- getBhattacharyya(x = ppc.i, y = ppc.j)
      bc <- bhat$bc

      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          absolute.d = absolute.d.mean,
                          absolute.d.L = absolute.d.L,
                          absolute.d.H = absolute.d.H,
                          absolute.d.hdi = absolute.d.hdi,
                          bc = bc)
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
    f.local <- "inst/extdata/continuous.stan"
    f.pkg <- system.file("extdata", "continuous.stan", package = "genphen")
  }
  else if(phenotype.type == "dichotomous") {
    f.local <- "inst/extdata/dichotomous.stan"
    f.pkg <- system.file("extdata", "dichotomous.stan", package = "genphen")
  }

  if(file.exists(f.pkg)) {
    model.stan <- stan_model(file = f.pkg, model_name = "model")
  }
  if(file.exists(f.local)) {
    model.stan <- stan_model(file = f.local, model_name = "model")
  }

  return(model.stan)
}



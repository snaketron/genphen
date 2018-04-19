


# Description:
# Get genotype-phenotype data in format for stan
getGenphenData <- function(genotype, phenotype,
                           phenotype.type,
                           min.observations) {

  if(phenotype.type == "continuous") {
    out <- vector(mode = "list")
    out.counter <- 1
    E_mu <- 0
    E_sigma <- 0

    for(i in 1:ncol(genotype)) {
      js <- unique(genotype[, i])
      if(length(js) != 1) {
        X <- as.numeric(as.factor(genotype[, i]))
        k <- which(table(X) >= min.observations)
        if(length(k) != 1) {
          Y <- phenotype[X %in% as.numeric(names(k))]
          G <- genotype[X %in% as.numeric(names(k)), i]
          X <- as.numeric(as.factor(X[X %in% as.numeric(names(k))]))
          Ng <- numeric(length = length(unique(X))) # nr.of samples per genotype
          for(j in 1:max(X)) {
            temp.Y <- Y[X == j]
            Ng[j] <- length(temp.Y)
            E_mu <- c(E_mu, mean(temp.Y))
            E_sigma <- c(E_sigma, stats::sd(temp.Y))
          }

          l <- list(site = i, G = G, Y = Y, X = X, Ng = Ng,
                    Nx = length(Ng), Ny = length(Y))
          out[[out.counter]] <- l
          out.counter <- out.counter + 1
        }
      }
    }

    # if empty return NULL
    if(length(out) == 0) {
      return(NULL)
    }

    # empirical mean and SD
    E_mu <- mean(E_mu, na.rm = TRUE)
    E_sigma <- mean(E_sigma, na.rm = TRUE)
    for(i in 1:length(out)) {
      out[[i]]$E_mu <- E_mu
      out[[i]]$E_sigma <- E_sigma
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

          Ng <- numeric(length = length(unique(X))) # nr.of samples per genotype
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
# @Deprecated
# Get genotype-phenotype data in format for stan
getCompleteGenphenData <- function(genotype,
                                   phenotype,
                                   min.observations) {

  out.data <- c()
  for(i in 1:ncol(genotype)) {
    s.data <- data.frame(S = i, genotype = genotype[, i],
                         Y = phenotype, stringsAsFactors = FALSE)
    out.data <- rbind(out.data, s.data)
  }

  # id variable
  out.data$X <- paste(out.data$S, out.data$genotype, sep = '')

  # as factor -> as numeric
  out.data$X <- as.numeric(as.factor(out.data$X))
  out.data <- out.data[order(out.data$X, decreasing = FALSE), ]

  # I -> counter
  out.data$I <- 1
  check <- aggregate(I~X, data = out.data, FUN = sum)

  # keep only ids with count >= min observations
  out.data <- out.data[out.data$X %in% check$X[check$I >= min.observations], ]

  # II -> remove single sites
  sites.to.remove <- c()
  for(site in unique(out.data$S)) {
    gs <- unique(out.data$genotype[out.data$S == site])
    if(length(gs) == 1) {
      sites.to.remove <- c(sites.to.remove, site)
    }
  }
  out.data <- out.data[!out.data$S %in% sites.to.remove, ]

  # refactor in case some were removed in previous step
  out.data$X <- as.numeric(as.factor(out.data$X))
  out.data$S <- as.numeric(as.factor(out.data$S))

  # compute empirical mean/sd
  stats.mean <- aggregate(Y~X, data = out.data, FUN = mean)
  colnames(stats.mean) <- c("X", "mean.Y")
  stats.sd <- aggregate(Y~X, data = out.data, FUN = sd)
  colnames(stats.sd) <- c("X", "sd.Y")
  stats <- merge(x = stats.mean, y = stats.sd, by = "X")

  # final data for stan
  data.list <- list(S = out.data$S,
                    genotype = out.data$genotype,
                    Y = out.data$Y,
                    X = out.data$X,
                    Ny = length(out.data$Y),
                    Nx = length(unique(out.data$X)),
                    Ns = length(unique(out.data$S)),
                    E_mu = mean(stats$mean.Y),
                    E_sigma = mean(stats$sd.Y))

  return (data.list)
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
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInputPhyloBias <- function(input.kinship.matrix, genotype) {


  checkGenotype <- function(genotype) {
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
  }


  checkKinship <- function(input.kinship.matrix) {

    if(is.matrix(input.kinship.matrix) == FALSE) {
      stop("precomputed kinship matrix must be a numeric matrix.")
    }

    if(is.numeric(input.kinship.matrix) == FALSE) {
      stop("precomputed kinship matrix must be a numeric matrix.")
    }

    if(nrow(input.kinship.matrix) != ncol(input.kinship.matrix)) {
      stop("precomputed kinship matrix must be NxN numeric matrix.")
    }

    if(nrow(input.kinship.matrix) <= 0) {
      stop("at least two individuals needed for the analysis.")
    }
  }


  if((is.null(input.kinship.matrix) | missing(input.kinship.matrix))
     & (is.null(genotype) | missing(genotype))) {
    stop("arguments must be non-NULL/specified")
  }

  if(is.null(input.kinship.matrix) | missing(input.kinship.matrix)) {
    if(is.null(genotype) | missing(genotype)) {
      stop("arguments must be non-NULL/specified")
    }
  }
  else {
    checkKinship(input.kinship.matrix = input.kinship.matrix)
  }
  checkGenotype(genotype = genotype)
}







# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using random forests.
getRfCa <- function(data.list, cv.fold, cv.steps,
                    hdi.level, ntree, mcmc.cores) {


  # Description:
  # Performs the bootstrapping iteratively and breaks if convergence
  # is met before the number of steps is hit.
  getIncrementalLearning <- function(Y, X, cv.fold, ntree, I, e = 0.01) {


    # Description:
    # Perform N number of classifications and compute N number of:
    # - classification accuracy
    # - kappa statistics
    # - number of successful classifications (ideally == N)
    booter <- function(X, Y, cv.fold, ntree) {
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
          ca<-sum(as.character(test$Y)==as.character(pr$predictions))/nrow(test)

          # compute kappa statistics
          kappa <- getKappa(real = as.character(test$Y),
                            predicted = as.character(pr$predictions),
                            aas = unique(Y))
          return (list(ca = ca, kappa = kappa))
        }
      }
    }


    # Description:
    # Given a confusion matrix table(predicted, real), compute the Cohen's
    # kappa statistics. Cohen makes the following distinction between the
    # different kappa ranges.
    getKappa <- function(predicted, real, aas) {

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

      return (kappa)
    }


    # cv.steps == 100
    if(I == 1) {
      # ca bootstrap
      ca.obj <- (foreach(f = 1:100) %dopar% booter(Y = Y, X = X,
                                                   cv.fold = cv.fold,
                                                   ntree = 1000))

      # get cas and kappas
      ca.obj <- unlist(ca.obj)
      cas <- as.numeric(ca.obj[names(ca.obj) == "ca"])
      kappas <- as.numeric(ca.obj[names(ca.obj) == "kappa"])

      # get 95% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = cas, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      ca.hdi <- paste("(", round(x = ca.L, digits = 2), ", ",
                      round(x = ca.H, digits = 2), ")", sep = '')

      # build 95% HDI for the kappas
      kappa.hdi <- getHdi(vec = kappas, hdi.level = hdi.level)
      kappa.L <- as.numeric(kappa.hdi[1])
      kappa.H <- as.numeric(kappa.hdi[2])
      kappa.hdi <- paste("(", round(x = kappa.L, digits = 2), ", ",
                         round(x = kappa.H, digits = 2), ")", sep = '')

      return(list(ca = mean(cas, na.rm = TRUE),
                  ca.L = ca.L,
                  ca.H = ca.H,
                  ca.hdi = ca.hdi,
                  kappa = mean(kappas, na.rm = TRUE),
                  kappa.L = kappa.L,
                  kappa.H = kappa.H,
                  kappa.hdi = kappa.hdi,
                  I = 1))
    }


    # keeptrack
    old.list <- list(kappa.L = NA, kappa.H = NA, ca.L = NA, ca.H = NA)
    updated.list <- list(kappa.L = NA, kappa.H = NA, ca.L = NA, ca.H = NA)

    cas <- c()
    kappas <- c()
    for(i in 1:I) {
      ca.obj <- (foreach(f = 1:100) %dopar% booter(Y = Y, X = X,
                                                   cv.fold = cv.fold,
                                                   ntree = 1000))

      # get cas and kappas
      ca.obj <- unlist(ca.obj)
      new.ca <- as.numeric(ca.obj[names(ca.obj) == "ca"])
      new.kappa <- as.numeric(ca.obj[names(ca.obj) == "kappa"])

      # Update parameters
      cas <- c(cas, new.ca)
      kappas <- c(kappas, new.kappa)

      # get 95% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = cas, hdi.level = hdi.level)
      updated.list[["ca.L"]] <- as.numeric(ca.hdi[1])
      updated.list[["ca.H"]] <- as.numeric(ca.hdi[2])

      # build 95% HDI for the kappas
      kappa.hdi <- getHdi(vec = kappas, hdi.level = hdi.level)
      updated.list[["kappa.L"]] <- as.numeric(kappa.hdi[1])
      updated.list[["kappa.H"]] <- as.numeric(kappa.hdi[2])

      if(i > 1) {
        # Error:
        errors <- c(abs(updated.list[["kappa.L"]]-old.list[["kappa.L"]]) <= e,
                    abs(updated.list[["kappa.H"]]-old.list[["kappa.H"]]) <= e,
                    abs(updated.list[["ca.L"]]-old.list[["ca.L"]]) <= e,
                    abs(updated.list[["ca.H"]]-old.list[["ca.H"]]) <= e)

        if(all(errors) == TRUE) {
          k.l <- round(x = updated.list[["kappa.L"]], digits = 2)
          k.h <- round(x = updated.list[["kappa.H"]], digits = 2)
          c.l <- round(x = updated.list[["ca.L"]], digits = 2)
          c.h <- round(x = updated.list[["ca.H"]], digits = 2)
          return(list(ca = mean(cas, na.rm = TRUE),
                      ca.L = updated.list[["ca.L"]],
                      ca.H = updated.list[["ca.H"]],
                      ca.hdi = paste("(", c.l, ", ", c.h, ")", sep = ''),
                      kappa = mean(kappas, na.rm = TRUE),
                      kappa.L = updated.list[["kappa.L"]],
                      kappa.H = updated.list[["kappa.H"]],
                      kappa.hdi = paste("(", k.l, ", ", k.h, ")", sep = ''),
                      I = i))
        }
      }

      # keep post
      old.list <- updated.list
    }


    # if no speedup is possible, return the result after cv.steps
    k.l <- round(x = updated.list[["kappa.L"]], digits = 2)
    k.h <- round(x = updated.list[["kappa.H"]], digits = 2)
    c.l <- round(x = updated.list[["ca.L"]], digits = 2)
    c.h <- round(x = updated.list[["ca.H"]], digits = 2)
    return(list(ca = mean(cas, na.rm = TRUE),
                ca.L = updated.list[["ca.L"]],
                ca.H = updated.list[["ca.H"]],
                ca.hdi = paste("(", c.l, ", ", c.h, ")", sep = ''),
                kappa = mean(kappas, na.rm = TRUE),
                kappa.L = updated.list[["kappa.L"]],
                kappa.H = updated.list[["kappa.H"]],
                kappa.hdi = paste("(", k.l, ", ", k.h, ")", sep = ''),
                I = I))
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


      # run incremental CA learnings (100 steps in each iteration)
      cl <- parallel::makeCluster(mcmc.cores)
      doParallel::registerDoParallel(cl)
      class.obj <- getIncrementalLearning(Y = Y,
                                          X = X,
                                          cv.fold = cv.fold,
                                          ntree = ntree,
                                          I = ceiling(x = cv.steps/100),
                                          e = 0.01)
      parallel::stopCluster(cl = cl)
      doParallel::stopImplicitCluster()

      # pack outputs
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          ca = class.obj$ca,
                          ca.L = class.obj$ca.L,
                          ca.H = class.obj$ca.H,
                          ca.hdi = class.obj$ca.hdi,
                          kappa = class.obj$kappa,
                          kappa.L = class.obj$kappa.L,
                          kappa.H = class.obj$kappa.H,
                          kappa.hdi = class.obj$kappa.hdi,
                          I = class.obj$I)
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
  # Performs the bootstrapping iteratively and breaks if convergence is met
  # before the number of steps is hit.
  getIncrementalLearning <- function(Y, X, cv.fold, I, e = 0.01) {

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

    # Description:
    # Given a confusion matrix table(predicted, real), compute the Cohen's
    # kappa statistics. Cohen makes the following distinction between the
    # different kappa ranges.
    getKappa <- function(predicted, real, aas) {

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

      return (kappa)
    }


    # cv.steps == 100
    if(I == 1) {
      # ca bootstrap
      ca.obj <- (foreach(f = 1:100) %dopar% booter(Y = Y, X = X,
                                                   cv.fold = cv.fold))

      # get cas and kappas
      ca.obj <- unlist(ca.obj)
      cas <- as.numeric(ca.obj[names(ca.obj) == "ca"])
      kappas <- as.numeric(ca.obj[names(ca.obj) == "kappa"])

      # get 95% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = cas, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      ca.hdi <- paste("(", round(x = ca.L, digits = 2), ", ",
                      round(x = ca.H, digits = 2), ")", sep = '')

      # build 95% HDI for the kappas
      kappa.hdi <- getHdi(vec = kappas, hdi.level = hdi.level)
      kappa.L <- as.numeric(kappa.hdi[1])
      kappa.H <- as.numeric(kappa.hdi[2])
      kappa.hdi <- paste("(", round(x = kappa.L, digits = 2), ", ",
                         round(x = kappa.H, digits = 2), ")", sep = '')

      return(list(ca = mean(cas, na.rm = TRUE),
                  ca.L = ca.L,
                  ca.H = ca.H,
                  ca.hdi = ca.hdi,
                  kappa = mean(kappas, na.rm = TRUE),
                  kappa.L = kappa.L,
                  kappa.H = kappa.H,
                  kappa.hdi = kappa.hdi,
                  I = 1))
    }


    # keeptrack
    old.list <- list(kappa.L = NA, kappa.H = NA, ca.L = NA, ca.H = NA)
    updated.list <- list(kappa.L = NA, kappa.H = NA, ca.L = NA, ca.H = NA)

    cas <- c()
    kappas <- c()
    for(i in 1:I) {
      ca.obj <- (foreach(f = 1:100) %dopar% booter(Y = Y, X = X,
                                                   cv.fold = cv.fold))

      # get cas and kappas
      ca.obj <- unlist(ca.obj)
      new.ca <- as.numeric(ca.obj[names(ca.obj) == "ca"])
      new.kappa <- as.numeric(ca.obj[names(ca.obj) == "kappa"])

      # Update parameters
      cas <- c(cas, new.ca)
      kappas <- c(kappas, new.kappa)

      # get 95% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = cas, hdi.level = hdi.level)
      updated.list[["ca.L"]] <- as.numeric(ca.hdi[1])
      updated.list[["ca.H"]] <- as.numeric(ca.hdi[2])

      # build 95% HDI for the kappas
      kappa.hdi <- getHdi(vec = kappas, hdi.level = hdi.level)
      updated.list[["kappa.L"]] <- as.numeric(kappa.hdi[1])
      updated.list[["kappa.H"]] <- as.numeric(kappa.hdi[2])

      if(i > 1) {
        # Error:
        errors <- c(abs(updated.list[["kappa.L"]]-old.list[["kappa.L"]]) <= e,
                    abs(updated.list[["kappa.H"]]-old.list[["kappa.H"]]) <= e,
                    abs(updated.list[["ca.L"]]-old.list[["ca.L"]]) <= e,
                    abs(updated.list[["ca.H"]]-old.list[["ca.H"]]) <= e)

        if(all(errors) == TRUE) {
          k.l <- round(x = updated.list[["kappa.L"]], digits = 2)
          k.h <- round(x = updated.list[["kappa.H"]], digits = 2)
          c.l <- round(x = updated.list[["ca.L"]], digits = 2)
          c.h <- round(x = updated.list[["ca.H"]], digits = 2)
          return(list(ca = mean(cas, na.rm = TRUE),
                      ca.L = updated.list[["ca.L"]],
                      ca.H = updated.list[["ca.H"]],
                      ca.hdi = paste("(", c.l, ", ", c.h, ")", sep = ''),
                      kappa = mean(kappas, na.rm = TRUE),
                      kappa.L = updated.list[["kappa.L"]],
                      kappa.H = updated.list[["kappa.H"]],
                      kappa.hdi = paste("(", k.l, ", ", k.h, ")", sep = ''),
                      I = i))
        }
      }

      # keep post
      old.list <- updated.list
    }


    # if no speedup is possible, return the result after cv.steps
    k.l <- round(x = updated.list[["kappa.L"]], digits = 2)
    k.h <- round(x = updated.list[["kappa.H"]], digits = 2)
    c.l <- round(x = updated.list[["ca.L"]], digits = 2)
    c.h <- round(x = updated.list[["ca.H"]], digits = 2)
    return(list(ca = mean(cas, na.rm = TRUE),
                ca.L = updated.list[["ca.L"]],
                ca.H = updated.list[["ca.H"]],
                ca.hdi = paste("(", c.l, ", ", c.h, ")", sep = ''),
                kappa = mean(kappas, na.rm = TRUE),
                kappa.L = updated.list[["kappa.L"]],
                kappa.H = updated.list[["kappa.H"]],
                kappa.hdi = paste("(", k.l, ", ", k.h, ")", sep = ''),
                I = I))
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




      # run incremental CA learnings (100 steps in each iteration)
      cl <- parallel::makeCluster(mcmc.cores)
      doParallel::registerDoParallel(cl)
      class.obj <- getIncrementalLearning(Y = Y,
                                          X = X,
                                          cv.fold = cv.fold,
                                          I = ceiling(x = cv.steps/100),
                                          e = 0.01)
      parallel::stopCluster(cl = cl)
      doParallel::stopImplicitCluster()


      # pack outputs
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          ca = class.obj$ca,
                          ca.L = class.obj$ca.L,
                          ca.H = class.obj$ca.H,
                          ca.hdi = class.obj$ca.hdi,
                          kappa = class.obj$kappa,
                          kappa.L = class.obj$kappa.L,
                          kappa.H = class.obj$kappa.H,
                          kappa.hdi = class.obj$kappa.hdi,
                          I = class.obj$I)
      ca.out <- rbind(ca.out, stats)
    }
  }

  return(ca.out)
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
runContinuous <- function(data.list, mcmc.chains, mcmc.iterations,
                          mcmc.warmup, mcmc.cores, hdi.level, model.stan,
                          with.rpa = FALSE, rpa.iterations = 0,
                          rpa.rope = 0) {

  # get initial parameter values
  posterior <- sampling(object = model.stan,
                        data = data.list,
                        pars = c("mu", "sigma", "nu"),
                        iter = mcmc.iterations,
                        warmup = mcmc.warmup,
                        chains = mcmc.chains,
                        cores = mcmc.cores,
                        control = list(adapt_delta = 0.95,
                                       max_treedepth = 10),
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
  ppc.out <- c()
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


      # Alternative Cohen's d based on Hedges 1981
      # pool.sd <- sqrt(((sigma.i^2)*(n.i-1)+(sigma.j^2)*(n.j-1))/(n.i+n.j-2))
      # cohens.d <- (mu.i - mu.j)/pool.sd
      cohens.d <- (mu.i - mu.j)/sqrt((sigma.i^2 + sigma.j^2)/2)
      cohens.d.mean <- mean(cohens.d)
      cohens.d.hdi <- getHdi(vec = cohens.d, hdi.level = hdi.level)
      cohens.d.L = cohens.d.hdi[1]
      cohens.d.H = cohens.d.hdi[2]
      cohens.d.hdi <- paste("(", round(x = cohens.d.L, digits = 2), ", ",
                            round(x = cohens.d.H, digits = 2), ")", sep = '')


      # Difference in sd
      sd.d <- sigma.i - sigma.j
      sd.d.mean <- mean(sd.d)
      sd.d.hdi <- getHdi(vec = sd.d, hdi.level = hdi.level)
      sd.d.L = sd.d.hdi[1]
      sd.d.H = sd.d.hdi[2]
      sd.d.hdi <- paste("(", round(x = sd.d.L, digits = 2), ", ",
                        round(x = sd.d.H, digits = 2), ")", sep = '')


      # Bhat coeff
      ppc.i <- mean(mu.i)+mean(sigma.i)*stats::rt(n = 10^6, df = mean(nu))
      ppc.j <- mean(mu.j)+mean(sigma.j)*stats::rt(n = 10^6, df = mean(nu))
      bhat <- getBhattacharyya(x = ppc.i, y = ppc.j)
      bc <- bhat$bc


      # predicted vs real means
      ppc.row <- data.frame(site = site,
                            general = general,
                            predicted.mu.i = mean(ppc.i),
                            predicted.mu.j = mean(ppc.j),
                            real.mu.i = mean(data.list$Y[data.list$X == i]),
                            real.mu.j = mean(data.list$Y[data.list$X == j]))
      ppc.out <- rbind(ppc.out, ppc.row)


      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          cohens.d = cohens.d.mean,
                          cohens.d.L = cohens.d.L,
                          cohens.d.H = cohens.d.H,
                          cohens.d.hdi = cohens.d.hdi,
                          bc = bc,
                          sd.d = sd.d.mean,
                          sd.d.L = sd.d.L,
                          sd.d.H = sd.d.H,
                          sd.d.hdi = sd.d.hdi)
      statistics.out <- rbind(statistics.out, stats)
    }
  }

  # special case for RPA
  rpa.out <- NULL
  if(with.rpa == TRUE) {
    rpa.out <- getRpaContinuous(data.list = data.list,
                                hdi.level = hdi.level,
                                rpa.iterations = rpa.iterations,
                                rpa.rope = rpa.rope,
                                posterior = posterior,
                                model.stan = model.stan,
                                mcmc.warmup = mcmc.warmup,
                                mcmc.iterations = mcmc.iterations,
                                mcmc.chains = mcmc.chains,
                                mcmc.cores = mcmc.cores)
  }

  return (list(statistics.out = statistics.out,
               convergence.out = convergence.out,
               rpa.out = rpa.out,
               ppc.out = ppc.out))
}




# Description:
# Computes a Bayesian odds-ratio test
runDichotomous <- function(data.list, mcmc.chains, mcmc.iterations,
                           mcmc.warmup, mcmc.cores, hdi.level, model.stan,
                           with.rpa = FALSE, rpa.iterations = 0,
                           rpa.rope = 0) {


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
  ppc.out <- c()
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
      ppc.i <- stats::rbinom(n = 10^6, prob = mean(mu.i), size = n.i)/n.i
      ppc.j <- stats::rbinom(n = 10^6, prob = mean(mu.j), size = n.j)/n.j
      bhat <- getBhattacharyya(x = ppc.i, y = ppc.j)
      bc <- bhat$bc


      # predicted vs real means
      ppc.row <- data.frame(site = site,
                            general = general,
                            predicted.mu.i = mean(ppc.i),
                            predicted.mu.j = mean(ppc.j),
                            real.mu.i = mean(data.list$Y[data.list$X == i]),
                            real.mu.j = mean(data.list$Y[data.list$X == j]))
      ppc.out <- rbind(ppc.out, ppc.row)


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

  # special case for RPA
  rpa.out <- NULL
  if(with.rpa == TRUE) {
    rpa.out <- getRpaDichotomous(data.list = data.list,
                                 hdi.level = hdi.level,
                                 rpa.iterations = rpa.iterations,
                                 rpa.rope = rpa.rope,
                                 posterior = posterior,
                                 model.stan = model.stan,
                                 mcmc.warmup = mcmc.warmup,
                                 mcmc.iterations = mcmc.iterations,
                                 mcmc.chains = mcmc.chains,
                                 mcmc.cores = mcmc.cores)
  }

  return (list(statistics.out = statistics.out,
               convergence.out = convergence.out,
               rpa.out = rpa.out,
               ppc.out = ppc.out))
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






# Description:
# Given a genotype dataset containing SNPs (columns) and N individuals (rows),
# the procedure computes a NxN kinship matrix for the individuals and estimates
# the phylogenetic bias related to each SNP.
getPhyloBias <- function(genotype, k.matrix) {
  phylo.bias <- c()

  # total mean phylogenetic distance
  mean.d.t <- mean(k.matrix[upper.tri(x = k.matrix, diag = FALSE)])

  for(i in 1:ncol(genotype)) {
    gs <- unique(genotype[, i])
    for(g in gs) {
      # feature mean phylogenetic distance
      mean.d.f <- mean(k.matrix[genotype[, i] == g, genotype[, i] == g])

      row <- data.frame(site = i, genotype = g, feature.dist = mean.d.f,
                        total.dist = mean.d.t, stringsAsFactors = FALSE)
      phylo.bias <- rbind(phylo.bias, row)
    }
  }

  return(phylo.bias)
}





# Description:
# Given an estimate mean (M) and HDI with low (L) and high (H) intervals, as
# well as a ROPE interval, check if the estimated interval overlaps with the
# ROPE interval.
getRopeTest <- function(ROPE, M, L, H) {
  # min ROPE = 0
  o <- ifelse(test = L <= 0 & H >= 0, yes = "fail", no = "pass")

  # for ROPE > 0
  if(ROPE > 0 & o == "pass") {
    abs.L <- ifelse(test = M > 0, yes = L, no = abs(H))
    o <- ifelse(test = abs.L <= ROPE, yes = "fail", no = "pass")
  }

  return (o)
}





# Description:
# RPA for continuous data
getRpaContinuous <- function(data.list, hdi.level, rpa.iterations,
                             rpa.rope, posterior, model.stan,
                             mcmc.iterations, mcmc.warmup,
                             mcmc.chains, mcmc.cores) {

  # get subset of posterior
  rpa.i <- sample(x = 1:nrow(posterior), size = rpa.iterations, replace = TRUE)
  posterior <- posterior[rpa.i, ]

  # posterior data
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
      rpa.counter <- 0

      for(k in 1:nrow(posterior)) {
        mu.i <- posterior[k, paste("mu.", i, sep = '')]
        sigma.i <- posterior[k, paste("sigma.", i, sep = '')]
        mu.j <- posterior[k, paste("mu.", j, sep = '')]
        sigma.j <- posterior[k, paste("sigma.", j, sep = '')]
        nu <- posterior[k, "nu"]

        # simulate new data
        ppc.i <- mu.i+sigma.i*stats::rt(n = n.i, df = nu)
        ppc.j <- mu.j+sigma.j*stats::rt(n = n.j, df = nu)


        # new data for stan
        phenotype <- c(ppc.i, ppc.j)
        genotype <- matrix(c(rep(x = g.i, times = length(ppc.i)),
                             rep(x = g.j, times = length(ppc.j))),
                           ncol = 1)
        ppc.data.list <- getGenphenData(genotype = genotype,
                                        phenotype = phenotype,
                                        phenotype.type = "continuous",
                                        min.observations = 3)

        # get initial parameter values
        ppc.posterior <- sampling(object = model.stan,
                                  data = ppc.data.list[[1]],
                                  pars = c("mu", "sigma", "nu"),
                                  iter = mcmc.iterations,
                                  warmup = mcmc.warmup,
                                  chains = mcmc.chains,
                                  cores = mcmc.cores,
                                  control = list(adapt_delta = 0.95,
                                                 max_treedepth = 10),
                                  verbose = FALSE,
                                  refresh = -1)

        # extract posterior
        ppc.posterior <- data.frame(extract(ppc.posterior))
        ppc.mu.i <- ppc.posterior[, "mu.1"]
        ppc.sigma.i <- ppc.posterior[, "sigma.1"]
        ppc.mu.j <- ppc.posterior[, "mu.2"]
        ppc.sigma.j <- ppc.posterior[, "sigma.2"]

        # cohen's d
        cohens.d <- (ppc.mu.i-ppc.mu.j)/sqrt((ppc.sigma.i^2+ppc.sigma.j^2)/2)
        M <- mean(cohens.d)
        cohens.d.hdi <- getHdi(vec = cohens.d, hdi.level = hdi.level)
        L = cohens.d.hdi[1]
        H = cohens.d.hdi[2]

        # compute if significant -> increment counter
        rpa <- getRopeTest(ROPE = rpa.rope, M = M, L = L, H = H)
        rpa.counter <- rpa.counter + sum(rpa == "pass")
      }


      # collect
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          rpa.counter = rpa.counter)
      statistics.out <- rbind(statistics.out, stats)
    }
  }

  return (statistics.out)
}




# Description:
# RPA for dichotomous data
getRpaDichotomous <- function(data.list, hdi.level, rpa.iterations,
                              rpa.rope, posterior, model.stan,
                              mcmc.iterations, mcmc.warmup,
                              mcmc.chains, mcmc.cores) {

  # get subset of posterior
  rpa.i <- sample(x = 1:nrow(posterior), size = rpa.iterations, replace = TRUE)
  posterior <- posterior[rpa.i, ]

  # posterior data
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
      rpa.counter <- 0

      for(k in 1:nrow(posterior)) {
        mu.i <- posterior[k, paste("mu.", i, sep = '')]
        mu.j <- posterior[k, paste("mu.", j, sep = '')]

        # simulate new data
        ppc.i <- stats::rbinom(n = n.i, prob = mu.i, size = 1)
        ppc.j <- stats::rbinom(n = n.j, prob = mu.j, size = 1)

        # new data for stan
        phenotype <- c(ppc.i, ppc.j)
        genotype <- matrix(c(rep(x = g.i, times = length(ppc.i)),
                             rep(x = g.j, times = length(ppc.j))),
                           ncol = 1)
        ppc.data.list <- getGenphenData(genotype = genotype,
                                        phenotype = phenotype,
                                        phenotype.type = "dichotomous",
                                        min.observations = 3)

        # get initial parameter values
        ppc.posterior <- sampling(object = model.stan,
                                  data = ppc.data.list[[1]],
                                  pars = c("mu"),
                                  iter = mcmc.iterations,
                                  warmup = mcmc.warmup,
                                  chains = mcmc.chains,
                                  cores = mcmc.cores,
                                  control = list(adapt_delta = 0.95,
                                                 max_treedepth = 10),
                                  verbose = FALSE,
                                  refresh = -1)

        # extract posterior
        ppc.posterior <- data.frame(extract(ppc.posterior))
        ppc.mu.i <- ppc.posterior[, "mu.1"]
        ppc.mu.j <- ppc.posterior[, "mu.2"]

        # absolute d
        absolute.d <- ppc.mu.i - ppc.mu.j
        absolute.d.hdi <- getHdi(vec = absolute.d, hdi.level = hdi.level)
        M <- mean(absolute.d)
        L <- absolute.d.hdi[1]
        H <- absolute.d.hdi[2]

        # compute if significant -> increment counter
        rpa <- getRopeTest(ROPE = rpa.rope, M = M, L = L, H = H)
        rpa.counter <- rpa.counter + sum(rpa == "pass")
      }


      # collect
      stats <- data.frame(site = site,
                          general = general,
                          mutation = mutation,
                          rpa.counter = rpa.counter)
      statistics.out <- rbind(statistics.out, stats)
    }
  }

  return (statistics.out)
}



# Description:
# Posterior predictive check
getPpc <- function() {

}




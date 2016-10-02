# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using random forests.
getRfClassification <- function(entry, cv.fold, cv.steps, hdi.level, ntree) {


  # Description:
  # Perform N number of classifications and compute N number of:
  # - classification accuracy
  # - kappa statistics
  # - number of successful classifications (ideally == N)
  booter <- function(entry, cv.fold, cv.steps, ntree) {
    # number of total data entries
    rows <- nrow(entry)

    # initialize result vectors
    accuracies <- c()
    kappas <- c()
    successful.boots <- 0

    for(b in 1:cv.steps) {
      # sample at random
      s <- sample(x = 1:rows, size = round(x = cv.fold * rows, digits = 0),
                  replace = FALSE)
      train <- entry[s, ]
      test <- entry[-s, ]

      # refactor (just in case)
      train$g <- as.factor(as.character(train$g))
      test$g <- as.factor(as.character(test$g))

      # only one type of predictor (no continous variable)
      if(length(unique(train$p)) <= 1) {
        accuracies <- c(accuracies, NA)
        kappas <- c(kappas, NA)
      }
      else {
        #train classification model (try condition to avoid errors in case
        # only one-level predictor is train data)
        rf.result <- try(randomForest::randomForest(g ~ p,
                                                    data = train,
                                                    ntree = ntree),
                         silent = TRUE)
        if(attr(rf.result, "class")[1] == "try-error") {
          accuracies <- c(accuracies, NA)
          kappas <- c(kappas, NA)
        }
        else {
          # test classification model
          prediction <- stats::predict(object = rf.result, newdata = test)

          # compute classification accuracy (1 - classification error)
          ac <- length(
            which(as.character(test$g)==as.character(prediction)))/nrow(test)
          accuracies <- c(accuracies, ac)


          # compute kappa statistics
          kappa <- getCohensKappa(real = as.character(test$g),
                                  predicted = as.character(prediction),
                                  aas = unique(entry$g))
          kappas <- c(kappas, kappa)


          # record bootstrap
          successful.boots <- successful.boots + 1
        }
      }
    }

    # build result list and return
    result <- list(accuracies = accuracies, kappas = kappas,
                   successful.boots = successful.boots)
    return (result)
  }

  # perform classification
  class.obj <- booter(entry = entry, cv.fold = cv.fold,
                      cv.steps = cv.steps, ntree = ntree)


  # get 95% HDI for the classification accuracy
  ca.hdi <- getHdi(vec = class.obj$accuracies, hdi.level = hdi.level)
  ca.hdi.low <- as.numeric(ca.hdi[1])
  ca.hdi.high <- as.numeric(ca.hdi[2])


  # build 95% HDI for the kappas
  kappa.hdi <- getHdi(vec = class.obj$kappas, hdi.level = hdi.level)
  kappa.hdi.low <- as.numeric(kappa.hdi[1])
  kappa.hdi.high <- as.numeric(kappa.hdi[2])


  # pack results in list
  result <- list(ca = mean(class.obj$accuracies, na.rm = TRUE),
                 ca.hdi.low = ca.hdi.low,
                 ca.hdi.high = ca.hdi.high,
                 ca.boots = class.obj$successful.boots,
                 kappa = mean(class.obj$kappas, na.rm = TRUE),
                 kappa.hdi.low = kappa.hdi.low,
                 kappa.hdi.high = kappa.hdi.high)


  # if problem return NAs
  if(length(result) != 7) {
    result <- list(ca = NA, ca.ci.low = NA, ca.ci.high = NA, ca.boots = NA,
                   kappa = NA, kappa.ci.low = NA, kappa.ci.high = NA)
  }
  return(result)
}







# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using linear SVM.
getSvmClassification <- function(entry, cv.fold, cv.steps, hdi.level) {


  # Description:
  # Perform N number of classifications and compute N number of:
  # - classification accuracies
  # - kappa statistics
  # - number of successful classifications (ideally == N)
  booter <- function(entry, cv.fold, cv.steps) {
    # number of total data entries
    rows <- nrow(entry)

    # initialize result vectors
    accuracies <- c()
    kappas <- c()
    successful.boots <- 0

    for(b in 1:cv.steps) {
      # sample at random
      s <- sample(x = 1:rows, size = round(x = cv.fold * rows, digits = 0),
                  replace = FALSE)
      train <- entry[s, ]
      test <- entry[-s, ]

      # refactor (just in case)
      train$g <- as.factor(as.character(train$g))
      test$g <- as.factor(as.character(test$g))

      # only one type of predictor (no continous variable)
      if(length(unique(train$p)) <= 1) {
        accuracies <- c(accuracies, NA)
        kappas <- c(kappas, NA)
      }
      else {
        #train classification model (try condition to avoid errors in case
        #only one-level predictor is train data)
        svm.result <- try(e1071::svm(g ~ p, data = train, kernel = "linear",
                                     type = "C-classification"),
                          silent = TRUE)
        if(attr(svm.result, "class")[1] == "try-error") {
          accuracies <- c(accuracies, NA)
          kappas <- c(kappas, NA)
        }
        else {
          # test classification model
          prediction <- stats::predict(object = svm.result, newdata = test)

          # compute classification accuracy (1 - classification accuracy)
          ac <- length(
            which(as.character(test$g)==as.character(prediction)))/nrow(test)
          accuracies <- c(accuracies, ac)


          # compute kappa statistics
          kappa <- getCohensKappa(real = as.character(test$g),
                                  predicted = as.character(prediction),
                                  aas = unique(entry$g))
          kappas <- c(kappas, kappa)


          # record bootstrap
          successful.boots <- successful.boots + 1
        }
      }
    }

    # build result list and return
    result <- list(accuracies = accuracies, kappas = kappas,
                   successful.boots = successful.boots)
    return (result)
  }


  # perform classification
  class.obj <- booter(entry = entry, cv.fold = cv.fold, cv.steps = cv.steps)


  # get 95% HDI for the classification accuracy
  ca.hdi <- getHdi(vec = class.obj$accuracies, hdi.level = hdi.level)
  ca.hdi.low <- as.numeric(ca.hdi[1])
  ca.hdi.high <- as.numeric(ca.hdi[2])


  # build 95% HDI for the kappas
  kappa.hdi <- getHdi(vec = class.obj$kappas, hdi.level = hdi.level)
  kappa.hdi.low <- as.numeric(kappa.hdi[1])
  kappa.hdi.high <- as.numeric(kappa.hdi[2])


  # pack results in list
  result <- list(ca = mean(class.obj$accuracies, na.rm = TRUE),
                 ca.hdi.low = ca.hdi.low,
                 ca.hdi.high = ca.hdi.high,
                 ca.boots = class.obj$successful.boots,
                 kappa = mean(class.obj$kappas, na.rm = TRUE),
                 kappa.hdi.low = kappa.hdi.low,
                 kappa.hdi.high = kappa.hdi.high)


  # if problem return NAs
  if(length(result) != 7) {
    result <- list(ca = NA, ca.ci.low = NA, ca.ci.high = NA, ca.boots = NA,
                   kappa = NA, kappa.ci.low = NA, kappa.ci.high = NA)
  }
  return(result)
}






# Description:
# Computes Bayesian oneway ANOVA between two vectors (genotype, phenotype)
# Return: contrasts between categories
getBayesianTtest <- function(entry, n.iter, n.chains, hdi.levels, model) {

  # Description:
  # Get data.list from entry
  getDataList <- function(entry) {
    xs <- as.character(entry$g)
    xs[xs == unique(entry$g)[1]] <- 1
    xs[xs == unique(entry$g)[2]] <- 2
    xs <- as.numeric(xs)
    ys <- entry$p

    mean.1 <- mean(ys)
    mean.2 <- mean(ys)
    sd.1 <- sd(ys)
    sd.2 <- sd(ys)
    # sd.1 <- sd.1 + 0.01
    # sd.2 <- sd.2 + 0.01

    data.list = list(x = xs,
                     y = ys,
                     Ny = length(ys),
                     meanY = c(mean.1, mean.2),
                     sdY = c(sd.1, sd.2))
    return(data.list)
  }


  # Description:
  # Computes the contrasts using a posterior data
  getContrast <- function(posterior, hdi.level) {
    # compute i
    post.i <- posterior[, "mu.1."]
    mu.i <- mean(post.i)
    hdi.i <- getHdi(vec = post.i, hdi.level = hdi.level)
    post.i.sd <- posterior[, "sigma.1."]
    sd.i <- mean(post.i.sd)

    # compute j
    post.j <- posterior[, "mu.2."]
    mu.j <- mean(post.j)
    hdi.j <- getHdi(vec = post.j, hdi.level = hdi.level)
    post.j.sd <- posterior[, "sigma.2."]
    sd.j <- mean(post.j.sd)

    # contract
    mu.diff <- mu.i - mu.j
    hdi.diff <- getHdi(vec = post.i - post.j, hdi.level = hdi.level)

    # pack as data frame
    result <- list(mu.i = mu.i, sd.i = sd.i,
                   mu.j = mu.j, sd.j = sd.j,
                   mu.diff = mu.diff,
                   hdi.diff.L = hdi.diff[1],
                   hdi.diff.H = hdi.diff[2])
    return (result)
  }


  # Description:
  # Computes the effective sampling size
  getEss <- function(posterior) {
    ess.mu.i <- as.numeric(effectiveSize(x = posterior$mu.1))
    ess.mu.j <- as.numeric(effectiveSize(x = posterior$mu.2))
    ess.sigma.i <- as.numeric(effectiveSize(x = posterior$sigma.1))
    ess.sigma.j <- as.numeric(effectiveSize(x = posterior$sigma.2))

    return(list(ess.mu.i = ess.mu.i, ess.mu.j = ess.mu.j,
                ess.sigma.i = ess.sigma.i, ess.sigma.j = ess.sigma.j))
  }


  # get model
  if(model == "tdist") {
    model.file <- system.file("extdata", "tdist.txt", package = "genphen")
  }
  else {
    model.file <- system.file("extdata", "ndist.txt", package = "genphen")
  }

  # TODO: depending on model get different data.list
  # TODO: depending on model sample different posterior variables
  # TODO: depending on model compute specific HDIs
  # TODO: what about multiple chains (average results)


  # data
  data.list <- getDataList(entry = entry)


  jags.model = suppressMessages(expr = jags.model(file = model.file,
                                                  data = data.list,
                                                  n.chains = n.chains,
                                                  n.adapt = n.iter,
                                                  quiet = TRUE))

  # feel d' burn
  suppressMessages(expr = update(jags.model, n.iter = n.iter))


  # sample posterior
  vars <- c("mu", "sigma")
  posterior <- suppressMessages(expr = coda.samples(model = jags.model,
                                                    thin = 1,
                                                    variable.names = vars,
                                                    n.iter = n.iter))
  posterior <- data.frame(do.call(rbind, posterior))


  # compute contrasts
  result <- vector(mode = "list", length = length(hdi.levels))
  names(result) <- hdi.levels
  for(hdi.level in hdi.levels) {
    result[[as.character(hdi.level)]] <- getContrast(posterior = posterior,
                                                     hdi.level = hdi.level)
  }

  # compute ess stats
  ess <- getEss(posterior = posterior)

  return(list(result = result, ess = ess))
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
getCohensKappa <- function(predicted, real, aas) {

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
# Given a vector of genotypes, the allele frequencies are computed.
# Return: 4-element vector (aa1, aa2, aa1.count, aa2.count)
getGeneralStats <- function(entry) {
  aa1 <- as.character(unique(entry$g)[1])
  aa2 <- as.character(unique(entry$g)[2])

  aa1.nr <- length(which(entry$g %in% aa1))
  aa2.nr <- length(which(entry$g %in% aa2))

  result <- list(aa1 = aa1, aa2 = aa2, aa1.nr = aa1.nr, aa2.nr = aa2.nr)
  return (result)
}






# Description:
# Using the effsize package, it computes the phenotype effect size between
# any two categories at a site.
# Return: 3-element vector (effect, ci.high, ci.low)
getCohensD <- function(entry, conf.level) {
  effect.obj <- effsize::cohen.d.formula(formula = p~g, data = entry,
                                         conf.level = conf.level)
  effect <- as.numeric(effect.obj$estimate)
  effect.ci.low <- as.numeric(effect.obj$conf.int[1])
  effect.ci.high <- as.numeric(effect.obj$conf.int[2])

  result <- list(effect = effect, effect.ci.low = effect.ci.low,
                 effect.ci.high = effect.ci.high)
  return (result)
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
# Computes the two-sample t-test for two vectors (genotype, phenotype)
# Return: p-value
getTTestScore <- function(entry) {
  ttest.summary <- t.test(formula = entry$p ~ entry$g)
  return(ttest.summary$p.value)
}





# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInputRf <- function(genotype, phenotype, cv.fold,
                         cv.steps, hdi.level, ntree) {

  if(is.null(genotype) | is.null(phenotype) | is.null(ntree)
     | is.null(cv.fold) | is.null(cv.steps) | is.null(hdi.level)) {
    stop("arguments must be non-NULL")
  }

  # is it matrix/dataframe or DAMultipleAlignment/AAMultipleAlignment
  if(is.null(attr(genotype, "class")) == FALSE) {
    if(!attr(genotype, "class") %in% c("AAMultipleAlignment",
                                       "DNAMultipleAlignment")) {
      stop("genotype can be one of the following structures: matrix/data.frame
           or AAMultipleAlignment/AAMultipleAlignment")
    }
    else {
      temp <- as.matrix(genotype)
      if(nrow(temp) < 2 | ncol(temp) == 0) {
        stop("the genotypes cannot have less than two observations/the number
             of genotypes cannot be 0.")
      }
    }
  }
  else { #matrix/data.frame/vector
    if(is.vector(genotype)) {
      genotype <- matrix(data = genotype, ncol = 1)
    }

    if(nrow(genotype) < 2 | ncol(genotype) == 0) {
      stop("the genotypes cannot have less than two observations/the number
             of genotypes cannot be 0.")
    }

    if(!is.vector(genotype) & !is.matrix(genotype) & !is.data.frame(genotype)) {
      stop("the genotype can be a vector, matrix, data.frame or
           DNAMultipleAlignment/AAMultipleAlignment object")
    }

    if(typeof(genotype) != "character") {
      stop("the genotype entries of character type.")
    }
  }

  if(!is.vector(phenotype)) {
    stop("the phenotype must be a vector.")
  }

  if(length(phenotype) < 2) {
    stop("phenotype cannot be a vector of less than 2 elements.")
  }

  if(!is.numeric(phenotype)) {
    stop("the phenotype must be of numeric type.")
  }

  if(nrow(genotype) != length(phenotype)) {
    stop("the lengths of the genotypes and the phenotypes are not equal. For
         each genotype there should be one phenotype.")
  }

  if(!is.numeric(cv.fold)) {
    stop("cv.fold must be numeric argument in range (0, 1)(default = 0.66).")
  }

  if(cv.fold >= 1 | cv.fold <= 0) {
    stop("cv.fold must be numeric argument in range (0, 1)(default = 0.66).")
  }

  if(!is.numeric(cv.steps)) {
    stop("cv.steps must be a numeric argument larger than 0 (default = 100).")
  }

  if(cv.steps == 0) {
    stop("cv.steps must be larger than 0 (default = 100).")
  }

  if(length(hdi.level) != 1) {
    stop("hdi.level is a numeric argument in range (0, 1)(default = 0.95).")
  }

  if(hdi.level >= 1 | hdi.level <= 0) {
    stop("hdi.level is a numeric argument in range (0, 1)(default = 0.95).")
  }

  if(!is.numeric(hdi.level)) {
    stop("hdi.level is a numeric argument in range (0, 1)(default = 0.95).")
  }


  if(length(ntree) != 1) {
    stop("ntree is a positive numeric argument in range (default = 1000).")
  }

  if(ntree <= 0) {
    stop("ntree is a positive numeric argument in range (default = 1000).")
  }

  if(!is.numeric(ntree)) {
    stop("ntree is a positive numeric argument in range (default = 1000).")
  }
}






# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInputSvm <- function(genotype, phenotype, cv.fold,
                          cv.steps, hdi.level) {

  if(is.null(genotype) | is.null(phenotype) | is.null(cv.fold)
     | is.null(cv.steps) | is.null(hdi.level)) {
    stop("arguments must be non-NULL")
  }

  # is it matrix/dataframe or DNAMultipleAlignment/AAMultipleAlignment
  if(is.null(attr(genotype, "class")) == FALSE) {
    if(!attr(genotype, "class") %in% c("AAMultipleAlignment",
                                       "DNAMultipleAlignment")) {
      stop("genotype can be one of the following structures: matrix/data.frame
           or DNAMultipleAlignment/AAMultipleAlignment")
    }
    else {
      temp <- as.matrix(genotype)
      if(nrow(temp) < 2 | ncol(temp) == 0) {
        stop("the genotypes cannot have less than two observations/the number
             of genotypes cannot be 0.")
      }
    }
  }
  else { #matrix/data.frame/vector
    if(is.vector(genotype)) {
      genotype <- matrix(data = genotype, ncol = 1)
    }

    if(nrow(genotype) < 2 | ncol(genotype) == 0) {
      stop("the genotypes cannot have less than two observations/the number
             of genotypes cannot be 0.")
    }

    if(!is.vector(genotype) & !is.matrix(genotype) & !is.data.frame(genotype)) {
      stop("the genotype can be a vector, matrix, data.frame or
           DNAMultipleAlignment/AAMultipleAlignment object")
    }

    if(typeof(genotype) != "character") {
      stop("the genotype entries of character type.")
    }
  }

  if(!is.vector(phenotype)) {
    stop("the phenotype must be a vector.")
  }

  if(length(phenotype) < 2) {
    stop("phenotype cannot be a vector of less than 2 elements.")
  }

  if(!is.numeric(phenotype)) {
    stop("the phenotype must be of numeric type.")
  }

  if(nrow(genotype) != length(phenotype)) {
    stop("the lengths of the genotypes and the phenotypes are not equal. For
         each genotype there should be one phenotype.")
  }

  if(!is.numeric(cv.fold)) {
    stop("cv.fold must be numeric argument in range (0, 1)(default = 0.66).")
  }

  if(cv.fold >= 1 | cv.fold <= 0) {
    stop("cv.fold must be numeric argument in range (0, 1)(default = 0.66).")
  }

  if(!is.numeric(cv.steps)) {
    stop("cv.steps must be a numeric argument larger than 0 (default = 100).")
  }

  if(cv.steps <= 0) {
    stop("cv.steps must be a positive integer (default = 100).")
  }

  if(length(hdi.level) != 1) {
    stop("hdi.level is a numeric argument in range (0, 1)(default = 0.95).")
  }

  if(hdi.level >= 1 | hdi.level <= 0) {
    stop("hdi.level is a numeric argument in range (0, 1)(default = 0.95).")
  }

  if(!is.numeric(hdi.level)) {
    stop("hdi.level is a numeric argument in range (0, 1)(default = 0.95).")
  }
}







# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInputBayes <- function(genotype, phenotype, model,
                            mcmc.iter, chain.nr, hdi.levels) {

  if(is.null(genotype) | is.null(phenotype) | is.null(model) |
     is.null(mcmc.iter) | is.null(chain.nr) | is.null(hdi.levels)) {
    stop("arguments must be non-NULL")
  }

  # CHECK: genotype
  if(is.null(attr(genotype, "class")) == FALSE) {
    if(!attr(genotype, "class") %in% c("AAMultipleAlignment",
                                       "DNAMultipleAlignment")) {
      stop("genotype can be one of the following structures: matrix/data.frame
           or DNAMultipleAlignment/AAMultipleAlignment")
    }
    else {
      temp <- as.matrix(genotype)
      if(nrow(temp) < 2 | ncol(temp) == 0) {
        stop("the genotypes cannot have less than two observations/the number
             of genotypes cannot be 0.")
      }
    }
  }
  else { #matrix/data.frame
    if(is.vector(genotype)) {
      genotype <- matrix(data = genotype, ncol = 1)
    }

    if(nrow(genotype) < 2 | ncol(genotype) == 0) {
      stop("the genotypes cannot have less than two observations/the number
             of genotypes cannot be 0.")
    }

    if(!is.matrix(genotype) & !is.data.frame(genotype)) {
      stop("genotype can be one of the following structures: matrix/data.frame
           or DNAMultipleAlignment/AAMultipleAlignment")
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
    stop("phenotype cannot be a vector of less than 2 elements.")
  }

  if(!is.numeric(phenotype)) {
    stop("the phenotype must be of numeric type.")
  }

  # CHECK: genotype & phenotype
  if(nrow(genotype) != length(phenotype)) {
    stop("the lengths of the genotypes and the phenotypes are not equal. For
         each genotype there should be one phenotype.")
  }

  # CHECK: mcmc.iter
  if(!is.numeric(mcmc.iter)) {
    stop("mcmc.iter must be a numeric argument (default = 10000).")
  }

  if(length(mcmc.iter) != 1) {
    stop("the mcmc.iter must be a number > 0 (default = 10000).")
  }

  if(mcmc.iter == 0) {
    stop("mcmc.iter must be larger than 0 (default = 10000).")
  }

  # CHECK: chain.nr
  if(!is.numeric(chain.nr)) {
    stop("chain.nr must be a numeric argument (default = 2).")
  }

  if(chain.nr == 0) {
    stop("chain.nr must be at least 1 (default = 2).")
  }

  if(length(chain.nr) != 1) {
    stop("chain.nr must be a number > 0 (default = 2).")
  }

  # CHECK: model
  if(!is.character(model)) {
    stop("model must be a character string: 'tdist' or 'pdist'.")
  }

  if(length(model) != 1) {
    stop("model must be a character string: 'tdist' or 'pdist'.")
  }

  if(!model %in% c("tdist", "ndist")) {
    stop("model must be a character string; 'tdist' or 'ndist'.")
  }

  # CHECK: HDI level
  if(length(hdi.levels) < 1) {
    stop("At least one HDI level must br given in the range (0, 1)
         (default = c(0.95, 0.99)).")
  }

  if(any(is.numeric(hdi.levels)) == FALSE) {
    stop("All HDI levels must be numeric in range (0, 1).")
  }

  if(any(hdi.levels >= 1 | hdi.levels <= 0) == TRUE) {
    stop("All HDI levels in range (0, 1).")
  }
}







# Description:
# Set types to each column of the resulting data.frame.
setResultTypesRfSvm <- function(results) {

  # all the col names
  cols <- c("site", "g.1", "g.2", "count.1", "count.2",
            "d", "d.CI.L", "d.CI.H",
            "ca", "ca.hdi.L", "ca.hdi.H", "cv.steps",
            "kappa", "kappa.hdi.L", "kappa.hdi.H",
            "t.test.pvalue")

  # numerical col names
  num.cols <- c( "count.1", "count.2",
                 "d", "d.CI.L", "d.CI.H",
                 "ca", "ca.hdi.L", "ca.hdi.H", "cv.steps",
                 "kappa", "kappa.hdi.L", "kappa.hdi.H",
                 "t.test.pvalue")

  results <- t(t(results))
  results <- data.frame(row.names = NULL, results, stringsAsFactors = FALSE)
  colnames(results) <- cols

  for(c in num.cols) {
    results[[c]] <- as.numeric(results[[c]])
  }

  return (results)
}







# Description:
# Set types to each column of the resulting data.frame.
setResultTypesBayesian <- function(results, hdi.levels) {

  effect.cols <- c()
  for(hdi.level in hdi.levels) {
    effect.cols <- c(effect.cols,
                     paste("mu.effect", hdi.level, "L", sep = '.'),
                     paste("mu.effect", hdi.level, "H", sep = '.'))
  }

  # all the col names
  cols <- c("site", "g.1", "g.2", "count.1", "count.2",
            "t.test.pvalue",
            "mu.1", "sd.1", "mu.2", "sd.2",
            "mu.effect", effect.cols,
            "ess.mu.1", "ess.mu.2")



  # numerical col names
  num.cols <- c("count.1", "count.2",
                "t.test.pvalue",
                "mu.1", "sd.1", "mu.2", "sd.2",
                "mu.effect", effect.cols,
                "ess.mu.1", "ess.mu.2")

  results <- t(t(results))
  results <- data.frame(row.names = NULL, results, stringsAsFactors = FALSE)
  colnames(results) <- cols

  for(c in num.cols) {
    results[[c]] <- as.numeric(results[[c]])
  }

  return (results)
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



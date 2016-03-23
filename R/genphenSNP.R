

# Description:
# Given a matrix of genotypes (SNPs as columns) and a corresponding vector of
# phenotypes, this function computes the association between each SNP and the
# phenotype.
runGenphenSnp <- function(genotype = NULL, phenotype = NULL, technique = NULL,
                          fold.cv = 0.66, boots = 100) {


  # Description:
  # Provided the input arguments, this function checks their validity. It
  # stops the execution if a problem is encountered and prints out warnings.
  checkInputs <- function(genotype, phenotype, technique, fold.cv, boots) {

    if(is.null(genotype) | is.null(phenotype) | is.null(technique)) {
      stop("the arguments genotype, phenotype and technique must be speified.")
    }

    # is it matrix/dataframe or DNAMultipleAlignment
    if(is.null(attr(genotype, "class")) == FALSE) {
      if(!attr(genotype, "class") %in% c("DNAMultipleAlignment")) {
        stop("genotype can be one of the following structures: matrix or
             data.frame or DNAMultipleAlignment")
      }
      else {
        temp <- as.matrix(genotype)
        if(nrow(temp) < 2 | ncol(temp) == 0) {
          stop("the genotype rows cannot be less than 2, and the number of
               columns cannot be 0.")
        }
      }
    }
    else { #matrix/data.frame
      if(nrow(genotype) < 2 | ncol(genotype) == 0) {
        stop("the genotype rows cannot be less than 2, and the number of
               columns cannot be 0.")
      }

      if(!is.matrix(genotype) & !is.data.frame(genotype)) {
        stop("genotype can be one of the following structures: matrix or
             data.frame or DNAMultipleAlignment")
      }

      if(typeof(genotype) != "character") {
        stop("if it is structured in matrix/data.frame the genotype must
             be of character type.")
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

    if(length(technique) != 1) {
      stop("the technique argument must be either rf or svm.")
    }

    if(!technique %in% c("rf", "svm")) {
      stop("the technique argument must be either rf or svm.")
    }

    if(!is.numeric(fold.cv)) {
      stop("fold.cv must be numeric argument in range (0, 1)(default = 0.66).")
    }

    if(fold.cv >= 1 | fold.cv <= 0) {
      stop("fold.cv must be numeric argument in range (0, 1)(default = 0.66).")
    }

    if(!is.numeric(boots)) {
      stop("boots must be larger than 0 (default = 100).")
    }

    if(boots == 0) {
      stop("boots must be larger than 0 (default = 100).")
    }

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
  # Computes the ANOVA score between two vectors (genotype, phenotype)
  # Return: p-value
  computeAnovaScore <- function(entry) {
    anova.summary <- summary(stats::aov(formula = entry$p ~ entry$g))
    if(length(anova.summary[[1]]) != 5) {
      return (NA)
    }
    else {
      return (anova.summary[[1]]$'Pr(>F)'[1])
    }
  }




  # Description:
  # Given two vectors, one dependent (genotype) and one independent
  # (phenotype), compute the classification accuracy of classifying the
  # genotype from the phenotype alone (and corresponding CIs). The
  # classification is non-linear and is computed using random forests.
  # Return: 3-element vector (mean.classification, ci.low, ci.high)
  computeRfClassification <- function(entry, fold.cv, boots) {

    # Description:
    # Perform N number of classifications and compute N number of:
    # - classification accuracy
    # - kappa statistics
    # - number of successful classifications (ideally == N)
    booter <- function(entry, fold.cv, boots) {
      # number of total data entries
      rows <- nrow(entry)

      # initialize result vectors
      accuracies <- c()
      kappas <- c()
      successful.boots <- 0

      for(b in 1:boots) {
        # sample at random
        s <- sample(x = 1:rows, size = round(x = fold.cv * rows, digits = 0),
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
          rf.result <- try(randomForest::randomForest(
            g ~ p, data = train, ntree = 1000), silent = TRUE)
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
            kappa <- computeCohensKappa(predicted = as.character(test$g),
                                        real = as.character(prediction),
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


    class.obj <- booter(entry = entry, fold.cv = fold.cv, boots = boots)


    # build 95% CI for the classification accuracy
    class.obj <- booter(entry = entry, fold.cv = fold.cv, boots = boots)
    ca.ci <- stats::quantile(x = class.obj$accuracies, probs = c(.025, .975),
                             na.rm = TRUE)
    ca.ci.low <- as.numeric(ca.ci[1])
    ca.ci.high <- as.numeric(ca.ci[2])


    # build 95% CI for the kappas
    kappa.ci <- stats::quantile(x = class.obj$kappas, probs = c(.025, .975),
                                na.rm = TRUE)
    kappa.ci.low <- as.numeric(kappa.ci[1])
    kappa.ci.high <- as.numeric(kappa.ci[2])



    # pack results in list
    result <- list(ca = mean(class.obj$accuracies, na.rm = TRUE),
                   ca.ci.low = ca.ci.low,
                   ca.ci.high = ca.ci.high,
                   ca.ci.length = abs(ca.ci.high - ca.ci.low),
                   ca.boots = class.obj$successful.boots,
                   kappa.mean = mean(class.obj$kappas, na.rm = TRUE),
                   kappa.ci.low = kappa.ci.low,
                   kappa.ci.high = kappa.ci.high,
                   kappa.ci.length = abs(kappa.ci.high - kappa.ci.low))


    # if problem return NAs
    if(length(result) != 9) {
      result <- list(ca = NA,
                     ca.ci.low = NA,
                     ca.ci.high = NA,
                     ca.ci.length = NA,
                     ca.boots = NA,
                     kappa.mean = NA,
                     kappa.ci.low = NA,
                     kappa.ci.high = NA,
                     kappa.ci.length = NA)
    }
    return(result)
  }





  # Description:
  # Given two vectors, one dependent (genotype) and one independent
  # (phenotype), compute the classification error of classifying the
  # genotype from the phenotype alone (and corresponding CIs). The
  # classification is linear and is computed using linear support
  # vector machines.
  # Return: 3-element vector (mean.classification, ci.low, ci.high)
  computeSvmClassification <- function(entry, fold.cv, boots) {

    # Description:
    # Perform N number of classifications and compute N number of:
    # - classification accuracies
    # - kappa statistics
    # - number of successful classifications (ideally == N)
    booter <- function(entry, fold.cv, boots) {
      # number of total data entries
      rows <- nrow(entry)

      # initialize result vectors
      accuracies <- c()
      kappas <- c()
      successful.boots <- 0

      for(b in 1:boots) {
        # sample at random
        s <- sample(x = 1:rows, size = round(x = fold.cv * rows, digits = 0),
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
          svm.result <- try(e1071::svm(g ~ p, data = train,
                                       kernel = "linear",
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
            kappa <- computeCohensKappa(predicted = as.character(test$g),
                                        real = as.character(prediction),
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


    class.obj <- booter(entry = entry, fold.cv = fold.cv, boots = boots)


    # build 95% CI for the classification accuracy
    class.obj <- booter(entry = entry, fold.cv = fold.cv, boots = boots)
    ca.ci <- stats::quantile(x = class.obj$accuracies, probs = c(.025, .975),
                             na.rm = TRUE)
    ca.ci.low <- as.numeric(ca.ci[1])
    ca.ci.high <- as.numeric(ca.ci[2])


    # build 95% CI for the kappas
    kappa.ci <- stats::quantile(x = class.obj$kappas, probs = c(.025, .975),
                                na.rm = TRUE)
    kappa.ci.low <- as.numeric(kappa.ci[1])
    kappa.ci.high <- as.numeric(kappa.ci[2])



    # pack results in list
    result <- list(ca = mean(class.obj$accuracies, na.rm = TRUE),
                   ca.ci.low = ca.ci.low,
                   ca.ci.high = ca.ci.high,
                   ca.ci.length = abs(ca.ci.high - ca.ci.low),
                   ca.boots = class.obj$successful.boots,
                   kappa.mean = mean(class.obj$kappas, na.rm = TRUE),
                   kappa.ci.low = kappa.ci.low,
                   kappa.ci.high = kappa.ci.high,
                   kappa.ci.length = abs(kappa.ci.high - kappa.ci.low))


    # if problem return NAs
    if(length(result) != 9) {
      result <- list(ca = NA,
                     ca.ci.low = NA,
                     ca.ci.high = NA,
                     ca.ci.length = NA,
                     ca.boots = NA,
                     kappa.mean = NA,
                     kappa.ci.low = NA,
                     kappa.ci.high = NA,
                     kappa.ci.length = NA)
    }

    return(result)
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
  computeCohensKappa <- function(predicted, real, aas) {

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
    expected.accuracy <- (sum(conf[1, ])/sum(conf) * sum(conf[, 1])/sum(conf))
    + (sum(conf[2, ])/sum(conf) * sum(conf[, 2])/sum(conf))
    real.accuracy <- (conf[1, 1] + conf[2, 2])/sum(conf)
    kappa <- (real.accuracy - expected.accuracy)/(1 - expected.accuracy)

    return (kappa)
  }





  # Description:
  # Given a vector of genotypes, the allele frequencies are computed.
  # Return: 4-element vector (aa1, aa2, aa1.count, aa2.count)
  computeGeneralStats <- function(entry) {
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
  computeCohensD <- function(entry) {
    effect.obj <- effsize::cohen.d.formula(formula = p~g, data = entry)
    effect <- as.numeric(effect.obj$estimate)
    effect.ci.low <- as.numeric(effect.obj$conf.int[1])
    effect.ci.high <- as.numeric(effect.obj$conf.int[2])

    result <- list(effect = effect, effect.ci.low = effect.ci.low,
                   effect.ci.high = effect.ci.high)
    return (result)
  }





  # Description:
  # Set types to each column of the resulting data.frame.
  setResultTypes <- function(results) {
    # all the col names
    cols <- c("site", "allele1", "allele2", "count.allele1", "count.allele2",
              "effect.size", "effect.CI.low", "effect.CI.high",
              "ca", "ca.CI.low", "ca.CI.high", "ca.CI.length", "ca.boots",
              "kappa", "kappa.CI.low", "kappa.CI.high", "kappa.CI.length",
              "anova.score")

    # numerical col names
    num.cols <- c("count.allele1", "count.allele2",
                  "effect.size", "effect.CI.low", "effect.CI.high",
                  "ca", "ca.CI.low", "ca.CI.high", "ca.CI.length", "ca.boots",
                  "kappa", "kappa.CI.low", "kappa.CI.high", "kappa.CI.length",
                  "anova.score")


    results <- t(t(results))
    results <- data.frame(row.names = NULL, results, stringsAsFactors = FALSE)
    colnames(results) <- cols

    for(c in num.cols) {
      results[[c]] <- as.numeric(results[[c]])
    }

    return (results)
  }






  # Description:
  # The main method.
  main <- function(genotype, phenotype, technique, fold.cv, boots) {

    results <- c()
    for(g in 1:ncol(genotype)) {

      # build predictor-response entry
      entry <- data.frame(g = genotype[, g], p = phenotype)

      # the current site position
      site.id <- g

      # in case of no variation (single SNP)
      if(length(unique(entry$g)) <= 1) {
        row <- c(site.id, rep(x = NA, times = 17))
        results <- rbind(results, row)
      }
      else {

        # compute allele counts
        general <- computeGeneralStats(entry = entry)

        if(any(table(entry$g) == 1)) {
          # append empty result with only the site id non-empty
          row <- c(site.id, rep(x = NA, times = 17))
          results <- rbind(results, row)
        }
        else {
          # compute scores
          anova.score <- computeAnovaScore(entry = entry)
          effect <- computeCohensD(entry = entry)

          # classifications
          if(technique == "rf") {
            class <- computeRfClassification(entry = entry, boots = boots,
                                             fold.cv = fold.cv)
          }
          else if(technique == "svm") {
            class <- computeSvmClassification(entry = entry, boots = boots,
                                              fold.cv = fold.cv)
          }

          # append result 1+4+3+1+4+1+4
          row <- c(site.id,
                   general$aa1, general$aa2, general$aa1.nr, general$aa2.nr,
                   effect$effect, effect$effect.ci.low, effect$effect.ci.high,
                   class$ca, class$ca.ci.low, class$ca.ci.high,
                   class$ca.ci.length, class$ca.boots,
                   class$kappa.mean, class$kappa.ci.low, class$kappa.ci.high,
                   class$kappa.ci.length,
                   anova.score)
          results <- rbind(results, row)
        }
      }
    }

    # set correct type of result data
    results <- setResultTypes(results = results)

    return (results)
  }



  # check inputs
  checkInputs(genotype = genotype, phenotype =  phenotype,
              technique = technique, fold.cv = fold.cv,
              boots = boots)



  # convert DNAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)



  # run
  results <- suppressWarnings(expr = main(genotype = genotype,
                                          phenotype = phenotype,
                                          technique = technique,
                                          fold.cv = fold.cv,
                                          boots = boots))

  return (results)
}




runGenphen <- function(genotype = NULL,
                       phenotype = NULL,
                       phenotype.type = NULL,
                       mcmc.chains = 2,
                       mcmc.iterations = 1000,
                       mcmc.warmup = 500,
                       mcmc.cores = 1,
                       hdi.level = 0.95,
                       stat.learn.method = "rf",
                       cv.iterations = 1000) {

  # check inputs
  checkInput(genotype = genotype,
             phenotype = phenotype,
             phenotype.type = phenotype.type,
             mcmc.chains = mcmc.chains,
             mcmc.iterations = mcmc.iterations,
             mcmc.warmup = mcmc.warmup,
             mcmc.cores = mcmc.cores,
             hdi.level = hdi.level,
             stat.learn.method = stat.learn.method,
             cv.iterations = cv.iterations)


  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)


  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }


  # genphen.data
  genphen.data <- getGenphenData(genotype = genotype,
                                 phenotype = phenotype,
                                 phenotype.type = phenotype.type,
                                 min.observations = 3)


  # compile model
  model.stan <- compileModel(phenotype.type = phenotype.type)


  convergence <- NULL
  results <- NULL
  cas <- NULL
  for(s in 1:length(genphen.data)) {
    if(phenotype.type == "continuous") {

      progress.indicator <- round(s/length(genphen.data)*100, digits = 2)
      cat("============================================================= \n")
      cat("======== Main Analysis Progress: ", progress.indicator, "% ",
          "site = ", genphen.data[[s]]$site, "======== \n")
      cat("============================================================= \n")
      o <- runContinuous(data.list = genphen.data[[s]],
                         mcmc.chains = mcmc.chains,
                         mcmc.iterations = mcmc.iterations,
                         mcmc.warmup = mcmc.warmup,
                         mcmc.cores = mcmc.cores,
                         hdi.level = hdi.level,
                         model.stan = model.stan)
    }
    else if(phenotype.type == "dichotomous") {
      progress.indicator <- round(s/length(genphen.data)*100, digits = 2)
      cat("============================================================= \n")
      cat("======== Main Analysis Progress: ", progress.indicator, "% ",
          "site = ", genphen.data[[s]]$site, "======== \n")
      cat("============================================================= \n")
      o <- runDichotomous(data.list = genphen.data[[s]],
                          mcmc.chains = mcmc.chains,
                          mcmc.iterations = mcmc.iterations,
                          mcmc.warmup = mcmc.warmup,
                          mcmc.cores = mcmc.cores,
                          hdi.level = hdi.level,
                          model.stan = model.stan)
    }



    cat("============================================================= \n")
    cat("=================== Statistical Learning ==================== \n")
    cat("============================================================= \n")

    # CA
    if(stat.learn.method == "none") {
      ca <- getNoneCa(data.list = genphen.data[[s]])
    }
    else if(stat.learn.method == "rf") {
      ca <- getRfCa(data.list = genphen.data[[s]],
                    cv.fold = 0.66,
                    cv.steps = cv.iterations,
                    hdi.level = hdi.level,
                    ntree = 1000,
                    mcmc.cores = mcmc.cores)
    }
    else if(stat.learn.method == "svm") {
      ca <- getSvmCa(data.list = genphen.data[[s]],
                     cv.fold = 0.66,
                     cv.steps = cv.iterations,
                     hdi.level = hdi.level,
                     mcmc.cores = mcmc.cores)
    }


    # append results
    results <- rbind(results, o$statistics.out)
    cas <- rbind(cas, ca)
    convergence <- rbind(convergence, o$convergence.out)
  }



  # merge effect sizes and cas
  scores <- merge(x = results, y = cas, all = TRUE,
                  by = c("site", "mutation", "general"))
  if(phenotype.type == "continuous") {
    nice.scores <- scores[, c("site", "mutation", "general",
                              "cohens.d", "cohens.d.L", "cohens.d.H",
                              "ca", "ca.L", "ca.H",
                              "kappa", "kappa.L", "kappa.H",
                              "bc")]
  }
  else if(phenotype.type == "dichotomous") {
    nice.scores <- scores[, c("site", "mutation", "general",
                              "absolute.d", "absolute.d.L", "absolute.d.H",
                              "ca", "ca.L", "ca.H",
                              "kappa", "kappa.L", "kappa.H",
                              "bc")]
  }
  return(list(scores = nice.scores,
              convergence = convergence,
              debug.scores = scores))
}


runDiagnostics <- function(genotype = NULL,
                           phenotype = NULL,
                           phenotype.type = NULL,
                           rf.importance.trees = 50000,
                           with.anchor.points = FALSE,
                           mcmc.chains = 2,
                           mcmc.iterations = 1000,
                           mcmc.warmup = 500,
                           mcmc.cores = 1,
                           hdi.level = 0.95,
                           anchor.points = c(1:5)) {

  # check inputs
  checkInput(genotype = genotype,
             phenotype = phenotype,
             phenotype.type = phenotype.type,
             mcmc.chains = mcmc.chains,
             mcmc.iterations = mcmc.iterations,
             mcmc.warmup = mcmc.warmup,
             mcmc.cores = mcmc.cores,
             hdi.level = hdi.level,
             stat.learn.method = "none",
             cv.iterations = 0)

  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)


  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }


  # check input diagnostics
  checkInputDiagnostics(genotype = genotype,
                        anchor.points = anchor.points,
                        with.anchor.points = with.anchor.points,
                        rf.importance.trees = rf.importance.trees)


  # find importances: prepare for ranger
  rf.data <- data.frame(genotype, stringsAsFactors = FALSE)
  if(phenotype.type == "continuous") {
    rf.data$phenotype <- phenotype
  }
  if(phenotype.type == "dichotomous") {
    rf.data$phenotype <- as.factor(phenotype)
  }


  # ranger: importance dataset
  rf.out <- ranger::ranger(dependent.variable.name = "phenotype",
                           data = rf.data, importance = "impurity",
                           num.trees = rf.importance.trees)
  rf.out <- data.frame(site = 1:length(rf.out$variable.importance),
                       importance = rf.out$variable.importance,
                       stringsAsFactors = FALSE)
  genotype <- genotype[, order(rf.out$importance, decreasing = TRUE)]
  rf.out <- rf.out[order(rf.out$importance, decreasing = TRUE), ]
  rf.out$importance.rank <- 1:nrow(rf.out)


  # if only RF analysis asked, then return importances only
  if(with.anchor.points == FALSE) {
    return (list(scores = NA, importance.scores = rf.out))
  }


  # compile model
  model.stan <- compileModel(phenotype.type = phenotype.type)


  results <- NULL
  cas <- NULL
  for(p in anchor.points) {
    genotype.data <- matrix(data = genotype[, p], ncol = 1)
    phenotype.data <- phenotype


    # genphen.data
    genphen.data <- getGenphenData(genotype = genotype.data,
                                   phenotype = phenotype.data,
                                   phenotype.type = phenotype.type,
                                   min.observations = 3)


    if(length(genphen.data) != 0) {
      for(s in 1:length(genphen.data)) {
        cat("============================================================= \n")
        cat("=========== Anchor point:", p, " progress: ",
            round(s/length(genphen.data)*100, digits = 2), " % =========== \n")
        cat("============================================================= \n")


        if(phenotype.type == "continuous") {
          o <- runContinuous(data.list = genphen.data[[s]],
                             mcmc.chains = mcmc.chains,
                             mcmc.iterations = mcmc.iterations,
                             mcmc.warmup = mcmc.warmup,
                             mcmc.cores = mcmc.cores,
                             hdi.level = hdi.level,
                             model.stan = model.stan)
        }
        else if(phenotype.type == "dichotomous") {
          o <- runDichotomous(data.list = genphen.data[[s]],
                              mcmc.chains = mcmc.chains,
                              mcmc.iterations = mcmc.iterations,
                              mcmc.warmup = mcmc.warmup,
                              mcmc.cores = mcmc.cores,
                              hdi.level = hdi.level,
                              model.stan = model.stan)
        }


        # add marker for diagnostics
        o$statistics.out$anchor.point <- p


        # append results
        results <- rbind(results, o$statistics.out)
      }
    }
    else {
      if(phenotype.type == "continuous") {
        dummy.stats <- data.frame(site = NA, general = NA, mutation = NA,
                                  cohens.d = NA, cohens.d.L = NA,
                                  cohens.d.H = NA, cohens.d.hdi = NA,
                                  bc = 1, anchor.point = p,
                                  stringsAsFactors = FALSE)
      }
      else if(phenotype.type == "dichotomous") {
        dummy.stats <- data.frame(site = NA, general = NA, mutation = NA,
                                  absolute.d = NA, absolute.d.L = NA,
                                  absolute.d.H = NA, absolute.d.hdi = NA,
                                  bc = 1, anchor.point = p,
                                  stringsAsFactors = FALSE)
      }
      # append results
      results <- rbind(results, dummy.stats)
    }
  }

  scores <- results
  scores <- scores[order(scores$anchor.point, decreasing = FALSE), ]
  return(list(scores = scores, importance.scores = rf.out))
}


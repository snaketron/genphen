

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
             cv.iterations = cv.iterations,
             diagnostics.points = 10,
             diagnostics.samples = 10,
             diagnostics.rf.trees = 50000)


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
    if(stat.learn.method == "rf") {
      ca <- getRfCa(data.list = genphen.data[[s]],
                    cv.fold = 0.66,
                    cv.steps = cv.iterations,
                    hdi.level = hdi.level,
                    ntree = 1000)
    }
    else if(stat.learn.method == "svm") {
      ca <- getSvmCa(data.list = genphen.data[[s]],
                     cv.fold = 0.66,
                     cv.steps = cv.iterations,
                     hdi.level = hdi.level)
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
                              "b.coef")]
  }
  else if(phenotype.type == "dichotomous") {
    nice.scores <- scores[, c("site", "mutation", "general",
                              "absolute.d", "absolute.d.L", "absolute.d.H",
                              "ca", "ca.L", "ca.H",
                              "b.coef")]
  }
  return(list(scores = nice.scores,
              convergence = convergence,
              debug.scores = scores))
}




runDiagnostics <- function(genotype = NULL,
                           phenotype = NULL,
                           phenotype.type = NULL,
                           mcmc.chains = 2,
                           mcmc.iterations = 1000,
                           mcmc.warmup = 500,
                           mcmc.cores = 1,
                           hdi.level = 0.95,
                           cv.iterations = 1000,
                           diagnostics.points = 10,
                           diagnostics.samples = 10,
                           diagnostics.rf.trees = 50000) {

  # check inputs
  checkInput(genotype = genotype,
             phenotype = phenotype,
             phenotype.type = phenotype.type,
             mcmc.chains = mcmc.chains,
             mcmc.iterations = mcmc.iterations,
             mcmc.warmup = mcmc.warmup,
             mcmc.cores = mcmc.cores,
             hdi.level = hdi.level,
             stat.learn.method = "rf",
             cv.iterations = cv.iterations,
             diagnostics.points = diagnostics.points,
             diagnostics.samples = diagnostics.samples,
             diagnostics.rf.trees = diagnostics.rf.trees)


  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)


  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }


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
                           num.trees = diagnostics.rf.trees)
  rf.out <- data.frame(site = 1:length(rf.out$variable.importance),
                       importance = rf.out$variable.importance,
                       diagnostics.point = NA,
                       stringsAsFactors = FALSE)
  genotype <- genotype[, order(rf.out$importance, decreasing = TRUE)]


  # compile model
  model.stan <- compileModel(phenotype.type = phenotype.type)


  # now only get the genotypes as used in genphen
  points <- ceiling(seq(from = 1, to = ncol(genotype)-diagnostics.samples,
                        length.out = diagnostics.points))

  # set diagnostcs.points
  rf.out <- rf.out[order(rf.out$importance, decreasing = TRUE), ]
  rf.out$diagnostics.point[points] <- 1:length(points)

  results <- NULL
  cas <- NULL
  for(p in 1:length(points)) {
    genotype.data <- genotype[, points[p]:(points[p]+diagnostics.samples-1)]
    phenotype.data <- phenotype

    # genphen.data
    genphen.data <- getGenphenData(genotype = genotype.data,
                                   phenotype = phenotype.data,
                                   phenotype.type = phenotype.type,
                                   min.observations = 3)

    if(length(genphen.data) != 0) {
      for(s in 1:length(genphen.data)) {
        cat("============================================================= \n")
        cat("========= Diagnostics point:", p, " progress: ",
            round(s/length(genphen.data)*100, digits = 2), " % ========= \n")
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


        # CA
        ca <- getRfCa(data.list = genphen.data[[s]],
                      cv.fold = 0.66,
                      cv.steps = cv.iterations,
                      hdi.level = hdi.level,
                      ntree = 500)


        # add marker for diagnostics
        o$statistics.out$diagnostics.point <- p
        ca$diagnostics.point <- p


        # append results
        results <- rbind(results, o$statistics.out)
        cas <- rbind(cas, ca)
      }
    }
  }



  # merge effect sizes and cas
  scores <- merge(x = results, y = cas, all = TRUE,
                  by = c("site", "mutation", "general", "diagnostics.point"))

  scores <- scores[order(scores$diagnostics.point, decreasing = FALSE), ]
  scores$diagnostics.point <- factor(scores$diagnostics.point,
                                     levels = 1:length(points))

  return(list(scores = scores, rf.out = rf.out))
}








runGenphen <- function(genotype = NULL,
                       phenotype = NULL,
                       phenotype.type = NULL,
                       mcmc.chains = 2,
                       mcmc.iterations = 1000,
                       mcmc.warmup = 500,
                       mcmc.cores = 1,
                       hdi.level = 0.95,
                       stat.learn.method = "rf",
                       cv.iterations = 1000,
                       with.rpa = FALSE,
                       rpa.iterations = 0,
                       rpa.rope = 0) {

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
  rpa <- NULL
  ppc <- NULL
  for(s in 1:length(genphen.data)) {
    if(phenotype.type == "continuous") {

      progress.indicator <- round(s/length(genphen.data)*100, digits = 2)
      cat("======== Main Analysis Progress: ", progress.indicator, "% ",
          "site = ", genphen.data[[s]]$site, "======== \n")
      o <- runContinuous(data.list = genphen.data[[s]],
                         mcmc.chains = mcmc.chains,
                         mcmc.iterations = mcmc.iterations,
                         mcmc.warmup = mcmc.warmup,
                         mcmc.cores = mcmc.cores,
                         hdi.level = hdi.level,
                         model.stan = model.stan,
                         with.rpa = with.rpa,
                         rpa.iterations = rpa.iterations,
                         rpa.rope = rpa.rope)
    }
    else if(phenotype.type == "dichotomous") {
      progress.indicator <- round(s/length(genphen.data)*100, digits = 2)
      cat("======== Main Analysis Progress: ", progress.indicator, "% ",
          "site = ", genphen.data[[s]]$site, "======== \n")
      o <- runDichotomous(data.list = genphen.data[[s]],
                          mcmc.chains = mcmc.chains,
                          mcmc.iterations = mcmc.iterations,
                          mcmc.warmup = mcmc.warmup,
                          mcmc.cores = mcmc.cores,
                          hdi.level = hdi.level,
                          model.stan = model.stan,
                          with.rpa = with.rpa,
                          rpa.iterations = rpa.iterations,
                          rpa.rope = rpa.rope)
    }



    cat("=================== Statistical Learning ==================== \n")
    # CA
    if(stat.learn.method == "none") {
      ca <- getNoneCa(data.list = genphen.data[[s]])
    }
    else if(stat.learn.method == "rf") {
      ca <- getRfCa(data.list = genphen.data[[s]],
                    cv.fold = 0.7,
                    cv.steps = cv.iterations,
                    hdi.level = hdi.level,
                    ntree = 1000,
                    mcmc.cores = mcmc.cores)
    }
    else if(stat.learn.method == "svm") {
      ca <- getSvmCa(data.list = genphen.data[[s]],
                     cv.fold = 0.7,
                     cv.steps = cv.iterations,
                     hdi.level = hdi.level,
                     mcmc.cores = mcmc.cores)
    }


    # append results
    results <- rbind(results, o$statistics.out)
    cas <- rbind(cas, ca)
    convergence <- rbind(convergence, o$convergence.out)
    rpa <- rbind(rpa, o$rpa.out)
    ppc <- rbind(ppc, o$ppc.out)
  }



  # merge effect sizes and cas
  scores <- merge(x = results, y = cas, all = TRUE,
                  by = c("site", "mutation", "general"))
  if(phenotype.type == "continuous") {
    nice.scores <- scores[, c("site", "mutation", "general",
                              "cohens.d", "cohens.d.L", "cohens.d.H",
                              "ca", "ca.L", "ca.H",
                              "kappa", "kappa.L", "kappa.H",
                              "bc",
                              "sd.d", "sd.d.L", "sd.d.H")]
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
              debug.scores = scores,
              rpa = rpa,
              ppc = ppc))
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
                             model.stan = model.stan,
                             with.rpa = FALSE,
                             rpa.iterations = 0,
                             rpa.rope = 0)
        }
        else if(phenotype.type == "dichotomous") {
          o <- runDichotomous(data.list = genphen.data[[s]],
                              mcmc.chains = mcmc.chains,
                              mcmc.iterations = mcmc.iterations,
                              mcmc.warmup = mcmc.warmup,
                              mcmc.cores = mcmc.cores,
                              hdi.level = hdi.level,
                              model.stan = model.stan,
                              with.rpa = FALSE,
                              rpa.iterations = 0,
                              rpa.rope = 0)
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




runPhyloBiasCheck <- function(input.kinship.matrix = NULL,
                              genotype = NULL) {

  # check params
  checkInputPhyloBias(input.kinship.matrix = input.kinship.matrix,
                      genotype = genotype)

  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)

  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }

  # compute kinship if needed
  if(is.null(input.kinship.matrix) | missing(input.kinship.matrix)) {
    kinship.matrix <- e1071::hamming.distance(genotype)
  }
  else {
    kinship.matrix <- input.kinship.matrix
  }

  # compute bias
  bias <- getPhyloBias(genotype = genotype, k.matrix = kinship.matrix)

  # bias = 1-dist(feature)/dist(total)
  bias$bias <- 1-bias$feature.dist/bias$total.dist

  # get mutations
  bias.mutations <- c()
  sites <- unique(bias$site)
  for(s in sites) {
    s.bias <- bias[bias$site == s, ]
    genotypes <- sort(unique(s.bias$genotype))
    if(length(genotypes) != 1) {
      for(g1 in 1:(length(genotypes) - 1)) {
        bias1 <- s.bias$bias[s.bias$genotype == genotypes[g1]]
        for(g2 in (g1 + 1):length(genotypes)) {
          bias2 <- s.bias$bias[s.bias$genotype == genotypes[g2]]
          mutation <- paste(genotypes[g1], genotypes[g2], sep = "->")
          mutation.row <- data.frame(site = s, mutation = mutation,
                                     bias = max(bias1, bias2))
          bias.mutations <- rbind(bias.mutations, mutation.row)
        }
      }
    }
  }


  return (list(bias = bias[, c("site", "genotype", "bias")],
               kinship.matrix = kinship.matrix,
               bias.mutations = bias.mutations))
}



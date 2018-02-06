

runGenphen <- function(genotype = NULL,
                       phenotype = NULL,
                       phenotype.type = NULL,
                       mcmc.chains = 2,
                       mcmc.iterations = 1000,
                       mcmc.warmup = 500,
                       mcmc.cores = 1,
                       hdi.level = 0.95,
                       with.ppc = FALSE,
                       stat.learn.method = "rf") {

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
             with.ppc = with.ppc)


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


  convergence <- NULL
  results <- NULL
  ppc <- NULL
  cas <- NULL
  for(s in 1:length(genphen.data)) {
    if(phenotype.type == "continuous") {
      # model.stan <- get(load(file = "inst/extdata/continuous.stan.RData"))
      model.stan <- get(load(system.file("extdata", "continuous.stan.RData",
                                         package = "genphen")))
      cat("============================================================= \n")
      cat("Progress: ", round(s/length(genphen.data)*100, digits = 2), "% ",
          "site = ", genphen.data[[s]]$site, "\n")
      cat("============================================================= \n")


      o <- runContinuous(data.list = genphen.data[[s]],
                         mcmc.chains = mcmc.chains,
                         mcmc.iterations = mcmc.iterations,
                         mcmc.warmup = mcmc.warmup,
                         mcmc.cores = mcmc.cores,
                         hdi.level = hdi.level,
                         model.stan = model.stan,
                         with.ppc = with.ppc)
      rm(model.stan)
      gc(verbose = FALSE)
    }
    else if(phenotype.type == "dichotomous") {
      # model.stan <- get(load(file = "inst/extdata/dichotomous.stan.RData"))
      model.stan <- get(load(system.file("extdata", "dichotomous.stan.RData",
                                         package = "genphen")))
      cat("============================================================= \n")
      cat("Progress: ", round(s/length(genphen.data)*100, digits = 2), "% ",
          "site = ", genphen.data[[s]]$site, "\n")
      cat("============================================================= \n")

      o <- runDichotomous(data.list = genphen.data[[s]],
                          mcmc.chains = mcmc.chains,
                          mcmc.iterations = mcmc.iterations,
                          mcmc.warmup = mcmc.warmup,
                          mcmc.cores = mcmc.cores,
                          hdi.level = hdi.level,
                          model.stan = model.stan,
                          with.ppc = with.ppc)
      rm(model.stan)
      gc(verbose = FALSE)
    }


    # CA and Kappa
    if(stat.learn.method == "rf") {
      ca <- getRfCa(data.list = genphen.data[[s]],
                    cv.fold = 0.66,
                    cv.steps = 1000,
                    hdi.level = hdi.level,
                    ntree = 1000)
    }
    else if(stat.learn.method == "svm") {
      ca <- getSvmCa(data.list = genphen.data[[s]],
                     cv.fold = 0.66,
                     cv.steps = 1000,
                     hdi.level = hdi.level)
    }


    # append results
    results <- rbind(results, o$statistics.out)
    cas <- rbind(cas, ca)
    convergence <- rbind(convergence, o$convergence.out)
    ppc <- rbind(ppc, o$ppc)
  }


  # merge effect sizes and cas
  scores <- merge(x = results, y = cas, all = TRUE,
                  by = c("site", "mutation", "general"))
  if(phenotype.type == "continuous") {
    nice.scores <- scores[, c("site", "mutation", "general",
                              "cohens.d", "cohens.d.L", "cohens.d.H",
                              "ca", "ca.L", "ca.H",
                              "kappa", "kappa.L", "kappa.H")]
  }
  else if(phenotype.type == "dichotomous") {
    nice.scores <- scores[, c("site", "mutation", "general",
                              "absolute.d", "absolute.d.L", "absolute.d.H",
                              "ca", "ca.L", "ca.H",
                              "kappa", "kappa.L", "kappa.H")]
  }
  return(list(scores = nice.scores,
              ppc = ppc,
              convergence = convergence,
              debug.scores = scores))
}









runGenphenRf <- function(genotype = NULL, phenotype = NULL,
                        cv.fold = 0.66, cv.steps = 100,
                        hdi.level = 0.99, ntree = 1000) {



  # Description:
  # The main procedure which executes all the necessary procedures.
  main <- function(genotype, phenotype, cv.fold, cv.steps, hdi.level, ntree) {

    #store all results
    results <- c()
    # loop through each genotype (MSA) site
    for(g in 1:ncol(genotype)) {

      # get the current amino acid site
      site <- as.character(genotype[, g])

      # the current site position
      site.id <- g

      # get the unique amino acid states at a site
      states <- unique(site)

      if(length(states) <= 1) {
        # append empty result with only the site id non-empty
        row <- c(site.id, rep(x = NA, times = 15))
        results <- rbind(results, row)
      }
      else {
        # substitution pair (s1, s2) analysis
        for(s1 in 1:(length(states) - 1)) {
          for(s2 in (s1 + 1):length(states)) {

            # build entry
            index <- which(site %in% c(states[s1], states[s2]))
            entry <- data.frame(g = site[index], p = phenotype[index])


            # compute AA counts
            general <- getGeneralStats(entry = entry)


            if(any(table(entry$g) == 1)) {
              # append result
              row <- c(site.id, general$aa1, general$aa2, general$aa1.nr,
                       general$aa2.nr, rep(x = NA, times = 11))
              results <- rbind(results, row)
            }
            else {
              # compute scores
              ttest.score <- getTTestScore(entry = entry)
              effect <- getCohensD(entry = entry, conf.level = hdi.level)

              class <- getRfClassification(entry = entry,
                                           cv.steps = cv.steps,
                                           cv.fold = cv.fold,
                                           hdi.level = hdi.level,
                                           ntree = ntree)

              # append result 16 entries
              row <- c(site.id,
                       general$aa1, general$aa2, general$aa1.nr, general$aa2.nr,
                       effect$effect, effect$effect.ci.low,
                       effect$effect.ci.high,
                       class$ca, class$ca.hdi.low, class$ca.hdi.high,
                       class$ca.boots,
                       class$kappa, class$kappa.hdi.low, class$kappa.hdi.high,
                       ttest.score)
              results <- rbind(results, row)
            }
          }
        }
      }
    }

    # set correct types to each results column
    results <- setResultTypesRfSvm(results = results)

    return (results)
  }


  # check inputs
  checkInputRf(genotype = genotype, phenotype = phenotype,
               cv.fold = cv.fold, cv.steps = cv.steps,
               hdi.level = hdi.level, ntree = ntree)


  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)


  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }

  # run
  results <- suppressWarnings(expr = main(genotype = genotype,
                                          phenotype = phenotype,
                                          cv.fold = cv.fold,
                                          cv.steps = cv.steps,
                                          hdi.level = hdi.level,
                                          ntree = ntree))

  return (results)
}







runGenphenSvm <- function(genotype = NULL, phenotype = NULL,
                         cv.fold = 0.66, cv.steps = 100,
                         hdi.level = 0.99) {



  # Description:
  # The main procedure which executes all the necessary procedures.
  main <- function(genotype, phenotype, cv.fold, cv.steps, hdi.level) {

    #store all results
    results <- c()
    # loop through each genotype (MSA) site
    for(g in 1:ncol(genotype)) {

      # get the current amino acid site
      site <- as.character(genotype[, g])

      # the current site position
      site.id <- g

      # get the unique amino acid states at a site
      states <- unique(site)

      if(length(states) <= 1) {
        # append empty result with only the site id non-empty
        row <- c(site.id, rep(x = NA, times = 15))
        results <- rbind(results, row)
      }
      else {
        # substitution pair (s1, s2) analysis
        for(s1 in 1:(length(states) - 1)) {
          for(s2 in (s1 + 1):length(states)) {

            # build entry
            index <- which(site %in% c(states[s1], states[s2]))
            entry <- data.frame(g = site[index], p = phenotype[index])


            # compute AA counts
            general <- getGeneralStats(entry = entry)


            if(any(table(entry$g) == 1)) {
              # append result
              row <- c(site.id, general$aa1, general$aa2, general$aa1.nr,
                       general$aa2.nr, rep(x = NA, times = 11))
              results <- rbind(results, row)
            }
            else {
              # compute scores
              ttest.score <- getTTestScore(entry = entry)
              effect <- getCohensD(entry = entry, conf.level = hdi.level)

              class <- getSvmClassification(entry = entry,
                                            cv.steps = cv.steps,
                                            cv.fold = cv.fold,
                                            hdi.level = hdi.level)

              # append result
              row <- c(site.id,
                       general$aa1, general$aa2, general$aa1.nr,
                       general$aa2.nr,
                       effect$effect, effect$effect.ci.low,
                       effect$effect.ci.high,
                       class$ca, class$ca.hdi.low, class$ca.hdi.high,
                       class$ca.boots,
                       class$kappa, class$kappa.hdi.low, class$kappa.hdi.high,
                       ttest.score)
              results <- rbind(results, row)
            }
          }
        }
      }
    }

    # set correct types to each results column
    results <- setResultTypesRfSvm(results = results)

    return (results)
  }


  # check inputs
  checkInputSvm(genotype = genotype, phenotype =  phenotype,
                cv.fold = cv.fold, cv.steps = cv.steps,
                hdi.level = hdi.level)


  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)

  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }

  # run
  results <- suppressWarnings(expr = main(genotype = genotype,
                                          phenotype = phenotype,
                                          cv.fold = cv.fold,
                                          cv.steps = cv.steps,
                                          hdi.level = hdi.level))

  return (results)
}







runGenphenBayes <- function(genotype = NULL, phenotype = NULL,
                            chain.nr = 2, mcmc.iter = 20000,
                            model = "ndist",
                            hdi.levels = c(0.95, 0.99)) {


  # Description:
  # The main procedure which executes all the necessary procedures.
  main <- function(genotype, phenotype,
                   mcmc.iter, chain.nr,
                   model, hdi.levels) {

    #store all results
    results <- c()
    # loop through each genotype (MSA) site
    for(g in 1:ncol(genotype)) {

      # get the current amino acid site
      site <- as.character(genotype[, g])

      # the current site position
      site.id <- g

      # get the unique amino acid states at a site
      states <- unique(site)

      if(length(states) <= 1) {
        # append empty result with only the site id non-empty
        row <- c(site.id, rep(x = NA, times = 12),
                 rep(x = NA, times = length(hdi.levels) * 2))
        results <- rbind(results, row)
      }
      else {
        # substitution pair (s1, s2) analysis
        for(s1 in 1:(length(states) - 1)) {
          for(s2 in (s1 + 1):length(states)) {

            # build entry
            index <- which(site %in% c(states[s1], states[s2]))
            entry <- data.frame(g = site[index], p = phenotype[index])

            # compute AA counts
            general <- getGeneralStats(entry = entry)

            if(any(table(entry$g) == 1)) {
              # append result 1+4+3+1+4+1+4
              row <- c(site.id, general$aa1, general$aa2, general$aa1.nr,
                       general$aa2.nr, rep(x = NA, times = 8),
                       rep(x = NA, times = length(hdi.levels) * 2))
              results <- rbind(results, row)
            }
            else {
              # compute scores
              ttest.score <- getTTestScore(entry = entry)


              btest <- getBayesianTtest(entry = entry,
                                        n.iter = mcmc.iter,
                                        n.chains = chain.nr,
                                        hdi.levels = hdi.levels,
                                        model = model)

              bres <- c(btest$result[[1]]$mu.i, btest$result[[1]]$sd.i,
                        btest$result[[1]]$mu.j, btest$result[[1]]$sd.j,
                        btest$result[[1]]$mu.diff)
              for(hdi.level in hdi.levels) {
                temp <- btest$result[[as.character(hdi.level)]]
                bres <- c(bres, temp$hdi.diff.L, temp$hdi.diff.H)
              }
              bess <- btest$ess

              # append result 1+4+3+1+4+1+4+11+2
              row <- c(site.id,
                       general$aa1, general$aa2, general$aa1.nr,
                       general$aa2.nr,
                       ttest.score,
                       bres,
                       bess$ess.mu.i, bess$ess.mu.j)
              results <- rbind(results, row)
            }
          }
        }
      }
    }

    # set correct types to each results column
    results <- setResultTypesBayesian(results = results,
                                      hdi.levels = hdi.levels)

    return (results)
  }


  # check inputs
  checkInputBayes(genotype = genotype, phenotype =  phenotype,
                  mcmc.iter = mcmc.iter, chain.nr = chain.nr,
                  model = model, hdi.levels = hdi.levels)


  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)

  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }

  # run
  results <- suppressWarnings(expr = main(genotype = genotype,
                                          phenotype = phenotype,
                                          mcmc.iter = mcmc.iter,
                                          chain.nr = chain.nr,
                                          model = model,
                                          hdi.levels = hdi.levels))

  return(results)
}











runGenphen <- function(genotype = NULL,
                       phenotype = NULL,
                       phenotype.type = NULL,
                       model.type = NULL,
                       mcmc.chains = 2,
                       mcmc.iterations = 2500,
                       mcmc.warmup = 500,
                       cores = 1,
                       hdi.level = 0.95,
                       stat.learn.method = "rf",
                       cv.iterations = 1000,
                       rpa.iterations = 0,
                       with.stan.obj = FALSE) {
  
  
  # check inputs
  checkInput(genotype = genotype,
             phenotype = phenotype,
             phenotype.type = phenotype.type,
             model.type = model.type,
             mcmc.chains = mcmc.chains,
             mcmc.iterations = mcmc.iterations,
             mcmc.warmup = mcmc.warmup,
             cores = cores,
             hdi.level = hdi.level,
             stat.learn.method = stat.learn.method,
             cv.iterations = cv.iterations,
             rpa.iterations = rpa.iterations,
             with.stan.obj = with.stan.obj)
  
  
  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)
  
  
  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }
  
  
  # if phenotype = dichotomous => convert to 1s and 0s
  if(phenotype.type == "dichotomous") {
    phenotype.new <- as.numeric(as.factor(phenotype))-1
    cat("Phenotype mapping to 1s and 0s: \n")
    print(table(phenotype.new, phenotype))
  }
  
  
  # genphen.data and final check for input
  genphen.data <- getGenphenData(genotype = genotype,
                                 phenotype = phenotype,
                                 cores = cores)
  if(is.null(genphen.data)) {
    stop("No genphen input data found.")
  }
  
  
  # compile model
  cat("======== Compiling Main Model ======== \n")
  model.stan <- compileModel(phenotype.type = phenotype.type, 
                             model.type = model.type)
  
  
  j <- NULL
  if(model.type == "univariate") {
    cat("======== Main Analysis Running (Univariate) ======== \n")
    # register parallel
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    if(phenotype.type == "continuous") {
      o <- (foreach(j = 1:max(genphen.data$J),
                    .export = c("runContU", "getHdi", "getGenphenData",
                                "getBhattacharyya", "getRpaCont"), 
                    .packages = c("rstan"))
            %dopar% runContU(genphen.data = genphen.data[genphen.data$J == j, ],
                             mcmc.chains = mcmc.chains,
                             mcmc.iterations = mcmc.iterations,
                             mcmc.warmup = mcmc.warmup,
                             cores = 1,
                             hdi.level = hdi.level,
                             model.stan = model.stan,
                             rpa.iterations = rpa.iterations,
                             with.stan.obj = with.stan.obj))
    }
    else if(phenotype.type == "dichotomous") {
      o <- (foreach(j = 1:max(genphen.data$J),
                    .export = c("runDichU", "getHdi", "getGenphenData",
                                "getBhattacharyya", "getRpaDich"), 
                    .packages = c("rstan"))
            %dopar% runDichU(genphen.data = genphen.data[genphen.data$J == j, ],
                             mcmc.chains = mcmc.chains,
                             mcmc.iterations = mcmc.iterations,
                             mcmc.warmup = mcmc.warmup,
                             cores = 1,
                             hdi.level = hdi.level,
                             model.stan = model.stan,
                             rpa.iterations = rpa.iterations,
                             with.stan.obj = with.stan.obj))
    }
    # stop cluster
    parallel::stopCluster(cl = cl)
    doParallel::stopImplicitCluster()
  }
  
  if(model.type == "hierarchical") {
    cat("======== Main Analysis Running (Hierarchical) ======== \n")
    if(phenotype.type == "continuous") {
      o <- runContH(genphen.data = genphen.data,
                    mcmc.chains = mcmc.chains,
                    mcmc.iterations = mcmc.iterations,
                    mcmc.warmup = mcmc.warmup,
                    cores = cores,
                    hdi.level = hdi.level,
                    model.stan = model.stan,
                    rpa.iterations = rpa.iterations,
                    with.stan.obj = with.stan.obj)
    }
    else if(phenotype.type == "dichotomous") {
      o <- runDichH(genphen.data = genphen.data,
                    mcmc.chains = mcmc.chains,
                    mcmc.iterations = mcmc.iterations,
                    mcmc.warmup = mcmc.warmup,
                    cores = cores,
                    hdi.level = hdi.level,
                    model.stan = model.stan,
                    rpa.iterations = rpa.iterations,
                    with.stan.obj = with.stan.obj)
    }
  }
  
  
  cat("======== Statistical Learning ======== \n")
  # register parallel
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  if(stat.learn.method == "none") {
    cas <- (foreach(j = 1:max(genphen.data$J),
                    .export = c("getNoneCa"),
                    .packages = c("ranger")) %dopar%
              getNoneCa(genphen.data = genphen.data[genphen.data$J == j, ]))
  }
  else if(stat.learn.method == "rf") {
    cas <- (foreach(j = 1:max(genphen.data$J),
                    .export = c("getRfCa", "getHdi", "getKappa"),
                    .packages = c("ranger")) %dopar%
              getRfCa(genphen.data = genphen.data[genphen.data$J == j, ],
                      cv.fold = 0.66,
                      cv.steps = cv.iterations,
                      hdi.level = hdi.level,
                      ntree = 1000))
  }
  else if(stat.learn.method == "svm") {
    cas <- (foreach(j = 1:max(genphen.data$J),
                    .export = c("getSvmCa", "getHdi", "getKappa"), 
                    .packages = c("e1071")) %dopar% 
              getSvmCa(genphen.data = genphen.data[genphen.data$J == j, ],
                       cv.fold = 0.66,
                       cv.steps = cv.iterations,
                       hdi.level = hdi.level))
  }
  # stop cluster
  parallel::stopCluster(cl = cl)
  doParallel::stopImplicitCluster()
  
  
  # assemble results - univariate
  if(model.type == "univariate") {
    getO <- function(x, y) {
      return(x[[y]])
    }
    getS <- function(x, y) {
      return(x[[y]])
    }
    
    
    results <- do.call(rbind, lapply(X = o, FUN = getO, y="statistics.out"))
    convergence <- do.call(rbind, lapply(X = o, FUN = getO,y="convergence.out"))
    if(rpa.iterations == 0) {
      rpa <- NULL
    }
    else {
      rpa <- do.call(rbind, lapply(X = o, FUN = getO, y="rpa.out"))
    }
    ppc <- do.call(rbind, lapply(X = o, FUN = getO, y="ppc.out"))
    cas <- do.call(rbind, cas)
    stan.obj <- lapply(X = o, FUN = getS, y = "stan.obj")
  }
  
  
  # assemble results - hierarchical
  if(model.type == "hierarchical") {
    results <- o[["statistics.out"]]
    convergence <- o[["convergence.out"]]
    rpa <- o[["rpa.out"]]
    ppc <- o[["ppc.out"]]
    cas <- do.call(rbind, cas)
    stan.obj <- o[["stan.obj"]]
  }
  
  
  # merge effect sizes and cas
  scores <- merge(x = results, y = cas, all = TRUE,
                  by = c("site", "g1", "g0", "n1", "n0"))
  
  
  if(phenotype.type == "continuous") {
    # format scores nicely
    nice.scores <- scores[, c("site", "g1", "g0",
                              "beta.mean", "beta.L", "beta.H", "beta.sd",
                              "alpha.mean", "alpha.L", "alpha.H", "alpha.sd",
                              "sigma.mean", "sigma.L", "sigma.H", "sigma.sd",
                              "nu.mean", "nu.L", "nu.H", "nu.sd",
                              "ca", "ca.L", "ca.H",
                              "kappa", "kappa.L", "kappa.H",
                              "bc")]
  }
  if(phenotype.type == "dichotomous") {
    # format scores nicely
    nice.scores <- scores[, c("site", "g1", "g0",
                              "beta.mean", "beta.L", "beta.H", "beta.sd",
                              "alpha.mean", "alpha.L", "alpha.H", "alpha.sd",
                              "ca", "ca.L", "ca.H",
                              "kappa", "kappa.L", "kappa.H",
                              "bc")]
  }
  
  
  # order by genotype site
  nice.scores <- nice.scores[order(nice.scores$site, decreasing = FALSE), ]
  convergence <- convergence[order(convergence$site, decreasing = FALSE), ]
  scores <- scores[order(scores$site, decreasing = FALSE), ]
  if(is.null(rpa) == FALSE) {
    rpa <- rpa[order(rpa$site, decreasing = FALSE), ]
  }
  
  # return
  return(list(scores = nice.scores,
              convergence = convergence,
              debug.scores = scores,
              rpa = rpa,
              ppc = ppc,
              stan.obj = stan.obj))
}




runDiagnostics <- function(genotype = NULL,
                           phenotype = NULL,
                           phenotype.type = NULL,
                           rf.trees = 5000,
                           mcmc.chains = 2,
                           mcmc.iterations = 2500,
                           mcmc.warmup = 500,
                           cores = 1,
                           hdi.level = 0.95,
                           diagnostic.points = NULL) {
  
  
  # check inputs
  checkInput(genotype = genotype,
             phenotype = phenotype,
             phenotype.type = phenotype.type,
             model.type = "univariate",
             mcmc.chains = mcmc.chains,
             mcmc.iterations = mcmc.iterations,
             mcmc.warmup = mcmc.warmup,
             cores = cores,
             hdi.level = hdi.level,
             stat.learn.method = "none",
             cv.iterations = 0,
             rpa.iterations = 0,
             with.stan.obj = FALSE)
  
  
  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)
  
  
  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }
  
  
  # check input diagnostics
  checkInputDiagnostics(genotype = genotype,
                        diagnostic.points = diagnostic.points,
                        rf.trees = rf.trees)
  
  
  # find importances: prepare for ranger
  rf.data <- data.frame(genotype)
  if(phenotype.type == "continuous") {
    rf.data$phenotype <- phenotype
  }
  if(phenotype.type == "dichotomous") {
    rf.data$phenotype <- as.factor(phenotype)
  }
  
  
  # ranger: importance dataset
  cat("======== RF diagnostics ======== \n")
  rf.out <- ranger::ranger(dependent.variable.name = "phenotype",
                           importance = "impurity",
                           data = rf.data, 
                           num.trees = rf.trees)
  rf.out <- data.frame(site = 1:length(rf.out$variable.importance),
                       importance = rf.out$variable.importance,
                       stringsAsFactors = FALSE)
  rf.out <- rf.out[order(rf.out$importance, decreasing = TRUE), ]
  rf.out$importance.rank <- 1:nrow(rf.out)
  
  
  
  # if diagnostic points are not provided, return the RF scores
  if(is.null(diagnostic.points) == TRUE ||
     length(diagnostic.points) == 0 ||
     is.na(diagnostic.points) == TRUE) {
    cat("No diagnostic points provided, only importance analysis performed.\n")
    return(list(scores = NA, 
                importance.scores = rf.out))
  }
  
  
  
  # genphen.data and final check for input
  genphen.data <- getGenphenData(genotype = genotype,
                                 phenotype = phenotype,
                                 cores = cores)
  
  
  if(is.null(genphen.data)) {
    stop("No genphen input data found with these diagnostic points.")
  }
  
  # in case of duplicate diagnostic points
  diagnostic.points <- unique(diagnostic.points)
  # actual anchor points
  anchors <- data.frame(real.anchors = rf.out$site[diagnostic.points],
                        diagnostic.points = diagnostic.points,
                        stringsAsFactors = FALSE)
  # check if given anchor points are too conserved to be used for the analysis
  anchors$miss <- ifelse(test = !anchors$real.anchors %in% genphen.data$S,
                         yes = TRUE, no = FALSE)
  if(sum(anchors$miss) != 0) {
    warning(paste(sum(anchors$miss), " anchor points are conserved, cannot be ",
                  "analyzed: ", paste(anchors$real.anchors[anchors$miss==TRUE], 
                                      collapse = ','), sep =''))
  }
  
  
  # final genphen data
  real.anchors <- anchors$real.anchors[anchors$miss == FALSE]
  genphen.data <- genphen.data[genphen.data$S %in% real.anchors, ]
  if(is.null(genphen.data)) {
    stop("No genphen input data found with these anchor points.")
  }
  genphen.data <- merge(x = genphen.data, y = anchors,
                        by.x = "S", by.y = "real.anchors")
  
  
  # compile model
  cat("======== Compiling Diagnostic Model ======== \n")
  model.stan <- compileModel(phenotype.type = phenotype.type, 
                             model.type = "univariate")
  
  
  cat("======== Main Analysis Running ======== \n")
  # register parallel
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  Js <- unique(genphen.data$J)
  j <- NULL
  if(phenotype.type == "continuous") {
    o <- (foreach(j = 1:length(Js),
                  .export = c("runContU", "getHdi", "getGenphenData",
                              "getBhattacharyya", "getRpaCont"), 
                  .packages = c("rstan"))
          %dopar% runContU(genphen.data = genphen.data[genphen.data$J==Js[j],],
                           mcmc.chains = mcmc.chains,
                           mcmc.iterations = mcmc.iterations,
                           mcmc.warmup = mcmc.warmup,
                           cores = 1,
                           hdi.level = hdi.level,
                           model.stan = model.stan,
                           rpa.iterations = 0))
  }
  else if(phenotype.type == "dichotomous") {
    o <- (foreach(j = 1:length(Js),
                  .export = c("runDichU", "getHdi", "getGenphenData",
                              "getBhattacharyya", "getRpaDich"), 
                  .packages = c("rstan"))
          %dopar% runDichU(genphen.data = genphen.data[genphen.data$J==Js[j],],
                           mcmc.chains = mcmc.chains,
                           mcmc.iterations = mcmc.iterations,
                           mcmc.warmup = mcmc.warmup,
                           cores = 1,
                           hdi.level = hdi.level,
                           model.stan = model.stan,
                           rpa.iterations = 0))
  }
  # stop cluster
  parallel::stopCluster(cl = cl)
  doParallel::stopImplicitCluster()
  
  
  # assemble results - univariate
  getO <- function(x, y) {
    return(x[[y]])
  }
  scores <- do.call(rbind, lapply(X = o, FUN = getO, y="statistics.out"))
  
  
  # merge
  scores <- merge(x = scores, 
                  y = anchors[, c("real.anchors", "diagnostic.points")], 
                  by.x = "site", by.y = "real.anchors", all.x = TRUE)
  
  
  
  if(phenotype.type == "continuous") {
    # format scores nicely
    nice.scores <- scores[, c("site", "g1", "g0",
                              "beta.mean", "beta.L", "beta.H", "beta.sd",
                              "alpha.mean", "alpha.L", "alpha.H", "alpha.sd",
                              "sigma.mean", "sigma.L", "sigma.H", "sigma.sd",
                              "nu.mean", "nu.L", "nu.H", "nu.sd",
                              "bc", "diagnostic.points")]
  }
  if(phenotype.type == "dichotomous") {
    # format scores nicely
    nice.scores <- scores[, c("site", "g1", "g0",
                              "beta.mean", "beta.L", "beta.H", "beta.sd",
                              "alpha.mean", "alpha.L", "alpha.H", "alpha.sd",
                              "bc", "diagnostic.points")]
  }
  
  
  return(list(scores = nice.scores, 
              importance.scores = rf.out))
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
  
  # gen.data and final check for input
  gen.data <- getGenSummary(genotype = genotype)
  if(is.null(gen.data)) {
    stop("No genphen input data found.")
  }
  
  # append bias to each SNP
  gen.data$bias.g1 <- NA
  gen.data$bias.g0 <- NA
  gen.data$bias <- NA
  for(i in 1:nrow(gen.data)) {
    bias.g1 <- bias[bias$site == gen.data$site[i] 
                    & bias$genotype == gen.data$g1[i], ]
    bias.g0 <- bias[bias$site == gen.data$site[i] 
                    & bias$genotype == gen.data$g0[i], ]
    gen.data$bias.g1[i] <- bias.g1$bias[1]
    gen.data$bias.g0[i] <- bias.g0$bias[1]
    gen.data$bias[i] <- max(bias.g1$bias[1], bias.g0$bias[1])
  }
  
  
  # sort by site
  bias <- gen.data[, c("site", "g1", "g0", "bias.g1", "bias.g0", "bias")]
  bias <- bias[order(bias$site, decreasing = FALSE), ]
  
  return (list(bias = bias, kinship.matrix = kinship.matrix))
}



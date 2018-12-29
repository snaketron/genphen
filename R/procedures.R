


# Description:
# Get genotype-phenotype data in a multicore fashion
getGenphenData <- function(genotype, 
                           phenotype, 
                           cores) {
  
  
  # set site id
  colnames(genotype) <- 1:ncol(genotype)
  
  
  # Description:
  # Get genotype-phenotype data
  getData <- function(genotype, phenotype) {
    # min.obs <- 3
    min.obs <- 1
    
    sites <- colnames(genotype)
    snp.counter <- 1
    genphen.data <- NULL
    temp.genphen.data <- NULL
    for(i in 1:ncol(genotype)) {
      gs <- unique(genotype[, i])
      
      # if not completely conserved position
      if(length(gs) != 1) {
        
        # loop through pairs of residues
        for(j in 1:(length(gs) - 1)) {
          for(k in (j+1):length(gs)) {
            
            # 1
            hit.j <- which(genotype[, i] == gs[j])
            
            # 0
            hit.k <- which(genotype[, i] == gs[k])
            
            if(length(hit.j) >= min.obs & length(hit.k) >= min.obs) {
              # create rows
              row.j <- data.frame(X = rep(x = 1, length = length(hit.j)),
                                  J = snp.counter, S = sites[i],
                                  Y = phenotype[hit.j], G = gs[j],
                                  stringsAsFactors = FALSE)
              row.k <- data.frame(X = rep(x = 0, length = length(hit.k)),
                                  J = snp.counter, S = sites[i],
                                  Y = phenotype[hit.k], G = gs[k],
                                  stringsAsFactors = FALSE)
              
              
              temp.genphen.data <- rbind(temp.genphen.data, 
                                         rbind(row.j, row.k))
              # genphen data
              if(snp.counter %% 50 == 1) {
                genphen.data <- rbind(genphen.data, temp.genphen.data)
                temp.genphen.data <- NULL
              }
              snp.counter <- snp.counter + 1
            }
          }
        }
      }
    }
    
    # remaining genphen data
    if(is.null(temp.genphen.data) == FALSE) {
      genphen.data <- rbind(genphen.data, temp.genphen.data)
      temp.genphen.data <- NULL
    }
    
    # if no genphen data, return null
    if(is.null(genphen.data) == TRUE
       ||length(genphen.data) == 0
       ||nrow(genphen.data) == 0) {
      return (NULL)
    }
    
    return (genphen.data)
  }
  
  
  # simple data (no-multicore)
  if(cores >= ncol(genotype)/2) {
    genphen.data <- getData(genotype, phenotype) 
  }
  else {
    # split data for mc
    is <- ceiling(seq(from = 1, to = ncol(genotype), length.out = cores + 1))
    is[1] <- 0
    
    
    # register cluster and go through all snps to extract stats, ppc, rpa
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    i <- NULL
    genphen.data <- (foreach(i = 1:(length(is) - 1)) %dopar% 
                       getData(genotype = genotype[, (is[i]+1):is[i+1]], 
                               phenotype = phenotype))
    # stop cluster
    parallel::stopCluster(cl = cl)
    doParallel::stopImplicitCluster()
    
    
    # rbind
    genphen.data <- do.call(rbind, genphen.data)
  }
  
  # form SNP ids
  genphen.data$J <- as.numeric(as.factor(paste(genphen.data$S, 
                                               genphen.data$J, 
                                               sep = '.')))
  
  return (genphen.data)
}




# Description:
# Get genotype summary
getGenSummary <- function(genotype) {
  # min.obs <- 3
  min.obs <- 1
  
  snp.counter <- 1
  gen.data <- NULL
  for(i in 1:ncol(genotype)) {
    gs <- unique(genotype[, i])
    
    # if not completely conserved position
    if(length(gs) != 1) {
      
      # loop through pairs of residues
      for(j in 1:(length(gs) - 1)) {
        for(k in (j+1):length(gs)) {
          
          row <- data.frame(g1 = gs[j], 
                            g0 = gs[k], 
                            site = i,
                            snp.id = snp.counter,
                            n1 = sum(genotype[, i] == gs[j]),
                            n0 = sum(genotype[, i] == gs[k]),
                            stringsAsFactors = FALSE)
          
          if(row$n1 >= min.obs & row$n0 >= min.obs) {
            # gen data
            gen.data <- rbind(gen.data, row)
            snp.counter <- snp.counter + 1
          }
        }
      }
    }
  }
  
  # if no genphen data, return null
  if(is.null(gen.data) || length(gen.data) == 0 || nrow(gen.data) == 0) {
    return (NULL)
  }
  
  return (gen.data)
}





# Description:
# Check the optional (...) parameters, if provided.
checkDotParameters <- function(...) {
  
  checkRpaSignificant <- function(rpa.significant) {
    if(length(rpa.significant) != 1) {
      stop("rpa.significant must be logical (TRUE/FALSE).")
    }
    
    if(is.logical(rpa.significant) == FALSE) {
      stop("rpa.significant must be logical.")
    }
  }
  
  checkAdaptDelta <- function(adapt_delta) {
    if(length(adapt_delta) != 1) {
      stop("adapt_delta must be in range (0, 1) (default = 0.8).")
    }
    
    if(is.numeric(adapt_delta) == FALSE) {
      stop("adapt_delta must be in range (0, 1)")
    }
    
    if(adapt_delta >= 1 | adapt_delta <= 0) {
      stop("adapt_delta must be in range (0, 1)")
    }
  }
  
  checkMaxTreedepth <- function(max_treedepth) {
    if(length(max_treedepth) != 1) {
      stop("max_treedepth is numeric parameter.")
    }
    
    if(is.numeric(max_treedepth) == FALSE) {
      stop("max_treedepth is numeric parameter.")
    }
    
    if(max_treedepth < 5) {
      stop("max_treedepth >= 5 (default = 10).")
    }
  }
  
  checkCvFold <- function(cv.fold) {
    if(length(cv.fold) != 1) {
      stop("cv.fold must be in range (0, 1) (default = 0.66).")
    }
    
    if(is.numeric(cv.fold) == FALSE) {
      stop("cv.fold must be in range (0, 1)")
    }
    
    if(cv.fold >= 1 | cv.fold <= 0) {
      stop("cv.fold must be in range (0, 1)")
    }
  }
  
  checkNtree <- function(ntree) {
    if(length(ntree) != 1) {
      stop("ntree is numeric parameter.")
    }
    
    if(is.numeric(ntree) == FALSE) {
      stop("ntree is numeric parameter.")
    }
    
    if(ntree < 100) {
      stop("ntree >= 100 (default = 500).")
    }
  }
  
  checkVerbose <- function(verbose) {
    if(length(verbose) != 1) {
      stop("verbose is a logical parameter.")
    }
    
    if(is.logical(verbose) == FALSE) {
      stop("verbose is a logical parameter.")
    }
  }
  
  checkRefresh <- function(refresh) {
    if(length(refresh) != 1) {
      stop("refresh is numeric parameter.")
    }
    
    if(is.numeric(refresh) == FALSE) {
      stop("refresh is a numeric parameter.")
    }
    
    return (refresh)
  }
  
  available.names <- c("adapt_delta", 
                       "max_treedepth", 
                       "ntree", 
                       "cv.fold", 
                       "rpa.significant",
                       "refresh",
                       "verbose")
  default.values <- list(adapt_delta = 0.9, 
                         max_treedepth = 10,
                         ntree = 1000, 
                         cv.fold = 0.66, 
                         rpa.significant = TRUE,
                         refresh = 250,
                         verbose = TRUE)
  
  # get the optional parameters
  dot.names <- names(list(...))
  
  if(length(dot.names) > 0) {
    if(any(dot.names %in% available.names) == FALSE) {
      wrong.names <- dot.names[!dot.names %in% available.names]
      stop(paste("Unknown optional parameter were provided! The following 
                 optional parameters are available:", dot.names, sep = ' '))
    }
  }
  
  # check each parameter
  for(p in dot.names) {
    if(is.null(list(...)[[p]]) || is.na(list(...)[[p]])) {
      stop(paste("optional parameter ", p, " can't be NULL", sep = ''))
    }
    if(p == "adapt_delta") {
      checkAdaptDelta(adapt_delta = list(...)[[p]])
      default.values[["adapt_delta"]] <- list(...)[[p]]
    }
    if(p == "rpa.significant") {
      checkRpaSignificant(rpa.significant = list(...)[[p]])
      default.values[["rpa.significant"]] <- list(...)[[p]]
    }
    if(p == "max_treedepth") {
      checkMaxTreedepth(max_treedepth = list(...)[[p]])
      default.values[["max_treedepth"]] <- list(...)[[p]]
    }
    if(p == "cv.fold") {
      checkCvFold(cv.fold = list(...)[[p]])
      default.values[["cv.fold"]] <- list(...)[[p]]
    }
    if(p == "ntree") {
      checkNtree(ntree = list(...)[[p]])
      default.values[["ntree"]] <- list(...)[[p]]
    }
    if(p == "refresh") {
      checkRefresh(refresh = list(...)[[p]])
      default.values[["refresh"]] <- list(...)[[p]]
    }
    if(p == "verbose") {
      checkVerbose(verbose = list(...)[[p]])
      default.values[["verbose"]] <- list(...)[[p]]
    }
  }
  
  return (default.values)
}





# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInput <- function(genotype, 
                       phenotype, 
                       phenotype.type, 
                       model.type,
                       mcmc.chains, 
                       mcmc.iterations, 
                       mcmc.warmup, 
                       cores, 
                       hdi.level, 
                       stat.learn.method, 
                       cv.iterations, 
                       rpa.iterations, 
                       with.stan.obj) {
  
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
  
  checkModelType <- function(model.type) {
    # CHECK: model.type
    if(length(model.type) != 1) {
      stop("model.type must be a string (default = 'continuous')")
    }
    
    if(!is.character(model.type)) {
      stop("model.type must be a string: 'univariate' or 'hierarchical'")
    }
    
    if(!model.type %in% c("univariate", "hierarchical")) {
      stop("phenotype.type must be a string: 'univariate' or 'hierarchical'")
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
  
  checkCores <- function(cores) {
    # CHECK: cores
    if(length(cores) != 1) {
      stop("cores is numeric parameter.")
    }
    
    if(is.numeric(cores) == FALSE) {
      stop("cores is numeric parameter.")
    }
    
    if(cores <= 0) {
      stop("cores is numeric parameter >=1.")
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
  
  checkRpa <- function(rpa.iterations) {
    if(is.numeric(rpa.iterations) == FALSE) {
      stop("rpa.iterations must be an integer")
    }
    
    if(length(rpa.iterations) != 1) {
      stop("rpa.iterations must be an integer.")
    }
    
    if(rpa.iterations < 0) {
      stop("rpa.iterations must be a positive integer (default = 10).")
    }
  }
  
  checkStanObj <- function(with.stan.obj) {
    if(is.logical(with.stan.obj) == FALSE) {
      stop("with.stan.obj must be logica.")
    }
    
    if(length(with.stan.obj) != 1) {
      stop("with.stan.obj must be logical (TRUE/FALSE).")
    }
  }
  
  if(is.null(genotype) | missing(genotype) |
     is.null(phenotype) | missing(phenotype) |
     is.null(phenotype.type) | missing(phenotype.type) |
     is.null(model.type) | missing(model.type) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.iterations) | missing(mcmc.iterations) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(cores) | missing(cores) |
     is.null(hdi.level) | missing(hdi.level) |
     is.null(stat.learn.method) | missing(stat.learn.method) |
     is.null(cv.iterations) | missing(cv.iterations) |
     is.null(rpa.iterations) | missing(rpa.iterations) |
     is.null(with.stan.obj) | missing(with.stan.obj)) {
    stop("arguments must be non-NULL/specified")
  }
  
  checkGenotypePhenotype(genotype = genotype, phenotype = phenotype)
  checkPhenotypeType(phenotype.type = phenotype.type)
  checkPhenotypeValidity(phenotype = phenotype, phenotype.type = phenotype.type)
  checkModelType(model.type = model.type)
  checkMcmcIterations(mcmc.iterations = mcmc.iterations)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkCores(cores = cores)
  checkHdi(hdi.level = hdi.level)
  checkMlMethod(stat.learn.method = stat.learn.method)
  checkCv(stat.learn.method = stat.learn.method, cv.iterations = cv.iterations)
  checkRpa(rpa.iterations = rpa.iterations)
  checkStanObj(with.stan.obj = with.stan.obj)
}





# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInputDiagnostics <- function(genotype, 
                                  diagnostic.points, 
                                  rf.trees){
  
  checkDiagnosticPoints <- function(diagnostic.points, rf.trees) {
    
    if(length(rf.trees) != 1) {
      stop("diagnostics.samples must be a number (default = 5,000).")
    }
    
    if(is.numeric(rf.trees) == FALSE) {
      stop("rf.trees must be a number (default = 5,000).")
    }
    
    if(rf.trees < 1000) {
      stop("rf.trees >= 1,000 accepted (default = 5,000).")
    }
  }
  
  if(is.null(diagnostic.points) | missing(diagnostic.points) | 
     is.null(rf.trees) | missing(rf.trees)) {
    stop("arguments must be non-NULL/specified")
  }
  
  checkDiagnosticPoints(diagnostic.points = diagnostic.points, 
                        rf.trees = rf.trees)
  
  
  if(length(diagnostic.points) > 0) {
    if(is.numeric(diagnostic.points) == FALSE) {
      stop("diagnostic.points must be a numeric vector.")
    }
    else {
      if(all(diagnostic.points %in% 1:ncol(genotype)) == FALSE) {
        stop("The diagnostic points must lie in the genotype space:",
             "[", 1, "-", ncol(genotype), "] \n", sep = '')
      }
    }
  }
}






# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInputPhyloBias <- function(input.kinship.matrix, 
                                genotype) {
  
  
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
# The classification is computed using random forest.
getRfCa <- function(genphen.data, 
                    cv.fold, 
                    cv.steps, 
                    hdi.level, 
                    ntree) {
  
  
  
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
      s <- sample(x = 1:rows, size = ceiling(x=cv.fold*rows), replace = TRUE)
      train <- data.frame(Y = Y[s], X = X[s], stringsAsFactors = FALSE)
      test <- data.frame(Y = Y[-s], X = X[-s], stringsAsFactors = FALSE)
      
      
      # only one type of predictor (no continous variable)
      ca <- NA
      kappa <- NA
      if(length(unique(train$X)) > 1) {
        
        
        # train classification model (try condition to avoid errors in case
        # only one-level predictor is train data)
        rf.out <- try(ranger::ranger(as.factor(Y) ~ X, data = train,
                                     num.trees = ntree), silent = TRUE)
        
        # test classification model
        pr <- stats::predict(object = rf.out, data = test)
        
        
        # compute classification accuracy (1 - classification error)
        ca<-sum(as.character(test$Y)==as.character(pr$predictions))/nrow(test)
        
        
        # compute kappa statistics
        kappa <- getKappa(real = as.character(test$Y),
                          predicted = as.character(pr$predictions),
                          aas = unique(Y))
      }
      
      return (list(ca = ca, kappa = kappa))
    }
    
    
    
    # cv.steps == 100
    if(I == 1) {
      
      
      # ca bootstrap
      ca.obj <- replicate(n = 100, expr = booter(Y = Y, X = X, 
                                                 ntree = 1000, 
                                                 cv.fold = cv.fold))
      
      
      # get cas and kappas
      cas <- unlist(ca.obj[which(rownames(ca.obj) == "ca"), ])
      kappas <- unlist(ca.obj[which(rownames(ca.obj) == "kappa"), ])
      
      
      # get X% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = cas, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      
      
      # build X% HDI for the kappas
      kappa.hdi <- getHdi(vec = kappas, hdi.level = hdi.level)
      kappa.L <- as.numeric(kappa.hdi[1])
      kappa.H <- as.numeric(kappa.hdi[2])
      
      return(list(ca = mean(cas, na.rm = TRUE),
                  ca.L = ca.L,
                  ca.H = ca.H,
                  kappa = mean(kappas, na.rm = TRUE),
                  kappa.L = kappa.L,
                  kappa.H = kappa.H,
                  I = 1))
    }
    
    
    # keeptrack
    old.list <- list(kappa.L = NA, kappa.H = NA, ca.L = NA, ca.H = NA)
    updated.list <- list(kappa.L = NA, kappa.H = NA, ca.L = NA, ca.H = NA)
    
    cas <- c()
    kappas <- c()
    for(i in 1:I) {
      
      
      # ca bootstrap
      ca.obj <- replicate(n = 100, expr = booter(Y = Y, X = X, 
                                                 ntree = 1000,
                                                 cv.fold = cv.fold))
      
      
      # get new cas and kappas
      new.ca <- unlist(ca.obj[which(rownames(ca.obj) == "ca"), ])
      new.kappa <- unlist(ca.obj[which(rownames(ca.obj) == "kappa"), ])
      
      
      # update parameters
      cas <- c(cas, new.ca)
      kappas <- c(kappas, new.kappa)
      
      
      # get X% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = cas, hdi.level = hdi.level)
      updated.list[["ca.L"]] <- as.numeric(ca.hdi[1])
      updated.list[["ca.H"]] <- as.numeric(ca.hdi[2])
      
      
      # build X% HDI for the kappas
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
                      kappa = mean(kappas, na.rm = TRUE),
                      kappa.L = updated.list[["kappa.L"]],
                      kappa.H = updated.list[["kappa.H"]],
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
                kappa = mean(kappas, na.rm = TRUE),
                kappa.L = updated.list[["kappa.L"]],
                kappa.H = updated.list[["kappa.H"]],
                I = I))
  }
  
  
  # general data
  general.data <- list(site = genphen.data$S[1],
                       snp.id = genphen.data$J[1],
                       g1 = unique(genphen.data$G[genphen.data$X == 1]),
                       g0 = unique(genphen.data$G[genphen.data$X == 0]),
                       n1 = sum(genphen.data$X == 1),
                       n0 = sum(genphen.data$X == 0),
                       y1 = sum(genphen.data$Y[genphen.data$X == 1]),
                       y0 = sum(genphen.data$Y[genphen.data$X == 0]))
  
  
  # subset predictor/response
  X <- genphen.data$Y
  Y <- genphen.data$X
  
  
  min.obs <- 3
  if(length(unique(X)) < 2 | all(table(Y) < min.obs)) {
    # pack dummy output, as at least a two-category predictor is needed
    # to run the incremental learning procedure
    class.obj <- list(ca = NA, ca.L = NA, ca.H = NA, 
                      kappa = NA, kappa.L = NA, kappa.H = NA, I = NA)
  }
  else {
    # run incremental CA learnings (100 steps in each iteration)
    class.obj <- getIncrementalLearning(Y = Y,
                                        X = X,
                                        cv.fold = cv.fold,
                                        ntree = ntree,
                                        I = ceiling(x = cv.steps/100),
                                        e = 0.01)
  }
  
  
  # collect the stats
  ca.out <- data.frame(site = general.data$site,
                       g1 = general.data$g1,
                       g0 = general.data$g0,
                       n1 = general.data$n1,
                       n0 = general.data$n0,
                       y1 = general.data$y1,
                       y0 = general.data$y0,
                       ca = class.obj$ca,
                       ca.L = class.obj$ca.L,
                       ca.H = class.obj$ca.H,
                       kappa = class.obj$kappa,
                       kappa.L = class.obj$kappa.L,
                       kappa.H = class.obj$kappa.H,
                       I = class.obj$I)
  
  
  return(ca.out)
}






# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using support vector machines.
getSvmCa <- function(genphen.data, 
                     cv.fold, 
                     cv.steps, 
                     hdi.level) {
  
  
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
      s <- sample(x = 1:rows, size = ceiling(x=cv.fold*rows), replace = TRUE)
      train <- data.frame(Y = Y[s], X = X[s], stringsAsFactors = FALSE)
      test <- data.frame(Y = Y[-s], X = X[-s], stringsAsFactors = FALSE)
      
      
      # only one type of predictor (no continous variable)
      ca <- NA
      kappa <- NA
      if(length(unique(train$X)) > 1) {
        
        
        # train classification model (try condition to avoid errors in case
        # only one-level predictor is train data)
        svm.out <- try(e1071::svm(as.factor(Y) ~ X, data = train,
                                  type = "C-classification"), silent = TRUE)
        
        # test classification model
        pr <- stats::predict(object = svm.out, data = test)
        
        
        # compute classification accuracy (1 - classification error)
        ca <- sum(as.character(test$Y)==as.character(pr))/nrow(test)
        
        
        # compute kappa statistics
        kappa <- getKappa(real = as.character(test$Y),
                          predicted = as.character(pr),
                          aas = unique(Y))
      }
      
      return (list(ca = ca, kappa = kappa))
    }
    
    
    
    # cv.steps == 100
    if(I == 1) {
      
      # ca bootstrap
      ca.obj<-replicate(n = 100, expr = booter(Y = Y, X = X, cv.fold = cv.fold))
      
      # get cas and kappas
      cas <- unlist(ca.obj[which(rownames(ca.obj) == "ca"), ])
      kappas <- unlist(ca.obj[which(rownames(ca.obj) == "kappa"), ])
      
      # get X% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = cas, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      
      # build X% HDI for the kappas
      kappa.hdi <- getHdi(vec = kappas, hdi.level = hdi.level)
      kappa.L <- as.numeric(kappa.hdi[1])
      kappa.H <- as.numeric(kappa.hdi[2])
      
      return(list(ca = mean(cas, na.rm = TRUE),
                  ca.L = ca.L,
                  ca.H = ca.H,
                  kappa = mean(kappas, na.rm = TRUE),
                  kappa.L = kappa.L,
                  kappa.H = kappa.H,
                  I = 1))
    }
    
    
    # keeptrack
    old.list <- list(kappa.L = NA, kappa.H = NA, ca.L = NA, ca.H = NA)
    updated.list <- list(kappa.L = NA, kappa.H = NA, ca.L = NA, ca.H = NA)
    
    cas <- c()
    kappas <- c()
    for(i in 1:I) {
      ca.obj<-replicate(n = 100, expr = booter(Y = Y, X = X, cv.fold = cv.fold))
      
      
      # get new cas and kappas
      new.ca <- unlist(ca.obj[which(rownames(ca.obj) == "ca"), ])
      new.kappa <- unlist(ca.obj[which(rownames(ca.obj) == "kappa"), ])
      
      
      # Update parameters
      cas <- c(cas, new.ca)
      kappas <- c(kappas, new.kappa)
      
      
      # get X% HDI for the classification accuracy
      ca.hdi <- getHdi(vec = cas, hdi.level = hdi.level)
      updated.list[["ca.L"]] <- as.numeric(ca.hdi[1])
      updated.list[["ca.H"]] <- as.numeric(ca.hdi[2])
      
      
      # build X% HDI for the kappas
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
                      kappa = mean(kappas, na.rm = TRUE),
                      kappa.L = updated.list[["kappa.L"]],
                      kappa.H = updated.list[["kappa.H"]],
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
                kappa = mean(kappas, na.rm = TRUE),
                kappa.L = updated.list[["kappa.L"]],
                kappa.H = updated.list[["kappa.H"]],
                I = I))
  }
  
  
  # general data
  general.data <- list(site = genphen.data$S[1],
                       snp.id = genphen.data$J[1],
                       g1 = unique(genphen.data$G[genphen.data$X == 1]),
                       g0 = unique(genphen.data$G[genphen.data$X == 0]),
                       n1 = sum(genphen.data$X == 1),
                       n0 = sum(genphen.data$X == 0),
                       y1 = sum(genphen.data$Y[genphen.data$X == 1]),
                       y0 = sum(genphen.data$Y[genphen.data$X == 0]))
  
  
  # subset predictor/response
  X <- genphen.data$Y
  Y <- genphen.data$X
  
  min.obs <- 3
  if(length(unique(X)) < 2 | all(table(Y) < min.obs)) {
    # pack dummy output, as at least a two-category predictor is needed
    # to run the incremental learning procedure
    class.obj <- list(ca = NA, ca.L = NA, ca.H = NA, kappa = NA, 
                      kappa.L = NA, kappa.H = NA, I = NA)
  }
  else {
    # run incremental CA learnings (100 steps in each iteration)
    class.obj <- getIncrementalLearning(Y = Y,
                                        X = X,
                                        cv.fold = cv.fold,
                                        I = ceiling(x = cv.steps/100),
                                        e = 0.01)
  }
  
  
  
  
  # collect the stats
  ca.out <- data.frame(site = general.data$site,
                       g1 = general.data$g1,
                       g0 = general.data$g0,
                       n1 = general.data$n1,
                       n0 = general.data$n0,
                       y1 = general.data$y1,
                       y0 = general.data$y0,
                       ca = class.obj$ca,
                       ca.L = class.obj$ca.L,
                       ca.H = class.obj$ca.H,
                       kappa = class.obj$kappa,
                       kappa.L = class.obj$kappa.L,
                       kappa.H = class.obj$kappa.H,
                       I = class.obj$I)
  
  
  return(ca.out)
}






# Description:
# Dummy CA output
getNoneCa <- function(genphen.data) {
  
  # general data
  general.data <- list(site = genphen.data$S[1],
                       snp.id = genphen.data$J[1],
                       g1 = unique(genphen.data$G[genphen.data$X == 1]),
                       g0 = unique(genphen.data$G[genphen.data$X == 0]),
                       n1 = sum(genphen.data$X == 1),
                       n0 = sum(genphen.data$X == 0),
                       y1 = sum(genphen.data$Y[genphen.data$X == 1]),
                       y0 = sum(genphen.data$Y[genphen.data$X == 0]))
  
  
  # collect the stats
  ca.out <- data.frame(site = general.data$site,
                       g1 = general.data$g1,
                       g0 = general.data$g0,
                       n1 = general.data$n1,
                       n0 = general.data$n0,
                       y1 = general.data$y1,
                       y0 = general.data$y0,
                       ca = NA, 
                       ca.L = NA, 
                       ca.H = NA,
                       kappa = NA, 
                       kappa.L = NA, 
                       kappa.H = NA,
                       I = NA)
  
  return (ca.out)
}






# Description:
# Bayesian univariate GLM with continuous outcome and a two-factor predictor.
runContU <- function(genphen.data, 
                     mcmc.chains, 
                     mcmc.iterations, 
                     mcmc.warmup, 
                     cores, 
                     hdi.level, 
                     model.stan, 
                     rpa.iterations, 
                     with.stan.obj, 
                     ...) {
  
  
  data.list <- list(X = genphen.data$X,
                    Y = genphen.data$Y,
                    Z = nrow(genphen.data),
                    site = genphen.data$S[1],
                    snp.id = genphen.data$J[1],
                    G = genphen.data$G)
  
  
  # get initial parameter values
  control <- list(adapt_delta = list(...)[["adapt_delta"]],
                 max_treedepth = list(...)[["max_treedepth"]])
  verbose <- list(...)[["verbose"]]
  refresh <- list(...)[["refresh"]]
  posterior <- rstan::sampling(object = model.stan,
                               data = data.list,
                               pars = c("alpha", "beta", "sigma", "nu"),
                               iter = mcmc.iterations,
                               warmup = mcmc.warmup,
                               chains = mcmc.chains,
                               cores = cores,
                               control = control,
                               verbose = verbose,
                               refresh = refresh)
  
  
  # if with.model == TRUE, keep the stan object
  stan.obj <- NULL
  if(with.stan.obj == TRUE) {
    stan.obj <- posterior
  }
  
  
  # general data
  general.data <- list(site = data.list$site,
                       snp.id = data.list$snp.id,
                       g1 = unique(data.list$G[data.list$X == 1]),
                       g0 = unique(data.list$G[data.list$X == 0]),
                       n1 = sum(data.list$X == 1),
                       n0 = sum(data.list$X == 0))
  
  
  # compute posterior summary
  hdi.L <- (1-hdi.level)/2
  hdi.H <- 1-(1-hdi.level)/2
  stats <- rstan::summary(object = posterior, 
                          pars = c("alpha", "beta", "sigma", "nu"), 
                          prob = c(hdi.L, hdi.H))$summary
  
  
  # convergence data
  convergence.out <- data.frame(stats[, c("Rhat", "n_eff")],
                                stringsAsFactors = FALSE)
  convergence.out$par <- rownames(convergence.out)
  convergence.out$site <- general.data$site
  convergence.out$g1 <- general.data$g1
  convergence.out$g0 <- general.data$g0
  convergence.out$n1 <- general.data$n1
  convergence.out$n0 <- general.data$n0
  rownames(convergence.out) <- NULL
  
  
  # posterior data
  posterior <- data.frame(rstan::extract(posterior))
  
  
  # ppc
  ppc.fun <- function(p, x, n) {
    return(p[1]+p[2]*x+p[3]*stats::rt(n = 1, df = p[4]))
  }
  ppc.1 <- apply(X = posterior[, c("alpha", "beta", "sigma", "nu")], MARGIN = 1, 
                 FUN = ppc.fun, x = 1, n = general.data$n1)
  ppc.0 <- apply(X = posterior[, c("alpha", "beta", "sigma", "nu")], MARGIN = 1, 
                 FUN = ppc.fun, x = 0, n = general.data$n0)
  bc <- getBhattacharyya(x = ppc.1, y = ppc.0)$bc
  ppc.1.hdi <- getHdi(vec = ppc.1, hdi.level = hdi.level)
  ppc.0.hdi <- getHdi(vec = ppc.0, hdi.level = hdi.level)
  
  
  # predicted vs real means
  ppc.out <- data.frame(site = general.data$site,
                        g1 = general.data$g1,
                        g0 = general.data$g0,
                        n1 = general.data$n1,
                        n0 = general.data$n0,
                        pred.mean.1 = mean(ppc.1),
                        pred.1.L = ppc.1.hdi[1],
                        pred.1.H = ppc.1.hdi[2],
                        pred.mean.0 = mean(ppc.0),
                        pred.0.L = ppc.0.hdi[1],
                        pred.0.H = ppc.0.hdi[2],
                        real.mean.1 = mean(data.list$Y[data.list$X == 1]),
                        real.mean.0 = mean(data.list$Y[data.list$X == 0]))
  
  
  # collect the stats
  stats.out <- data.frame(site = general.data$site,
                          g1 = general.data$g1,
                          g0 = general.data$g0,
                          n1 = general.data$n1,
                          n0 = general.data$n0,
                          beta.mean = stats["beta", "mean"],
                          beta.se = stats["beta", "se_mean"],
                          beta.sd = stats["beta", "sd"],
                          beta.L = stats["beta", paste(hdi.L*100,"%",sep='')],
                          beta.H = stats["beta", paste(hdi.H*100,"%",sep='')],
                          alpha.mean = stats["alpha", "mean"],
                          alpha.se = stats["alpha", "se_mean"],
                          alpha.sd = stats["alpha", "sd"],
                          alpha.L = stats["alpha", paste(hdi.L*100,"%",sep='')],
                          alpha.H = stats["alpha", paste(hdi.H*100,"%",sep='')],
                          sigma.mean = stats["sigma", "mean"],
                          sigma.se = stats["sigma", "se_mean"],
                          sigma.sd = stats["sigma", "sd"],
                          sigma.L = stats["sigma", paste(hdi.L*100,"%",sep='')],
                          sigma.H = stats["sigma", paste(hdi.H*100,"%",sep='')],
                          nu.mean = stats["nu", "mean"],
                          nu.se = stats["nu", "se_mean"],
                          nu.sd = stats["nu", "sd"],
                          nu.L = stats["nu", paste(hdi.L*100,"%",sep='')],
                          nu.H = stats["nu", paste(hdi.H*100,"%",sep='')],
                          bc = bc)
  
  # special case for RPA
  rpa.out <- NULL
  if(list(...)[["rpa.significant"]] == TRUE) {
    if((stats.out$beta.L <= 0 & stats.out$beta.H >= 0) == TRUE) {
      rpa.iterations <- 0
      rpa.out <- data.frame(site = general.data$site,
                            g1 = general.data$g1,
                            g0 = general.data$g0,
                            n1 = general.data$n1,
                            n0 = general.data$n0,
                            rpa.power.error = NA,
                            rpa.sign.error = NA,
                            rpa.beta.mean = NA,
                            rpa.beta.sd = NA,
                            rpa.N = NA,
                            stringsAsFactors = FALSE)
    }
  }
  
  if(rpa.iterations > 0) {
    rpa.out <- getRpaCont(posterior = posterior,
                          beta.mean = stats["beta", "mean"],
                          site = general.data$site,
                          n1 = general.data$n1,
                          n0 = general.data$n0,
                          g1 = general.data$g1, 
                          g0 = general.data$g0,
                          hdi.level = hdi.level, 
                          rpa.iterations = rpa.iterations,
                          model.stan = model.stan,
                          mcmc.iterations = mcmc.iterations,
                          mcmc.warmup = mcmc.warmup,
                          mcmc.chains = mcmc.chains,
                          adapt_delta = list(...)[["adapt_delta"]],
                          max_treedepth = list(...)[["max_treedepth"]])
  }
  
  # return
  return (list(statistics.out = stats.out,
               convergence.out = convergence.out,
               rpa.out = rpa.out,
               ppc.out = ppc.out,
               stan.obj = stan.obj))
}





# Description:
# Bayesian hierarchical GLM with continuous outcome and a two-factor predictor.
runContH <- function(genphen.data, 
                     mcmc.chains, 
                     mcmc.iterations, 
                     mcmc.warmup, 
                     cores, 
                     hdi.level, 
                     model.stan,
                     rpa.iterations, 
                     with.stan.obj, 
                     ...) {
  
  
  # Description:
  # Collection of results, ppc and rpa in case of hierarchical analysis
  getResults <- function(j, data.list, general.data, posterior, 
                         stats, hdi.level, mcmc.chains, mcmc.iterations, 
                         mcmc.warmup, model.stan, rpa.iterations, ...) {
    
    # ppc
    ppc.fun <- function(p, x, n) {
      return(p[1]+p[2]*x+p[3]*stats::rt(n = 1, df = p[4]))
    }
    
    # extract their stats including ppc, rpa
    ppc.1 <- apply(X = posterior, MARGIN = 1, FUN = ppc.fun, 
                   x = 1, n = general.data$n1)
    ppc.0 <- apply(X = posterior, MARGIN = 1, FUN = ppc.fun, 
                   x = 0, n = general.data$n0)
    bc <- getBhattacharyya(x = ppc.1, y = ppc.0)$bc
    ppc.1.hdi <- getHdi(vec = ppc.1, hdi.level = hdi.level)
    ppc.0.hdi <- getHdi(vec = ppc.0, hdi.level = hdi.level)
    
    
    # predicted vs real means
    ppc.out <- data.frame(site = general.data$site,
                          g1 = general.data$g1,
                          g0 = general.data$g0,
                          n1 = general.data$n1,
                          n0 = general.data$n0,
                          pred.mean.1 = mean(ppc.1),
                          pred.1.L = ppc.1.hdi[1],
                          pred.1.H = ppc.1.hdi[2],
                          pred.mean.0 = mean(ppc.0),
                          pred.0.L = ppc.0.hdi[1],
                          pred.0.H = ppc.0.hdi[2],
                          real.mean.1 = mean(data.list$Y[data.list$X == 1
                                                         & data.list$S == j]),
                          real.mean.0 = mean(data.list$Y[data.list$X == 0
                                                         & data.list$S == j]))
    
    
    # collect the stats
    a.key <- paste("alpha[", j, "]", sep = '')
    b.key <- paste("beta[", j, "]", sep = '')
    s.key <- "sigma"
    n.key <- "nu"
    stats.out <- data.frame(site = general.data$site,
                            g1 = general.data$g1,
                            g0 = general.data$g0,
                            n1 = general.data$n1,
                            n0 = general.data$n0,
                            beta.mean = stats[b.key, "mean"],
                            beta.se = stats[b.key, "se_mean"],
                            beta.sd = stats[b.key, "sd"],
                            beta.L = stats[b.key, paste(hdi.L*100,"%",sep='')],
                            beta.H = stats[b.key, paste(hdi.H*100,"%",sep='')],
                            alpha.mean = stats[a.key, "mean"],
                            alpha.se = stats[a.key, "se_mean"],
                            alpha.sd = stats[a.key, "sd"],
                            alpha.L = stats[a.key, paste(hdi.L*100,"%",sep='')],
                            alpha.H = stats[a.key, paste(hdi.H*100,"%",sep='')],
                            sigma.mean = stats[s.key, "mean"],
                            sigma.se = stats[s.key, "se_mean"],
                            sigma.sd = stats[s.key, "sd"],
                            sigma.L = stats[s.key, paste(hdi.L*100,"%",sep='')],
                            sigma.H = stats[s.key, paste(hdi.H*100,"%",sep='')],
                            nu.mean = stats[n.key, "mean"],
                            nu.se = stats[n.key, "se_mean"],
                            nu.sd = stats[n.key, "sd"],
                            nu.L = stats[n.key, paste(hdi.L*100,"%",sep='')],
                            nu.H = stats[n.key, paste(hdi.H*100,"%",sep='')],
                            bc = bc)
    
    # special case for RPA
    rpa.out <- NULL
    if(list(...)[["rpa.significant"]] == TRUE) {
      if((stats.out$beta.L <= 0 & stats.out$beta.H >= 0) == TRUE) {
        rpa.iterations <- 0
        rpa.out <- data.frame(site = general.data$site,
                              g1 = general.data$g1,
                              g0 = general.data$g0,
                              n1 = general.data$n1,
                              n0 = general.data$n0,
                              rpa.power.error = NA,
                              rpa.sign.error = NA,
                              rpa.beta.mean = NA,
                              rpa.beta.sd = NA,
                              rpa.N = NA,
                              stringsAsFactors = FALSE)
      }
    }
    
    
    if(rpa.iterations > 0) {
      # create a temporary posterior similar to the one in the univariate case
      a.key <- paste("alpha", j, sep = '.')
      b.key <- paste("beta", j, sep = '.')
      s.key <- "sigma"
      n.key <- "nu"
      colnames(posterior) <- c("alpha", "beta", "sigma", "nu")
      
      rpa.out <- getRpaCont(posterior = posterior,
                            beta.mean = stats.out$beta.mean,
                            site = general.data$site,
                            g1 = general.data$g1,
                            g0 = general.data$g0,
                            n1 = general.data$n1,
                            n0 = general.data$n0,
                            hdi.level = hdi.level, 
                            rpa.iterations = rpa.iterations,
                            model.stan = model.stan,
                            mcmc.iterations = mcmc.iterations,
                            mcmc.warmup = mcmc.warmup,
                            mcmc.chains = mcmc.chains,
                            adapt_delta = list(...)[["adapt_delta"]],
                            max_treedepth = list(...)[["max_treedepth"]])
    }
    
    # return
    return (list(statistics.out = stats.out,
                 rpa.out = rpa.out,
                 ppc.out = ppc.out))
  }
  
  
  data.list <- list(X = genphen.data$X,
                    Y = genphen.data$Y,
                    Z = nrow(genphen.data),
                    S = genphen.data$J, 
                    S_N = max(genphen.data$J), 
                    G = genphen.data$G,
                    site = genphen.data$S)
  
  
  # get initial parameter values
  control <- list(adapt_delta = list(...)[["adapt_delta"]],
                 max_treedepth = list(...)[["max_treedepth"]])
  refresh <- list(...)[["refresh"]]
  verbose <- list(...)[["verbose"]]
  posterior <- rstan::sampling(object = model.stan,
                               data = data.list,
                               iter = mcmc.iterations,
                               warmup = mcmc.warmup,
                               chains = mcmc.chains,
                               cores = cores,
                               control = control,
                               verbose = verbose,
                               refresh = refresh)
  
  
  # if with.model == TRUE, keep the stan object
  stan.obj <- NULL
  if(with.stan.obj == TRUE) {
    stan.obj <- posterior
  }
  
  
  # get general data
  general.data <- c()
  for(j in 1:max(data.list$S)) {
    r <- data.frame(site = unique(data.list$site[data.list$S == j]),
                    S = j,
                    g1 = data.list$G[data.list$S == j & data.list$X == 1][1],
                    g0 = data.list$G[data.list$S == j & data.list$X == 0][1],
                    n1 = sum(data.list$S == j & data.list$X == 1),
                    n0 = sum(data.list$S == j & data.list$X == 0),
                    stringsAsFactors = FALSE)
    general.data <- rbind(general.data, r)
  }
  rm(r, j)
  
  
  # compute posterior summary
  hdi.L <- (1-hdi.level)/2
  hdi.H <- 1-(1-hdi.level)/2
  stats <- rstan::summary(object = posterior, 
                          pars = c("alpha", "beta", "sigma", "nu",
                                   "mu_alpha", "mu_beta",
                                   "sigma_alpha", "sigma_beta",
                                   "nu_alpha", "nu_beta"), 
                          prob = c(hdi.L, hdi.H))$summary
  
  
  # convergence data
  convergence.out <- c()
  c.row <- data.frame(stats[c("sigma", "nu", 
                              "mu_alpha", "mu_beta", 
                              "sigma_alpha", "sigma_beta", 
                              "nu_alpha", "nu_beta"), 
                            c("Rhat", "n_eff")], stringsAsFactors = FALSE)
  c.row$par <- c("sigma", "nu", 
                 "mu_alpha", "mu_beta", 
                 "sigma_alpha", "sigma_beta",
                 "nu_alpha", "nu_beta")
  c.row$site <- ''
  c.row$g1 <- ''
  c.row$g0 <- ''
  c.row$n1 <- ''
  c.row$n0 <- ''
  convergence.out <- rbind(convergence.out, c.row)
  
  for(j in 1:nrow(general.data)) {
    c.row <- data.frame(stats[paste(c("alpha", "beta"), "[", j, "]", sep = ''), 
                              c("Rhat", "n_eff")], stringsAsFactors = FALSE)
    c.row$par <- c("alpha", "beta")
    c.row$site <- general.data$site[j]
    c.row$g1 <- general.data$g1[j]
    c.row$g0 <- general.data$g0[j]
    c.row$n1 <- general.data$n1[j]
    c.row$n0 <- general.data$n0[j]
    convergence.out <- rbind(convergence.out, c.row)
  }
  
  
  # posterior data
  posterior <- data.frame(rstan::extract(posterior))
  
  
  # special case for RPA, compile univariate stan model
  if(rpa.iterations > 0) {
    cat("======== Compiling RPA Model ======== \n")
    model.stan <- compileModel(phenotype.type = "continuous", 
                                   model.type = "univariate")
  }
  
  
  # register cluster and go through all snps to extract stats, ppc, rpa
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  results <- (foreach(j = 1:nrow(general.data),
                      .export = c("runContU", "getHdi", "getGenphenData",
                                  "getBhattacharyya", "getRpaCont"), 
                      .packages = c("rstan")) %dopar% 
                getResults(j = j, 
                           data.list = data.list, 
                           general.data = general.data[j, ], 
                           posterior = posterior[, c(paste(c("alpha", "beta"), 
                                                           j, sep = '.'), 
                                                     "sigma", "nu")], 
                           stats = stats, 
                           hdi.level = hdi.level,
                           mcmc.chains = mcmc.chains, 
                           mcmc.iterations = mcmc.iterations, 
                           mcmc.warmup = mcmc.warmup, 
                           model.stan = model.stan, 
                           rpa.iterations = rpa.iterations,
                           adapt_delta = list(...)[["adapt_delta"]],
                           max_treedepth = list(...)[["max_treedepth"]],
                           rpa.significant = list(...)[["rpa.significant"]]))
  # stop cluster
  parallel::stopCluster(cl = cl)
  doParallel::stopImplicitCluster()
  
  
  # rbind results
  getO <- function(x, y) {
    return(x[[y]])
  }
  stats.out <- do.call(rbind, lapply(X = results,FUN = getO,y="statistics.out"))
  if(rpa.iterations == 0) {
    rpa.out <- NULL
  }
  else {
    rpa.out <- do.call(rbind, lapply(X = results, FUN = getO, y="rpa.out"))
  }
  ppc.out <- do.call(rbind, lapply(X = results, FUN = getO, y="ppc.out"))
  
  
  # return
  return (list(statistics.out = stats.out,
               convergence.out = convergence.out,
               rpa.out = rpa.out,
               ppc.out = ppc.out,
               stan.obj = stan.obj))
}




# Description:
# Bayesian GLM with dichotomous outcome and a two-factor predictor.
runDichU <- function(genphen.data,
                     mcmc.chains, 
                     mcmc.iterations, 
                     mcmc.warmup, 
                     cores, 
                     hdi.level, 
                     model.stan, 
                     rpa.iterations, 
                     with.stan.obj, 
                     ...) {
  
  
  genphen.data$N <- 1
  d.N <- stats::aggregate(formula = N~X+S+J+G, data = genphen.data, FUN = sum)
  d.Y <- stats::aggregate(formula = Y~X+S+J+G, data = genphen.data, FUN = sum)
  d <- merge(x = d.N, y = d.Y, by = c("X", "S", "J", "G"))
  data.list <- list(X = d$X, Y = d$Y, N = d$N, Z = nrow(d), site = d$S[1], 
                    snp.id = d$J[1], G = d$G)
  rm(genphen.data, d.N, d.Y, d)
  
  
  # get initial parameter values
  control <- list(adapt_delta = list(...)[["adapt_delta"]],
                  max_treedepth = list(...)[["max_treedepth"]])
  refresh <- list(...)[["refresh"]]
  verbose <- list(...)[["verbose"]]
  posterior <- rstan::sampling(object = model.stan,
                               data = data.list,
                               pars = c("alpha", "beta"),
                               iter = mcmc.iterations,
                               warmup = mcmc.warmup,
                               chains = mcmc.chains,
                               cores = cores,
                               control = control,
                               verbose = verbose,
                               refresh = refresh)
  
  
  # if with.model == TRUE, keep the stan object
  stan.obj <- NULL
  if(with.stan.obj == TRUE) {
    stan.obj <- posterior
  }
  
  
  # general data
  general.data <- list(site = data.list$site,
                       snp.id = data.list$snp.id,
                       g1 = data.list$G[data.list$X == 1],
                       g0 = data.list$G[data.list$X == 0],
                       n1 = data.list$N[data.list$X == 1],
                       n0 = data.list$N[data.list$X == 0],
                       y1 = data.list$Y[data.list$X == 1],
                       y0 = data.list$Y[data.list$X == 0])
  
  
  # compute posterior summary
  hdi.L <- (1-hdi.level)/2
  hdi.H <- 1-(1-hdi.level)/2
  stats <- rstan::summary(object = posterior,
                          pars = c("alpha", "beta"),
                          prob = c(hdi.L, hdi.H))$summary
  
  
  # convergence data
  convergence.out <- data.frame(stats[, c("Rhat", "n_eff")],
                                stringsAsFactors = FALSE)
  convergence.out$par <- rownames(convergence.out)
  convergence.out$site <- general.data$site
  convergence.out$g1 <- general.data$g1
  convergence.out$g0 <- general.data$g0
  convergence.out$n1 <- general.data$n1
  convergence.out$n0 <- general.data$n0
  rownames(convergence.out) <- NULL
  
  
  # posterior data
  posterior <- data.frame(rstan::extract(posterior))
  
  
  # ppc
  ppc.fun <- function(p, x, n) {
    return(stats::rbinom(n = 1, prob = 1/(1+exp(-(p[1]+p[2]*x))), size = n))
  }
  
  ppc.1 <- apply(X = posterior[, c("alpha", "beta")], MARGIN = 1, 
                 FUN = ppc.fun, x = 1, n = general.data$n1)
  ppc.0 <- apply(X = posterior[, c("alpha", "beta")], MARGIN = 1, 
                 FUN = ppc.fun, x = 0, n = general.data$n0)
  bc <- getBhattacharyya(x = ppc.1/general.data$n1, 
                         y = ppc.0/general.data$n0)$bc
  ppc.1.hdi <- getHdi(vec = ppc.1/general.data$n1, 
                      hdi.level = hdi.level)
  ppc.0.hdi <- getHdi(vec = ppc.0/general.data$n0, 
                      hdi.level = hdi.level)
  
  
  # predicted vs real means
  ppc.out <- data.frame(site = general.data$site,
                        g1 = general.data$g1,
                        g0 = general.data$g0,
                        n1 = general.data$n1,
                        n0 = general.data$n0,
                        pred.mean.1 = mean(ppc.1/general.data$n1),
                        pred.1.L = ppc.1.hdi[1],
                        pred.1.H = ppc.1.hdi[2],
                        pred.mean.0 = mean(ppc.0/general.data$n0),
                        pred.0.L = ppc.0.hdi[1],
                        pred.0.H = ppc.0.hdi[2],
                        real.mean.1 = mean(data.list$Y[data.list$X == 1]),
                        real.mean.0 = mean(data.list$Y[data.list$X == 0]))
  
  
  # collect the stats
  stats.out <- data.frame(site = general.data$site,
                          g1 = general.data$g1,
                          g0 = general.data$g0,
                          n1 = general.data$n1,
                          n0 = general.data$n0,
                          beta.mean = stats["beta", "mean"],
                          beta.se = stats["beta", "se_mean"],
                          beta.sd = stats["beta", "sd"],
                          beta.L = stats["beta", paste(hdi.L*100,"%",sep='')],
                          beta.H = stats["beta", paste(hdi.H*100,"%",sep='')],
                          alpha.mean = stats["alpha", "mean"],
                          alpha.se = stats["alpha", "se_mean"],
                          alpha.sd = stats["alpha", "sd"],
                          alpha.L = stats["alpha", paste(hdi.L*100,"%",sep='')],
                          alpha.H = stats["alpha", paste(hdi.H*100,"%",sep='')],
                          bc = bc)
  
  
  # special case for RPA
  rpa.out <- NULL
  if(list(...)[["rpa.significant"]] == TRUE) {
    if((stats.out$beta.L <= 0 & stats.out$beta.H >= 0) == TRUE) {
      rpa.iterations <- 0
      rpa.out <- data.frame(site = general.data$site,
                            g1 = general.data$g1,
                            g0 = general.data$g0,
                            n1 = general.data$n1,
                            n0 = general.data$n0,
                            rpa.power.error = NA,
                            rpa.sign.error = NA,
                            rpa.beta.mean = NA,
                            rpa.beta.sd = NA,
                            rpa.N = NA,
                            stringsAsFactors = FALSE)
    }
  }
  
  if(rpa.iterations > 0) {
    colnames(posterior) <- c("alpha", "beta")
    rpa.out <- getRpaDich(posterior = posterior,
                          beta.mean = mean(posterior$beta),
                          site = general.data$site,
                          g1 = general.data$g1,
                          g0 = general.data$g0,
                          n1 = general.data$n1,
                          n0 = general.data$n0,
                          hdi.level = hdi.level, 
                          rpa.iterations = rpa.iterations,
                          model.stan = model.stan,
                          mcmc.iterations = mcmc.iterations,
                          mcmc.warmup = mcmc.warmup,
                          mcmc.chains = mcmc.chains,
                          adapt_delta = list(...)[["adapt_delta"]],
                          max_treedepth = list(...)[["max_treedepth"]])
  }
  
  
  # return
  return (list(statistics.out = stats.out,
               convergence.out = convergence.out,
               rpa.out = rpa.out,
               ppc.out = ppc.out,
               stan.obj = stan.obj))
}





# Description:
# Bayesian GLM with dichotomous outcome and a two-factor predictor.
runDichH <- function(genphen.data, 
                     mcmc.chains, 
                     mcmc.iterations, 
                     mcmc.warmup, 
                     cores, 
                     hdi.level, 
                     model.stan, 
                     rpa.iterations, 
                     with.stan.obj, 
                     ...) {
  
  
  # Description:
  # Collection of results, ppc and rpa in case of hierarchical analysis
  getResults <- function(j, data.list, general.data, posterior, 
                         stats, hdi.level, mcmc.chains, mcmc.iterations, 
                         mcmc.warmup, model.stan, rpa.iterations, ...) {
    
    # ppc
    ppc.fun <- function(p, x, n) {
      # y <- stats::rnorm(n = 1, mean = p[1]+p[2]*x, sd = p[3])
      return(stats::rbinom(n = 1, size = n, prob = 1/(1+exp(-(p[1]+p[2]*x)))))
    }
    
    ppc.1 <- apply(X = posterior, MARGIN = 1, FUN = ppc.fun, 
                   x = 1, n = general.data$n1)
    ppc.0 <- apply(X = posterior, MARGIN = 1, FUN = ppc.fun, 
                   x = 0, n = general.data$n0)
    bc <- getBhattacharyya(x=ppc.1/general.data$n1, y=ppc.0/general.data$n0)$bc
    ppc.1.hdi <- getHdi(vec = ppc.1/general.data$n1, hdi.level = hdi.level)
    ppc.0.hdi <- getHdi(vec = ppc.0/general.data$n0, hdi.level = hdi.level)
    
    
    # predicted vs real means
    ppc.out <- data.frame(site = general.data$site,
                          g1 = general.data$g1,
                          g0 = general.data$g0,
                          n1 = general.data$n1,
                          n0 = general.data$n0,
                          pred.mean.1 = mean(ppc.1/general.data$n1),
                          pred.1.L = ppc.1.hdi[1],
                          pred.1.H = ppc.1.hdi[2],
                          pred.mean.0 = mean(ppc.0/general.data$n0),
                          pred.0.L = ppc.0.hdi[1],
                          pred.0.H = ppc.0.hdi[2],
                          real.mean.1 = mean(data.list$Y[data.list$X == 1
                                                         & data.list$S == j]),
                          real.mean.0 = mean(data.list$Y[data.list$X == 0
                                                         & data.list$S == j]))
    
    
    # collect the stats
    a.key <- paste("alpha[", j, "]", sep = '')
    b.key <- paste("beta[", j, "]", sep = '')
    stats.out <- data.frame(site = general.data$site,
                            g1 = general.data$g1,
                            g0 = general.data$g0,
                            n1 = general.data$n1,
                            n0 = general.data$n0,
                            beta.mean = stats[b.key, "mean"],
                            beta.se = stats[b.key, "se_mean"],
                            beta.sd = stats[b.key, "sd"],
                            beta.L = stats[b.key,paste(hdi.L*100,"%",sep='')],
                            beta.H = stats[b.key,paste(hdi.H*100,"%",sep='')],
                            alpha.mean = stats[a.key, "mean"],
                            alpha.se = stats[a.key, "se_mean"],
                            alpha.sd = stats[a.key, "sd"],
                            alpha.L =stats[a.key,paste(hdi.L*100,"%",sep='')],
                            alpha.H =stats[a.key,paste(hdi.H*100,"%",sep='')],
                            bc = bc)
    
    
    # special case for RPA
    rpa.out <- NULL
    if(list(...)[["rpa.significant"]] == TRUE) {
      if((stats.out$beta.L <= 0 & stats.out$beta.H >= 0) == TRUE) {
        rpa.iterations <- 0
        rpa.out <- data.frame(site = general.data$site,
                              g1 = general.data$g1,
                              g0 = general.data$g0,
                              n1 = general.data$n1,
                              n0 = general.data$n0,
                              rpa.power.error = NA,
                              rpa.sign.error = NA,
                              rpa.beta.mean = NA,
                              rpa.beta.sd = NA,
                              rpa.N = NA,
                              stringsAsFactors = FALSE)
      }
    }
    
    if(rpa.iterations > 0) {
      colnames(posterior) <- c("alpha", "beta")
      rpa.out <- getRpaDich(posterior = posterior,
                            beta.mean = mean(posterior$beta),
                            site = general.data$site,
                            g1 = general.data$g1,
                            g0 = general.data$g0,
                            n1 = general.data$n1,
                            n0 = general.data$n0,
                            hdi.level = hdi.level, 
                            rpa.iterations = rpa.iterations,
                            model.stan = model.stan,
                            mcmc.iterations = mcmc.iterations,
                            mcmc.warmup = mcmc.warmup,
                            mcmc.chains = mcmc.chains,
                            adapt_delta = list(...)[["adapt_delta"]],
                            max_treedepth = list(...)[["max_treedepth"]])
    }
    
    
    # return
    return (list(statistics.out = stats.out,
                 rpa.out = rpa.out,
                 ppc.out = ppc.out))
  }
  
  
  genphen.data$N <- 1
  d.N <- stats::aggregate(formula = N~X+S+J+G, data = genphen.data, FUN = sum)
  d.Y <- stats::aggregate(formula = Y~X+S+J+G, data = genphen.data, FUN = sum)
  d <- merge(x = d.N, y = d.Y, by = c("X", "S", "J", "G"))
  data.list <- list(X = d$X, Y = d$Y, N = d$N, Z = nrow(d), site = d$S, 
                    S = d$J, S_N = max(d$J), G = d$G)
  rm(genphen.data, d.N, d.Y, d)
  
  
  # get initial parameter values
  control <- list(adapt_delta = list(...)[["adapt_delta"]],
                  max_treedepth = list(...)[["max_treedepth"]])
  refresh <- list(...)[["refresh"]]
  verbose <- list(...)[["verbose"]]
  posterior <- rstan::sampling(object = model.stan,
                               data = data.list,
                               iter = mcmc.iterations,
                               warmup = mcmc.warmup,
                               chains = mcmc.chains,
                               cores = cores,
                               control = control,
                               verbose = verbose,
                               refresh = refresh)
  
  # if with.model == TRUE, keep the stan object
  stan.obj <- NULL
  if(with.stan.obj == TRUE) {
    stan.obj <- posterior
  }
  
  
  # get general data
  general.data <- c()
  for(j in 1:max(data.list$S)) {
    r <- data.frame(site = data.list$site[data.list$S == j & data.list$X == 1],
                    S = data.list$S[data.list$S == j & data.list$X == 1],
                    g1 = data.list$G[data.list$S == j & data.list$X == 1],
                    g0 = data.list$G[data.list$S == j & data.list$X == 0],
                    n1 = data.list$N[data.list$S == j & data.list$X == 1],
                    n0 = data.list$N[data.list$S == j & data.list$X == 0],
                    y1 = data.list$Y[data.list$S == j & data.list$X == 1],
                    y0 = data.list$Y[data.list$S == j & data.list$X == 0],
                    stringsAsFactors = FALSE)
    general.data <- rbind(general.data, r)
  }
  rm(r, j)
  
  
  # compute posterior summary
  hdi.L <- (1-hdi.level)/2
  hdi.H <- 1-(1-hdi.level)/2
  stats <- rstan::summary(object = posterior, 
                          pars = c("alpha", "beta", 
                                   "mu_alpha", "mu_beta",
                                   "sigma_alpha", "sigma_beta",
                                   "nu_alpha", "nu_beta"), 
                          prob = c(hdi.L, hdi.H))$summary
  
  
  
  # convergence data
  convergence.out <- c()
  c.row <- data.frame(stats[c("mu_alpha", "mu_beta", 
                              "sigma_alpha", "sigma_beta",
                              "nu_alpha", "nu_beta"), 
                            c("Rhat", "n_eff")], 
                      stringsAsFactors = FALSE)
  c.row$par <- c("mu_alpha", "mu_beta", 
                 "sigma_alpha", "sigma_beta", 
                 "nu_alpha", "nu_beta")
  c.row$site <- ''
  c.row$g1 <- ''
  c.row$g0 <- ''
  c.row$n1 <- ''
  c.row$n0 <- ''
  convergence.out <- rbind(convergence.out, c.row)
  for(j in 1:nrow(general.data)) {
    c.row <- data.frame(stats[paste(c("alpha", "beta"), "[", j, "]", sep = ''), 
                              c("Rhat", "n_eff")], stringsAsFactors = FALSE)
    c.row$par <- c("alpha", "beta")
    c.row$site <- general.data$site[j]
    c.row$g1 <- general.data$g1[j]
    c.row$g0 <- general.data$g0[j]
    c.row$n1 <- general.data$n1[j]
    c.row$n0 <- general.data$n0[j]
    convergence.out <- rbind(convergence.out, c.row)
  }
  
  
  # posterior data
  posterior <- data.frame(rstan::extract(posterior))
  
  
  # special case for RPA, compile univariate stan model
  if(rpa.iterations > 0) {
    cat("======== Compiling RPA Model ======== \n")
    model.stan <- compileModel(phenotype.type = "dichotomous", 
                                   model.type = "univariate")
  }
  
  
  # register cluster and go through all snps to extract stats, ppc, rpa
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  results <- (foreach(j = 1:nrow(general.data),
                      .export = c("runDichU", "getHdi", "getGenphenData",
                                  "getBhattacharyya", "getRpaDich"),
                      .packages = c("rstan")) %dopar%
                getResults(j = j,
                           data.list = data.list,
                           general.data = general.data[j, ],
                           posterior = posterior[, c(paste(c("alpha", "beta"),
                                                           j, sep = '.'))],
                           stats = stats,
                           hdi.level = hdi.level,
                           mcmc.chains = mcmc.chains,
                           mcmc.iterations = mcmc.iterations,
                           mcmc.warmup = mcmc.warmup,
                           model.stan = model.stan,
                           rpa.iterations = rpa.iterations,
                           adapt_delta = list(...)[["adapt_delta"]],
                           max_treedepth = list(...)[["max_treedepth"]],
                           rpa.significant = list(...)[["rpa.significant"]]))
  # stop cluster
  parallel::stopCluster(cl = cl)
  doParallel::stopImplicitCluster()
  
  
  # rbind results
  getO <- function(x, y) {
    return(x[[y]])
  }
  stats.out <- do.call(rbind, lapply(X = results,FUN = getO,y="statistics.out"))
  if(rpa.iterations == 0) {
    rpa.out <- NULL
  }
  else {
    rpa.out <- do.call(rbind, lapply(X = results, FUN = getO, y="rpa.out"))
  }
  ppc.out <- do.call(rbind, lapply(X = results, FUN = getO, y="ppc.out"))
  
  
  # return
  return (list(statistics.out = stats.out,
               convergence.out = convergence.out,
               rpa.out = rpa.out,
               ppc.out = ppc.out,
               stan.obj = stan.obj))
}





# Description:
# Given a genotype dataset containing SNPs (columns) and N individuals (rows),
# the procedure computes a NxN kinship matrix for the individuals and estimates
# the phylogenetic bias related to each SNP.
getPhyloBias <- function(genotype, 
                         k.matrix) {
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
# RPA for continuous data
getRpaCont <- function(posterior, 
                       beta.mean, 
                       site, 
                       g1, 
                       g0, 
                       n1, 
                       n0, 
                       hdi.level, 
                       rpa.iterations, 
                       model.stan, 
                       mcmc.iterations, 
                       mcmc.warmup, 
                       mcmc.chains,
                       ...) {
  # ppc
  rpaRun <- function(p, n1, n0, model.stan, mcmc.iterations,
                     mcmc.warmup, mcmc.chains, hdi.level, ...) {
    
    y1 <- p[1]+p[2]*1+p[3]*stats::rt(n = n1, df = p[4])
    y0 <- p[1]+p[2]*0+p[3]*stats::rt(n = n0, df = p[4])
    data.list <- list(X = rep(x = c(1, 0), times = c(n1, n0)),
                      Y = c(y1, y0), Z = n0+n1)
    
    # get initial parameter values
    control <- list(adapt_delta = list(...)[["adapt_delta"]],
                    max_treedepth = list(...)[["max_treedepth"]])
    ppc.posterior <- rstan::sampling(object = model.stan,
                                     data = data.list,
                                     pars = c("alpha", "beta", "sigma", "nu"),
                                     iter = mcmc.iterations,
                                     warmup = mcmc.warmup,
                                     chains = mcmc.chains,
                                     cores = 1,
                                     control = control,
                                     verbose = FALSE,
                                     refresh = -1)
    
    # compute posterior summary
    hdi.L <- (1-hdi.level)/2
    hdi.H <- 1-(1-hdi.level)/2
    stats <- rstan::summary(object = ppc.posterior, 
                            pars = c("alpha", "beta", "sigma", "nu"), 
                            prob = c(hdi.L, hdi.H))$summary
    
    # return
    return(data.frame(beta.mean = stats["beta", "mean"],
                      beta.L = stats["beta", paste(hdi.L*100,"%",sep='')],
                      beta.H = stats["beta", paste(hdi.H*100,"%",sep='')],
                      stringsAsFactors = FALSE))
  }
  
  # get subset of posterior
  rpa.i <- sample(x = 1:nrow(posterior), size = rpa.iterations, replace = TRUE)
  posterior <- posterior[rpa.i, ]
  
  
  rpa <- apply(X = posterior[, c("alpha", "beta", "sigma", "nu")], MARGIN = 1, 
               FUN = rpaRun, n1 = n1, n0 = n0, hdi.level = hdi.level, 
               model.stan = model.stan, mcmc.iterations = mcmc.iterations, 
               mcmc.warmup = mcmc.warmup, mcmc.chains = mcmc.chains,
               adapt_delta = list(...)[["adapt_delta"]],
               max_treedepth = list(...)[["max_treedepth"]])
  rpa <- do.call(rbind, rpa)
  
  
  rpa$site <- site
  rpa$g1 <- g1
  rpa$g0 <- g0
  rpa$n1 <- n1
  rpa$n0 <- n0
  
  
  # compute stats
  # power error test
  power.error <- ifelse(test = (rpa$beta.L <= 0 & rpa$beta.H >= 0), 
                        yes = TRUE, no = FALSE)
  power.error <- sum(power.error)/length(power.error)
  
  
  # sign error test
  sign.error <- ifelse(test = (rpa$beta.mean <= 0 & beta.mean >= 0) 
                       | (rpa$beta.mean >= 0 & beta.mean <= 0), 
                       yes = TRUE, no = FALSE)
  sign.error <- sum(sign.error)/length(sign.error)
  
  
  # mean and sd of RPA-estimated effect
  rpa.beta.mean <- mean(rpa$beta.mean)
  rpa.beta.sd <- sd(rpa$beta.mean)
  
  # collect stats
  rpa.stats <- data.frame(site = site, 
                          g1 = g1, g0 = g0, 
                          n1 = n1, n0 = n0, 
                          rpa.power.error = power.error, 
                          rpa.sign.error = sign.error, 
                          rpa.beta.mean = rpa.beta.mean, 
                          rpa.beta.sd = rpa.beta.sd,
                          rpa.N = nrow(rpa),
                          stringsAsFactors = FALSE)
  
  return (rpa.stats)
}





# Description:
# RPA for dichotomous data
getRpaDich <- function(posterior, 
                       beta.mean, 
                       site, 
                       g1, 
                       g0, 
                       n1, 
                       n0, 
                       hdi.level, 
                       rpa.iterations, 
                       model.stan, 
                       mcmc.iterations, 
                       mcmc.warmup, 
                       mcmc.chains,
                       ...) {
  
  
  # ppc
  rpaRun <- function(p, n1, n0, model.stan, mcmc.iterations,
                     mcmc.warmup, mcmc.chains, hdi.level, ...) {
    y1 <- stats::rbinom(n = 1, prob = 1/(1+exp(-(p[1]+p[2]*1))), size = n1)
    y0 <- stats::rbinom(n = 1, prob = 1/(1+exp(-(p[1]+p[2]*0))), size = n0)
    data.list <- list(X = c(1, 0), Y = c(y1, y0), N = c(n1, n0), Z = 2)
    
    # get initial parameter values
    control <- list(adapt_delta = list(...)[["adapt_delta"]],
                    max_treedepth = list(...)[["max_treedepth"]])
    ppc.posterior <- rstan::sampling(object = model.stan,
                                     data = data.list,
                                     pars = c("alpha", "beta"),
                                     iter = mcmc.iterations,
                                     warmup = mcmc.warmup,
                                     chains = mcmc.chains,
                                     cores = 1,
                                     control = control,
                                     verbose = FALSE,
                                     refresh = -1)
    
    # compute posterior summary
    hdi.L <- (1-hdi.level)/2
    hdi.H <- 1-(1-hdi.level)/2
    stats <- rstan::summary(object = ppc.posterior, 
                            pars = c("alpha", "beta"), 
                            prob = c(hdi.L, hdi.H))$summary
    
    # return
    return(data.frame(beta.mean = stats["beta", "mean"],
                      beta.L = stats["beta", paste(hdi.L*100,"%",sep='')],
                      beta.H = stats["beta", paste(hdi.H*100,"%",sep='')],
                      stringsAsFactors = FALSE))
  }
  
  
  # get subset of posterior
  rpa.i <- sample(x = 1:nrow(posterior), size = rpa.iterations, replace = TRUE)
  posterior <- posterior[rpa.i, ]
  
  
  rpa <- apply(X = posterior[, c("alpha", "beta")], MARGIN = 1, 
               FUN = rpaRun, n1 = n1, n0 = n0, hdi.level = hdi.level, 
               model.stan = model.stan, mcmc.iterations = mcmc.iterations, 
               mcmc.warmup = mcmc.warmup, mcmc.chains = mcmc.chains,
               adapt_delta = list(...)[["adapt_delta"]],
               max_treedepth = list(...)[["max_treedepth"]])
  rpa <- do.call(rbind, rpa)
  rpa$site <- site
  rpa$g1 <- g1
  rpa$g0 <- g0
  rpa$n1 <- n1
  rpa$n0 <- n0
  
  
  # compute stats
  # power error test
  power.error <- ifelse(test = (rpa$beta.L <= 0 & rpa$beta.H >= 0), 
                        yes = TRUE, no = FALSE)
  power.error <- sum(power.error)/length(power.error)
  
  
  # sign error test
  sign.error <- ifelse(test = (rpa$beta.mean <= 0 & beta.mean >= 0) 
                       | (rpa$beta.mean >= 0 & beta.mean <= 0), 
                       yes = TRUE, no = FALSE)
  sign.error <- sum(sign.error)/length(sign.error)
  
  
  # mean and sd of RPA-estimated effect
  rpa.beta.mean <- mean(rpa$beta.mean)
  rpa.beta.sd <- sd(rpa$beta.mean)
  
  
  # collect stats
  rpa.stats <- data.frame(site = site, 
                          g1 = g1, g0 = g0, 
                          n1 = n1, n0 = n0, 
                          rpa.power.error = power.error, 
                          rpa.sign.error = sign.error, 
                          rpa.beta.mean = rpa.beta.mean, 
                          rpa.beta.sd = rpa.beta.sd,
                          rpa.N = nrow(rpa),
                          stringsAsFactors = FALSE)
  
  return (rpa.stats)
}





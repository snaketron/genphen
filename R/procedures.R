
# Description:
# Parse input data and format it for stan and statistical learning
getStanDataBackup <- function(genotype, 
                              phenotype, 
                              phenotype.type) {
  
  
  # Description:
  # Returns TRUE if the genotype data is composed of at most two-alleles 
  # at each SNP.
  isBiallelic <- function(genotype) {
    
    isBi <- function(x) {
      return(length(unique(x)) <= 2)
    }
    
    is.bi <- apply(X = genotype, MARGIN = 2, FUN = isBi)
    
    return(all(is.bi == TRUE))
  }
  
  
  # Description:
  # Get genotype-phenotype data
  getFormattedGenphen <- function(genotype, 
                                  phenotype, 
                                  phenotype.type) {
    
    
    # convert AAMultipleAlignment to matrix if needed
    genotype <- convertMsaToGenotype(genotype = genotype)
    
    
    # if vector genotype => matrix genotype
    if(is.vector(genotype)) {
      genotype <- matrix(data = genotype, ncol = 1)
    }
    
    
    # map the groups (genotypes) at each SNP to a number
    genotype.map <- matrix(data = 0, nrow = nrow(genotype),
                           ncol = ncol(genotype))
    for(i in 1:ncol(genotype)) {
      genotype.map[, i] <- as.numeric(as.factor(genotype[, i]))
    }
    
    
    
    # if vector genotype => matrix genotype
    if(is.vector(phenotype)) {
      phenotype <- matrix(data = phenotype, ncol = 1)
    }
    phenotype <- data.frame(phenotype)
    
    
    # TODO: check
    if(sum(phenotype.type == "D") != 0) {
      d <- which(phenotype.type == "D")
      for(i in 1:length(d)) {
        if(all(phenotype[, d] %in% c(1, 0)) == FALSE) {
          # mapping 1st element to 1, 2nd to 0
          u <- unique(phenotype[, d])
          phenotype[phenotype[, d] == u[1], d] <- "1"
          phenotype[phenotype[, d] == u[2], d] <- "0"
          cat("Mapping dichotomous phenotype:", i, 
              "(", u[1], "->1,", u[2], "->0) \n")
        }
        phenotype[, d] <- as.factor(as.character(phenotype[, d]))
      }
    }
    
    # return
    return (list(genotype = genotype, 
                 genotype.map = genotype.map,
                 phenotype = phenotype,
                 phenotype.type = phenotype.type))
  }
  
  
  # Description:
  # Convert genphen data to stan data
  getStanFormat <- function(f.data, 
                            is.bi) {
    
    # set flag
    f.data$is.bi <- TRUE
    
    if(is.bi == TRUE) {
      X <- f.data$genotype.map
      X[X == 2] <- -1
      f.data$Y <- f.data$phenotype
      f.data$X <- X
      
      f.data$Y <- f.data$phenotype
      f.data$Ns <- ncol(f.data$X)
      f.data$Ntq <- sum(f.data$phenotype.type == "Q")
      f.data$Ntd <- sum(f.data$phenotype.type == "D")
      f.data$N <- nrow(X)
      
      if(sum(phenotype.type == "Q") != 0) {
        f.data$Yq <- f.data$Y[, f.data$phenotype.type == "Q"]
      }
      if(sum(phenotype.type == "D") != 0) {
        f.data$Yd <- f.data$Y[, f.data$phenotype.type == "D"]
      }
      
      return (f.data)
    }
    else {
      X <- f.data$genotype.map
      if(ncol(X) > 1) {
        for(i in 2:ncol(X)) {
          X[, i] <- max(X[, (i-1)]) + X[, i]
        }
      }
      f.data$X <- X
      f.data$Y <- f.data$phenotype
      f.data$Ns <- ncol(f.data$X)
      f.data$Ntq <- sum(f.data$phenotype.type == "Q")
      f.data$Ntd <- sum(f.data$phenotype.type == "D")
      M <- data.frame(sk = as.vector(X), s = rep(x = 1:ncol(X), each = nrow(X)))
      M <- M[duplicated(M) == F, ]
      f.data$Ms <- M$s
      f.data$Msk <- M$sk
      f.data$Nsk <- max(M$sk)
      f.data$N <- nrow(X)
      
      if(sum(phenotype.type == "Q") != 0) {
        f.data$Yq <- as.matrix(f.data$Y[, f.data$phenotype.type == "Q"])
      }
      else {
        f.data$Yq <- matrix(data = 0, ncol = 0, nrow = f.data$N)
      }
      
      
      if(sum(phenotype.type == "D") != 0) {
        f.data$Yd <- as.matrix(f.data$Y[, f.data$phenotype.type == "D"])
      }
      else {
        f.data$Yd <- matrix(data = 0, ncol = 0, nrow = f.data$N)
      }
      return (f.data)
    }
  }
  
  
  
  # genphen data
  f.data <- getFormattedGenphen(genotype = genotype, 
                                phenotype = phenotype, 
                                phenotype.type = phenotype.type)
  
  
  # stan data
  stan.data <- getStanFormat(f.data = f.data, 
                             is.bi = isBiallelic(f.data$genotype.map))
  
  
  
  
  
  # return
  return (stan.data)
}





# Description:
# Parse input data and format it for stan and statistical learning
getStanData <- function(genotype, 
                        phenotype, 
                        phenotype.type) {
  
  
  
  # Description:
  # Returns TRUE if the genotype data is composed of at most two-alleles 
  # at each SNP.
  isBiallelic <- function(genotype) {
    
    isBi <- function(x) {
      return(length(unique(x)) <= 2)
    }
    
    is.bi <- apply(X = genotype, MARGIN = 2, FUN = isBi)
    
    return(all(is.bi == TRUE))
  }
  
  
  
  # Description:
  # Get genotype-phenotype data
  getFormattedGenphen <- function(genotype, 
                                  phenotype, 
                                  phenotype.type) {
    
    
    # convert AAMultipleAlignment to matrix if needed
    genotype <- convertMsaToGenotype(genotype = genotype)
    
    
    # if vector genotype => matrix genotype
    if(is.vector(genotype)) {
      genotype <- matrix(data = genotype, ncol = 1)
    }
    
    
    # map the groups (genotypes) at each SNP to a number
    genotype.map <- matrix(data = 0, nrow = nrow(genotype),
                           ncol = ncol(genotype))
    for(i in 1:ncol(genotype)) {
      genotype.map[, i] <- as.numeric(as.factor(genotype[, i]))
    }
    
    
    # if vector genotype => matrix genotype
    if(is.vector(phenotype)) {
      phenotype <- matrix(data = phenotype, ncol = 1)
    }
    phenotype <- data.frame(phenotype)
    
    
    # TODO: check
    if(sum(phenotype.type == "D") != 0) {
      d <- which(phenotype.type == "D")
      for(i in 1:length(d)) {
        if(all(phenotype[, d] %in% c(1, 0)) == FALSE) {
          # mapping 1st element to 1, 2nd to 0
          u <- unique(phenotype[, d])
          phenotype[phenotype[, d] == u[1], d] <- "1"
          phenotype[phenotype[, d] == u[2], d] <- "0"
          cat("Mapping dichotomous phenotype:", i, 
              "(", u[1], "->1,", u[2], "->0) \n")
        }
        phenotype[, d] <- as.factor(as.character(phenotype[, d]))
      }
    }
    
    # return
    return (list(genotype = genotype, 
                 genotype.map = genotype.map,
                 phenotype = phenotype,
                 phenotype.type = phenotype.type))
  }
  
  
  
  # Description:
  # Convert genphen data to stan data
  getStanFormat <- function(f.data, 
                            is.bi) {
    
    
    # Split SNP into pairs of genotypes
    getSplitX <- function(x) {
      ux <- unique(x)
      ns <- length(ux)
      nc <- choose(n = ns, k = 2)
      x.split <- matrix(data = 0, nrow = length(x), ncol = nc)
      
      if(ns == 1) {
        x.split[, 1] <- 1
        x.map <- data.frame(ref = ux[1], alt = NA, 
                            refN = sum(x == ux[1]), 
                            altN = 0,
                            stringsAsFactors = FALSE)
      }
      else {
        x.map <- c()
        counter <- 1
        for(i in 1:(ns-1)) {
          for(j in (i+1):ns) {
            x.split[x == ux[i], counter] <- 1
            x.split[x == ux[j], counter] <- -1
            counter <- counter + 1
            x.map <- rbind(x.map, data.frame(ref = ux[i], alt = ux[j], 
                                             refN = sum(x == ux[i]),
                                             altN = sum(x == ux[j]),
                                             stringsAsFactors = FALSE))
          }
        }
      }
      
      # return
      return(list(x.split = x.split, 
                  x.map = x.map))
    }
    
    
    # X.data
    getXData <- function(x) {
      return (x$x.split)
    }
    
    
    # X.map
    getXMap <- function(x) {
      return (x$x.map)
    }
    
    
    
    # make huge matrix with predictors
    X <- f.data$genotype
    colnames(X) <- 1:ncol(X)
    x.data <- apply(X = X, MARGIN = 2, FUN = getSplitX)
    for(i in 1:length(x.data)) {
      x.data[[i]]$x.map$site <- i
    }
    X <- do.call(cbind, lapply(X = x.data, FUN = getXData))
    x.map <- do.call(rbind, lapply(X = x.data, FUN = getXMap))
    
    if(sum(phenotype.type == "Q") != 0) {
      Yq <- as.matrix(f.data$phenotype[, f.data$phenotype.type == "Q"])
    }
    else {
      Yq <- matrix(data = 0, ncol = 0, nrow = nrow(X))
    }
    
    if(sum(phenotype.type == "D") != 0) {
      Yd <- as.matrix(f.data$phenotype[, f.data$phenotype.type == "D"])
    }
    else {
      Yd <- matrix(data = 0, ncol = 0, nrow = nrow(X))
    }
    
    s <- list(X = X,
              Y = f.data$phenotype,
              Yq = Yq,
              Yd = Yd,
              N = nrow(X),
              Ns = ncol(f.data$genotype),
              Nsk = ncol(X),
              Ntq = ncol(Yq),
              Ntd = ncol(Yd),
              is.bi = is.bi,
              xmap = x.map)
    
    return (s)
  }
  
  
  
  # genphen data
  f.data <- getFormattedGenphen(genotype = genotype, 
                                phenotype = phenotype, 
                                phenotype.type = phenotype.type)
  
  
  # stan data
  stan.data <- getStanFormat(f.data = f.data, 
                             is.bi = isBiallelic(f.data$genotype.map))
  
  
  
  # return
  return (stan.data)
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
                       mcmc.steps, 
                       mcmc.warmup, 
                       cores, 
                       hdi.level, 
                       stat.learn.method, 
                       cv.steps, 
                       rpa.iterations, 
                       with.stan.obj) {
  
  checkGenotypePhenotype <- function(genotype, phenotype, phenotype.type) {
    # CHECK: genotype
    if(is.null(attr(genotype, "class")) == FALSE) {
      if(!attr(genotype, "class") %in% c("AAMultipleAlignment",
                                         "DNAMultipleAlignment")) {
        stop("genotype can be one of the following structures: vector (for a 
             single SNP), matrix, data.frame or Biostrings structures such as
             DNAMultipleAlignment or AAMultipleAlignment")
      }
      else {
        temp <- as.matrix(genotype)
        if(nrow(temp) < 2 | ncol(temp) == 0) {
          stop("The genotypes cannot have less than two observations, or the
               number of genotypes cannot be 0.")
        }
      }
    }
    else {
      if(is.vector(genotype)) {
        genotype <- matrix(data = genotype, ncol = 1)
      }
      
      if(nrow(genotype) < 2 | ncol(genotype) == 0) {
        stop("The genotypes cannot have less than two observations or the
             number of genotypes cannot be 0.")
      }
      
      if(!is.matrix(genotype) & !is.data.frame(genotype)) {
        stop("genotype can be one of the following structures: vector (for a 
             single SNP), matrix, data.frame or Biostrings structures such as
             DNAMultipleAlignment or AAMultipleAlignment")
      }
      
      if(typeof(genotype) != "character") {
        stop("If it is structured as vector/matrix/data.frame, 
             the genotype have be of character type.")
      }
    }
    
    
    
    # CHECK: phenotype
    if(!is.vector(phenotype) & !is.matrix(phenotype)) {
      stop("The phenotype must be either a vector (single phenotype) or matrix 
           (with multiple phenotypes = columns), where the rows match the rows 
           of the genotype data")
    }
    
    # convert vector -> matrix 
    if(is.vector(phenotype) == TRUE) {
      phenotype <- matrix(data = phenotype, ncol = 1)
    }
    
    if(!is.numeric(phenotype)) {
      stop("The phenotype must be of numeric type.")
    }
    
    if(length(phenotype) < 3) {
      stop("The phenotype must contain at least 3 data points.")
    }
    
    if(nrow(genotype) != nrow(phenotype)) {
      stop("length(genotype) != length(phenotype),
           they must be equal in length.")
    }
    
    
    
    # CHECK: phenotype.type 
    if(is.vector(phenotype.type) == FALSE) {
      stop("phenotype.type must be vector. Each element in this vector refers 
           to the type of each phenotype, with 'Q' (for quantitative phenotypes) 
           or 'D' (for dichotomous) included as each columns in the phenotype 
           data. If a single phenotype is provided, then the phenotype type 
           should be a single 'Q' or 'D'.")
    }
    
    if(length(phenotype.type) == 0) {
      stop("The phenotype.type vector must contain at least 1 element.")
    }
    
    if(typeof(phenotype.type) != "character") {
      stop("phenotype.type must be character vector with elements 'Q' (for 
           quantitative phenotypes) or 'D' (for dichotomous)")
    }
    
    if(ncol(phenotype) != length(phenotype.type)) {
      stop("Number of phenotypes provided differs from phenotypes types.")
    }
    
    if(all(phenotype.type %in% c("Q", "D")) == FALSE) {
      stop("phenotype.type must be character vector with elements 'Q' (for 
           quantitative phenotypes) or 'D' (for dichotomous)")
    }
  }
  
  checkPhenotypeValidity <- function(phenotype, phenotype.type) {
    
    # convert vector -> matrix 
    if(is.vector(phenotype) == TRUE) {
      phenotype <- matrix(data = phenotype, ncol = 1)
    }
    
    # check phenotype types
    for(i in 1:length(phenotype.type)) {
      if(phenotype.type[i] == "D") {
        if(length(unique(phenotype[, i])) != 2) {
          stop("The dichotomous phenotypes must contains exactly two 
               categories (classes) \n")
        }
      }
      if(phenotype.type[i] == "Q") {
        if(length(unique(phenotype[, i])) <= 3) {
          warning("The quantitative phenotype/s contains 3 or less unique 
                  elements \n")
        }
      }
    }
  }
  
  checkModelType <- function(model.type) {
    # CHECK: model.type
    if(length(model.type) != 1) {
      stop("model.type must be a string (default = 'univariate')")
    }
    
    if(!is.character(model.type)) {
      stop("model.type must be a string: 'univariate' or 'hierarchical'")
    }
    
    if(!model.type %in% c("univariate", "hierarchical")) {
      stop("phenotype.type must be a string: 'univariate' or 'hierarchical'")
    }
  }
  
  checkMcmcIterations <- function(mcmc.steps) {
    # CHECK: mcmc.steps
    if(length(mcmc.steps) != 1) {
      stop("the mcmc.steps must be a number > 0 (default = 10000).")
    }
    
    if(!is.numeric(mcmc.steps)) {
      stop("mcmc.steps must be a numeric argument (default = 10000).")
    }
    
    if(mcmc.steps <= 0) {
      stop("mcmc.steps must be larger than 0 (default = 10000).")
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
      stop("stat.learn.method must be a string: 'rf' or 'svm'")
    }
    
    if(!is.character(stat.learn.method)) {
      stop("stat.learn.method must be a string: 'rf', 'svm' or 'none'")
    }
    
    if(!stat.learn.method %in% c("rf", "svm", "none")) {
      stop("stat.learn.method must be a string: 'rf', 'svm' or 'none'")
    }
  }
  
  checkCv <- function(stat.learn.method, cv.steps) {
    if(stat.learn.method %in% c("rf", "svm")) {
      if(length(cv.steps) != 1) {
        stop("cv.steps must be a number (default = 1,000).")
      }
      
      if(is.numeric(cv.steps) == FALSE) {
        stop("cv.steps must be a number (default = 1,000).")
      }
      
      if(cv.steps < 100) {
        stop("cv.steps >= 100 recomended (default = 1,000).")
      }
    }
  }
  
  if(is.null(genotype) | missing(genotype) |
     is.null(phenotype) | missing(phenotype) |
     is.null(phenotype.type) | missing(phenotype.type) |
     is.null(model.type) | missing(model.type) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.steps) | missing(mcmc.steps) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(cores) | missing(cores) |
     is.null(hdi.level) | missing(hdi.level) |
     is.null(stat.learn.method) | missing(stat.learn.method) |
     is.null(cv.steps) | missing(cv.steps)) {
    stop("arguments must be non-NULL/specified")
  }
  
  checkGenotypePhenotype(genotype = genotype, 
                         phenotype = phenotype, 
                         phenotype.type = phenotype.type)
  checkPhenotypeValidity(phenotype = phenotype, 
                         phenotype.type = phenotype.type)
  checkModelType(model.type = model.type)
  checkMcmcIterations(mcmc.steps = mcmc.steps)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkCores(cores = cores)
  checkHdi(hdi.level = hdi.level)
  checkMlMethod(stat.learn.method = stat.learn.method)
  checkCv(stat.learn.method = stat.learn.method, cv.steps = cv.steps)
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
runStatLearn <- function(genphen.data,
                         method,
                         cv.fold, 
                         cv.steps,
                         ntree,
                         hdi.level,
                         cores,
                         dot.param) {
  
  
  # RF analysis
  runRf <- function(X, Y, cv.fold, cv.steps, 
                    ntree, hdi.level, site) {
    
    
    getBoot <- function(D, cv.fold, cv.steps, ntree, hdi.level) {
      
      # posterior output (one extra for multi-trait)
      if(ncol(D) > 2) {
        posterior <- vector(mode = "list", length = ncol(D))
      }
      else {
        posterior <- vector(mode = "list", length = 1)
      }
      for(i in 1:length(posterior)) {
        posterior[[i]] <- matrix(data = NA, nrow = cv.steps, ncol = 2)
      }
      rm(i)
      
      
      D$Y <- as.character(D$Y)
      for(i in 1:cv.steps) {
        # sample at random
        s <- sample(x = 1:nrow(D), size = ceiling(x = cv.fold*nrow(D)), 
                    replace = FALSE)
        train <- D[s, ]
        test <- D[-s, ]
        
        
        # TODO: dummy (check if inference needed at all e.g. 1 class only)
        if(length(unique(train$Y)) == 1) {
          for(j in 1:length(posterior)) {
            posterior[[j]][i, ] <- c(NA, NA)
          }
        }
        else {
          for(j in 1:length(posterior)) {
            # train classification model (try condition to avoid errors in case
            # only one-level predictor is train data)
            if(j == length(posterior) & j != 1) {
              rf.out <- try(ranger::ranger(Y~., data = train,
                                           num.trees = ntree), 
                            silent = TRUE)
            }
            else {
              rf.out <- try(ranger::ranger(Y~., data = train[, c(1, j+1)],
                                           num.trees = ntree), 
                            silent = TRUE)
            }
            
            
            if(class(rf.out) == "try-error") {
              posterior[[j]][i, 1:2] <- c(NA, NA)
            }
            else {
              # test classification model
              pr <- stats::predict(object = rf.out, data = test)
              
              
              # compute classification accuracy (1 - classification error)
              test$Y <- as.character(test$Y)
              pr$predictions <- as.character(pr$predictions)
              ca <- sum(test$Y == pr$predictions)/nrow(test)
              
              
              # compute k statistics
              k <- getKappa(real = test$Y, 
                            predicted = pr$predictions, 
                            aas = unique(D$Y))
              
              posterior[[j]][i, 1:2] <- c(ca, k)
            }
          }
        }
      }
      
      # compute stats
      summary <- c()
      for(j in 1:length(posterior)) {
        ca <- posterior[[j]][, 1]
        ca <- ca[is.finite(ca)]
        ca.mean <- mean(x = ca)
        # get HDI
        ca.hdi <- getHdi(vec = posterior[[j]][,1], hdi.level = hdi.level)
        ca.L <- as.numeric(ca.hdi[1])
        ca.H <- as.numeric(ca.hdi[2])
        
        k <- posterior[[j]][, 2]
        k <- k[is.finite(k)]
        k.mean <- mean(x = k)
        # get HDI
        k.hdi <- getHdi(vec = posterior[[j]][, 2], hdi.level = hdi.level)
        k.L <- as.numeric(k.hdi[1])
        k.H <- as.numeric(k.hdi[2])
        
        # summary append
        row <- data.frame(ca = ca.mean, ca.L = ca.L, ca.H = ca.H,
                          k = k.mean, k.L = k.L, k.H = k.H, i = j)
        summary <- rbind(summary, row)
      }
      
      return (summary)
    }
    
    
    gs <- sort(unique(X), decreasing = T)
    result <- c()
    if(length(gs) == 1) {
      D <- data.frame(Y = X)
      D <- cbind(D, Y)
      
      for(i in 1:ncol(Y)) {
        out <- data.frame(ca = NA, ca.L = NA, ca.H = NA,
                          k = NA, k.L = NA, k.H = NA,
                          i = i, ref = gs, alt = NA, 
                          refN = nrow(D), altN = NA, 
                          site = site, stringsAsFactors = FALSE)
        result <- rbind(result, out)
      }
    }
    else {
      for(gi in 1:(length(gs)-1)) {
        for(gj in (gi+1):length(gs)) {
          
          i <- which(X %in% c(gs[gi], gs[gj]))
          D <- data.frame(Y = X)
          colnames(Y) <- paste("P", 1:ncol(Y), sep = '')
          D <- cbind(D, Y)
          D <- D[i, ]
          
          
          # if c.v. can be done with given cv.fold and # of data points
          if(ceiling(cv.fold*nrow(D)) == nrow(D)) {
            for(i in 1:ncol(Y)) {
              out <- data.frame(ca = NA, ca.L = NA, ca.H = NA,
                                k = NA, k.L = NA, k.H = NA,
                                i = i, ref = gs[gi], alt = gs[gj], 
                                refN = sum(X == gs[gi]), 
                                altN = sum(X == gs[gj]), 
                                site = site, stringsAsFactors = FALSE)
              result <- rbind(result, out)
            }
          }
          else {
            # run
            p <- getBoot(D = D, 
                         cv.fold = cv.fold, 
                         cv.steps = cv.steps, 
                         hdi.level = hdi.level,
                         ntree = ntree)
            
            # summarize
            p$ref <- gs[gi]
            p$alt <- gs[gj]
            p$refN = sum(X == gs[gi])
            p$altN = sum(X == gs[gj])
            p$site = site
            result <- rbind(result, p)
          }
        }
      }
    }
    return (result)
  }
  
  
  # SVM analysis
  runSvm <- function(X, Y, cv.fold, cv.steps, 
                     hdi.level, site) {
    
    
    getBoot <- function(D, cv.fold, cv.steps, hdi.level) {
      
      # posterior output (one extra for multi-trait)
      if(ncol(D) > 2) {
        posterior <- vector(mode = "list", length = ncol(D))
      }
      else {
        posterior <- vector(mode = "list", length = 1)
      }
      for(i in 1:length(posterior)) {
        posterior[[i]] <- matrix(data = NA, nrow = cv.steps, ncol = 2)
      }
      rm(i)
      
      
      D$Y <- as.character(D$Y)
      for(i in 1:cv.steps) {
        # sample at random
        s <- sample(x = 1:nrow(D), size = ceiling(x = cv.fold*nrow(D)), 
                    replace = FALSE)
        train <- D[s, ]
        test <- D[-s, ]
        
        
        # TODO: dummy (check if inference needed at all e.g. 1 class only)
        if(length(unique(train$Y)) == 1) {
          for(j in 1:length(posterior)) {
            posterior[[j]][i, ] <- c(NA, NA)
          }
        }
        else {
          for(j in 1:length(posterior)) {
            # train classification model (try condition to avoid errors in case
            # only one-level predictor is train data)
            if(j == length(posterior) & j != 1) {
              svm.out <- try(e1071::svm(as.factor(Y) ~ ., 
                                        data = train,
                                        type = "C-classification"), 
                             silent = TRUE)
            }
            else {
              svm.out <- try(e1071::svm(as.factor(Y) ~ ., 
                                        data = train[, c(1, j+1)],
                                        type = "C-classification"), 
                             silent = TRUE)
            }
            
            
            if(class(svm.out)[1] == "try-error") {
              posterior[[j]][i, 1:2] <- c(NA, NA)
            }
            else {
              # test classification model
              pr <- stats::predict(object = svm.out, newdata = test)
              
              
              # compute classification accuracy (1 - classification error)
              test$Y <- as.character(test$Y)
              pr <- as.character(pr)
              ca <- sum(test$Y == pr)/nrow(test)
              
              
              # compute k statistics
              k <- getKappa(real = test$Y, predicted = pr, aas = unique(D$Y))
              
              posterior[[j]][i, 1:2] <- c(ca, k)
            }
          }
        }
      }
      
      # compute stats
      summary <- c()
      for(j in 1:length(posterior)) {
        ca <- posterior[[j]][, 1]
        ca <- ca[is.finite(ca)]
        ca.mean <- mean(x = ca)
        # get HDI
        ca.hdi <- getHdi(vec = posterior[[j]][,1], hdi.level = hdi.level)
        ca.L <- as.numeric(ca.hdi[1])
        ca.H <- as.numeric(ca.hdi[2])
        
        k <- posterior[[j]][, 2]
        k <- k[is.finite(k)]
        k.mean <- mean(x = k)
        # get HDI
        k.hdi <- getHdi(vec = posterior[[j]][, 2], hdi.level = hdi.level)
        k.L <- as.numeric(k.hdi[1])
        k.H <- as.numeric(k.hdi[2])
        
        
        # summary append
        row <- data.frame(ca = ca.mean, ca.L = ca.L, ca.H = ca.H,
                          k = k.mean, k.L = k.L, k.H = k.H, i = j)
        summary <- rbind(summary, row)
      }
      
      return (summary)
    }
    
    
    gs <- sort(unique(X), decreasing = T)
    result <- c()
    if(length(gs) == 1) {
      D <- data.frame(Y = X)
      D <- cbind(D, Y)
      
      for(i in 1:ncol(Y)) {
        out <- data.frame(ca = NA, ca.L = NA, ca.H = NA,
                          k = NA, k.L = NA, k.H = NA,
                          i = i, ref = gs, alt = NA, 
                          refN = nrow(D), altN = NA,
                          site = site, stringsAsFactors = FALSE)
        result <- rbind(result, out)
      }
    }
    else {
      for(gi in 1:(length(gs)-1)) {
        for(gj in (gi+1):length(gs)) {
          
          
          i <- which(X %in% c(gs[gi], gs[gj]))
          D <- data.frame(Y = X)
          colnames(Y) <- paste("P", 1:ncol(Y), sep = '')
          D <- cbind(D, Y)
          D <- D[i, ]
          
          # if c.v. can be done with given cv.fold and # of data points
          if(ceiling(cv.fold*nrow(D)) == nrow(D)) {
            for(i in 1:ncol(Y)) {
              out <- data.frame(ca = NA, ca.L = NA, ca.H = NA,
                                k = NA, k.L = NA, k.H = NA,
                                i = i, ref = gs[gi], alt = gs[gj],
                                refN = sum(X == gs[gi]),
                                altN = sum(X == gs[gj]),
                                site = site, stringsAsFactors = FALSE)
              result <- rbind(result, out)
            }
          }
          else {
            # run
            p <- getBoot(D = D, 
                         cv.fold = cv.fold, 
                         cv.steps = cv.steps, 
                         hdi.level = hdi.level)
            
            # summarize
            p$ref <- gs[gi]
            p$alt <- gs[gj]
            p$refN = sum(X == gs[gi])
            p$altN = sum(X == gs[gj])
            p$site = site
            result <- rbind(result, p)
          }
        }
      }
    }
    return (result)
  }
  
  
  # multicore classification
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  if(method == "rf") {
    cas <- (foreach(j = 1:ncol(genphen.data$genotype),
                    .export = c("getHdi", "getKappa"),
                    .packages = c("ranger")) %dopar%
              runRf(X = as.matrix(genphen.data$genotype[, j]),
                    Y = genphen.data$Y,
                    cv.fold = dot.param[["cv.fold"]],
                    cv.steps = cv.steps,
                    hdi.level = hdi.level,
                    ntree = dot.param[["ntree"]],
                    site = j))
  }
  else if(method == "svm") {
    cas <- (foreach(j = 1:ncol(genphen.data$genotype),
                    .export = c("getHdi", "getKappa"),
                    .packages = c("e1071")) %dopar%
              runSvm(X = as.matrix(genphen.data$genotype[, j]),
                     Y = genphen.data$Y,
                     cv.fold = dot.param[["cv.fold"]],
                     cv.steps = cv.steps,
                     hdi.level = hdi.level,
                     site = j))
  }
  # stop cluster
  parallel::stopCluster(cl = cl)
  doParallel::stopImplicitCluster()
  
  return (cas)
}




# Description:
# Bayesian inference
runBayesianInference <- function(genphen.data,
                                 mcmc.chains,
                                 mcmc.steps,
                                 mcmc.warmup,
                                 cores,
                                 model.stan,
                                 ...) {
  
  # extra conversion
  Yd <- c()
  if(genphen.data$Ntd > 0) {
    for(i in 1:genphen.data$Ntd) {
      Yd <- cbind(Yd, as.numeric(genphen.data$Yd[, i]))
    }
  }
  genphen.data$Yd <- Yd
  
  
  data.list <- list(N = genphen.data$N, 
                    Ntq = genphen.data$Ntq, 
                    Ntd = genphen.data$Ntd, 
                    Ns = genphen.data$Ns, 
                    Nsk = genphen.data$Nsk, 
                    Yq = genphen.data$Yq, 
                    Yd = genphen.data$Yd, 
                    X = genphen.data$X, 
                    M = genphen.data$Ms)
  
  
  # get initial parameter values
  control <- list(adapt_delta = list(...)[["adapt_delta"]],
                  max_treedepth = list(...)[["max_treedepth"]])
  refresh <- list(...)[["refresh"]]
  verbose <- list(...)[["verbose"]]
  posterior <- rstan::sampling(object = model.stan,
                               data = data.list,
                               iter = mcmc.steps,
                               warmup = mcmc.warmup,
                               chains = mcmc.chains,
                               cores = cores,
                               control = control,
                               verbose = verbose,
                               refresh = refresh)
  
  # return
  return (list(posterior = posterior))
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












# Description:
# It plots the results of the different method runGenphen* methods. Each result
# entry is plotted as point with respect to its effect size and classification
# accuracy attributes.
# Return: ggplot object of the plot
plotGenphenRfSvm <- function(genphen.results = NULL) {

  # check for NA or NULL
  if(is.null(genphen.results)) {
    stop("the 'genphen.results' is NULL.")
  }

  # check for empty row
  if(all(is.na(genphen.results))) {
    stop("the 'genphen.results' contains only NAs")
  }

  # check for 0 row count
  if(nrow(genphen.results) == 0) {
    stop("the 'genphen.results' has 0 rows.")
  }

  # tricking R-CMD check
  d <- ca <- ca.hdi.high <- ca.hdi.low <- NULL

  g <- ggplot2::ggplot(data = genphen.results)+
    geom_errorbarh(aes(xmax = ca.hdi.H, xmin = ca.hdi.L,
                       x = ca, y = abs(d)), height = 0.025)+
    geom_point(aes(x = ca, y = abs(d), fill = kappa),
               col = "black", pch = 21, size = 3)+
    theme_bw()+
    xlab('Classification accuracy')+
    ylab('Effect size')+
    theme(text = element_text(size = 12))+
    xlim(c(0, 1))+
    scale_fill_gradientn(colours = terrain.colors(10))
  print(g)
  return(g)
}





# Description:
# It plots the results of the different method runGenphen* methods. Each result
# entry is plotted as point with respect to its effect size and classification
# accuracy attributes.
# Return: ggplot object of the plot
plotGenphenBayes <- function(genphen.results = NULL, hdi) {

  # check for NA or NULL
  if(is.null(genphen.results)) {
    stop("the 'genphen.results' is NULL.")
  }

  # check for empty row
  if(all(is.na(genphen.results))) {
    stop("the 'genphen.results' contains only NAs")
  }

  # check for 0 row count
  if(nrow(genphen.results) == 0) {
    stop("the 'genphen.results' has 0 rows.")
  }

  l <- which(regexpr(pattern = paste("mu.effect", hdi, "L", sep = '.'),
                     text = colnames(genphen.results)) != -1)
  if(length(l) != 1) {
    stop("select correct L hdi.")
  }

  h <- which(regexpr(pattern = paste("mu.effect", hdi, "H", sep = '.'),
                     text = colnames(genphen.results)) != -1)
  if(length(h) != 1) {
    stop("select correct H hdi.")
  }



  # tricking R-CMD check
  mu.effect <- NULL

  g <- ggplot2::ggplot(data = genphen.results)+
    geom_errorbarh(aes(xmax = genphen.results[, h],
                       xmin = genphen.results[, l],
                       x = mu.effect,
                       y = as.factor(paste(g.1, site, g.2, sep = ''))),
                   height = 0.025)+
    geom_point(aes(x = mu.effect,
                   y = as.factor(paste(g.1, site, g.2, sep = ''))),
               fill = "white", col = "black", pch = 21, size = 3)+
    theme_bw()+
    xlab('Bayesian effect size (mu.effect)')+
    ylab('Polymorphism')+
    theme(legend.position = "none")+
    theme(text = element_text(size = 12))+
    geom_vline(xintercept = 0, col = "red", linetype = "dashed")
  print(g)
  return(g)
}






# Description:
# Given a genotype data, a vector of the corresponding phenotypes and a
# genotype index, this method visualizes the phenotypic distribution as a
# function of the different genotype states of the indexed specific genotype.
# Return: ggplot object
plotSpecificGenotype <- function(genotype, phenotype, index) {


  # Description:
  # Provided the input arguments, this function checks their validity. It
  # stops the execution if a problem is encountered and prints out warnings.
  checkInputs <- function(genotype, phenotype, index = NULL) {

    if(is.null(genotype) | is.null(phenotype) | is.null(index)) {
      stop("the arguments genotype, phenotype and index must be specified.")
    }

    # is it matrix/dataframe or DNAMultipleAlignment
    if(is.null(attr(genotype, "class")) == FALSE) {
      if(!attr(genotype, "class") %in% c("DNAMultipleAlignment",
                                         "AAMultipleAlignment")) {
        stop("genotype can be one of the following structures: matrix or
             data.frame or DNAMultipleAlignment or AAMultipleAlignment")
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
             data.frame or DNAMultipleAlignment or AAMultipleAlignment")
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


    if(!is.vector(index)) {
      stop("the index must be a vector.")
    }

    if(length(index) > 1) {
      stop("index is one-element vector")
    }

    if(!is.numeric(index)) {
      stop("the index must be of numeric type.")
    }

    if(index > ncol(as.matrix(genotype))) {
      stop("index must be in range of [1, ncol(genotype)]")
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



  # check inputs
  checkInputs(genotype = genotype, phenotype =  phenotype, index = index)



  # convert DNAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)



  data <- data.frame(snp = as.factor(as.character(genotype[, index])),
                     phenotype = as.numeric(as.character(phenotype)))


  col.values <- c("#cc0000", "#009933", "#5a00e1", "#be6e00", "#ccff00",
                  "#007f52", "#28e7b2", "#aebf23", "#d03dc2", "#c1e1da",
                  "#738782", "#000000", "#E1C1C8", "#1e0099", "#e230dd")


  if(length(unique(data$snp)) > 15) {
    stop("too many factors in the genotype data (max = 15), change the
         col.values in the code if more values are needed.")
  }


  # tricking R-CMD check
  snp <- phenotype <- NULL

  g <- ggplot()+
    stat_boxplot(data = data, aes(x = snp, y = phenotype, fill = snp),
                 geom ='errorbar')+
    geom_boxplot(data = data, aes(x = snp, y = phenotype, fill = snp))+
    theme_bw()+
    ylab("Phenotype")+
    xlab("Genotype state")+
    theme(plot.title = element_text(face = "bold"))+
    scale_fill_manual(name = "state", values = col.values)+
    theme(legend.position = "none")+
    theme(text = element_text(size = 12))

  print(g)
  return(g)
}




# Description:
# It plots the results of the different method runGenphen* methods using the
# p-values obtained from the ANOVA tests. The p-values are first corrected
# with FDR, then transformed to -log10, and finaly plotted as the so-called
# Manhattan plots.
# Return: ggplot object of the plot
plotManhattan <- function(genphen.results = NULL) {

  # check for NA or NULL
  if(is.null(genphen.results)) {
    stop("the 'genphen.results' is NULL.")
  }

  # check for empty row
  if(all(is.na(genphen.results))) {
    stop("the 'genphen.results' contains only NAs")
  }

  # check for 0 row count
  if(nrow(genphen.results) == 0) {
    stop("the 'genphen.results' has 0 rows.")
  }

  g <- ggplot(data = genphen.results)+
    geom_point(aes(x = 1:nrow(genphen.results), y = abs(t.test.pvalue)),
               fill = "white", col = "black", pch = 21, size = 3)+
    theme_bw()+
    xlab('genotypes')+
    ylab('-log10 (P-values)')+
    theme(legend.position = "right")+
    theme(text = element_text(size = 12))
  print(g)
  return(g)
}




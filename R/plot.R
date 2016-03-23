


# Description:
# It plots the results of the different method runGenphen* methods. Each result
# entry is plotted as point with respect to its effect size and classification
# accuracy attributes.
# Return: ggplot object of the plot
plotGenphenResults <- function(genphen.results) {

  # check for empty row
  if(all(is.na(genphen.results))) {
    stop("the 'genphen.results' contains only NAs")
  }

  # tricking R-CMD check
  classification.accuracy <- effect.size <- ca <- NULL

  g <- ggplot2::ggplot()+
    geom_point(data = genphen.results,
               aes(x = ca, y = abs(effect.size)), size = 2.5, col = "black")+
    geom_point(data = genphen.results,
               aes(x = ca, y = abs(effect.size), col = ca), size = 2)+
    theme_bw()+
    scale_color_gradient(low = "#c7d9e8",  high = "#2a4e6c")+
    xlab('Classification accuracy')+
    ylab('Effect size')+
    theme(legend.position = "none")+
    theme(text = element_text(size = 12))+
    xlim(c(0, 1))
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
plotManhattan <- function(genphen.results) {

  # check for empty row
  if(all(is.na(genphen.results))) {
    stop("the 'genphen.results' contains only NAs")
  }

  # correct with FDR and convert to -log
  genphen.results$corrected.anova <-
    stats::p.adjust(p = genphen.results$anova.score, method = "fdr")
  genphen.results$corrected.anova <-
    -log(x = genphen.results$corrected.anova, base = 10)

  # tricking R-CMD check
  corrected.anova <- NULL

  g <- ggplot(data = genphen.results)+
    geom_point(aes(x = 1:nrow(genphen.results), y = corrected.anova),
               col = "black", size = 2.4)+
    geom_point(aes(x = 1:nrow(genphen.results), y = corrected.anova,
                   col = corrected.anova), size = 2)+
    theme_bw()+
    scale_color_gradient(low = "#c7d9e8",  high = "#2a4e6c")+
    xlab('genotypes')+
    ylab('-log10 (p-value)')+
    theme(legend.position = "right")+
    theme(text = element_text(size = 12))
  print(g)
  return(g)
}



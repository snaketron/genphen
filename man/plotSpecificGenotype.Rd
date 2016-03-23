\name{plotSpecificGenotype}



\alias{plotSpecificGenotype}



\title{Visualizing specific genotypes}



\usage{
plotSpecificGenotype(genotype, phenotype, index)
}



\description{
This procedure visualizes the phenotypic distribution linked to each of the
genetic states of a specific genotype.
}




\arguments{
  \item{genotype}{Character matrix or data frame, containing SNPs/SAAPs as
  columns or alternatively as a DNAMultipleAlignment/AAMultipleAlignment
  Biostrings object.}
  \item{phenotype}{Numerical vector whose elements correspond to the
  genotype.}
  \item{index}{Index (number) of the specific genotype column within the
  genotype data which is to be plotted.}
}




\examples{
#Example 1:
data(genotype.snp) #or data(genotype.snp.msa)
data(phenotype.snp)
specific.genotype.plot <- plotSpecificGenotype(genotype = genotype.snp,
phenotype = phenotype.snp, index = 1)

#Example 2:
data(genotype.saap) #or ata(genotype.saap.msa)
data(phenotype.saap)
specific.genotype.plot <- plotSpecificGenotype(genotype = genotype.saap,
phenotype = phenotype.saap, index = 3)
}



\seealso{
runGenphenSaap, runGenphenSnp, plotGenphenResults, plotManhattan
}



\details{
This procedure allows the user to inspect a specific genotype with respect to
the the phenotype. It uses a boxplot notation to plot the phenotypes as a
function of the states of that genotype. The resulting boxplot will visualize
whether the different states of the specific genotype are linked to different
and disjoint phenotypic distributions, which is a signature of a strong
association between the genotype and the phenotype.
}



\value{
\item{plot}{ggplot object}
}



\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


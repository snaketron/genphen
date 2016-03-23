\name{plotGenphenResults}



\alias{plotGenphenResults}



\title{Visualizing GWAS results}



\usage{
plotGenphenResults(genphen.results)
}



\arguments{
  \item{genphen.results}{Data.frame resulting from runGenphenSnp or
  runGenphenSaap.}
}



\description{
This procedure plots the results obtained using runGenphenSnp or runGenphenSaap.
}



\examples{
#Example 1:
data(genotype.snp)
#or data(genotype.snp.msa) in this case you cannot subset genotype.snp[, 1:5]
data(phenotype.snp)
genphen.results <- runGenphenSnp(genotype = genotype.snp[, 1:5],
phenotype = phenotype.snp, technique = "svm", fold.cv = 0.66, boots = 100)
genphen.plot <- plotGenphenResults(genphen.results = genphen.results)

#Example 2:
data(genotype.saap)
#or data(genotype.saap.msa) in this case you cannot subset genotype.saap[, 1:5]
data(phenotype.saap)
genphen.results <- runGenphenSaap(genotype = genotype.saap[, 1:5],
phenotype = phenotype.saap, technique = "svm", fold.cv = 0.66, boots = 100)
genphen.plot <- plotGenphenResults(genphen.results = genphen.results)
}



\seealso{
runGenphenSaap, runGenphenSnp, plotGenphenResults, plotSpecificGenotype
}



\details{
This procedure plots the results obtained using runGenphenSnp or runGenphenSaap.
Each result entry is plotted as a point with respect to its effect size and
classification accuracy attributes, whereby the color of the points is directly
proportional to the classification accuracy. The region in the top-right corner
of the plot is where the genotypes which have the strongest association with
the phenotype are found.
}



\value{
\item{plot}{ggplot plot object.}
}



\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



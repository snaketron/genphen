\name{plotGenphenRfSvm}



\alias{plotGenphenRfSvm}



\title{Visualizing the genphen results of runGenphenRf or runGenphenSvm}



\usage{
plotGenphenRfSvm(genphen.results)
}



\arguments{
  \item{genphen.results}{data.frame resulting from runGenphenRf/runGenphenSvm.}
}



\description{
This procedure visualizes the results obtained using the functions runGenphenRf
or runGenphenSvm.
}



\examples{
#Example 1 (RF results):
data(genotype.snp)
#or data(genotype.snp.msa) in this case you cannot subset genotype.snp[, 1:5]
data(phenotype.snp)
genphen.rf <- runGenphenRf(genotype = genotype.snp[, 1:5],
                           phenotype = phenotype.snp,
                           cv.fold = 0.66,
                           cv.steps = 100,
                           hdi.level = 0.99,
                           ntree = 1000)
plotGenphenRfSvm(genphen.results = genphen.rf)
}



\seealso{
runGenphenRf, runGenphenSvm, plotGenphenBayes,
plotSpecificGenotype, plotManhattan
}



\details{
This procedure plots the results of the functions runGenphenRf or runGenphenSvm.
Each association is shown as a point with respect to its effect size (Cohen's
d) and the classification accuracy (CA). The colors of the points correspond to
the Cohen's kappas. The confidence intervals estimated for the CA point
estimates are provided here as well.
}



\value{
\item{plot}{ggplot plot object.}
}



\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



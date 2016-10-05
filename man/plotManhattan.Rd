\name{plotManhattan}



\alias{plotManhattan}



\title{Visualizing genphen results with Manhattan plots}



\usage{
plotManhattan(genphen.results)
}



\arguments{
  \item{genphen.results}{Data.frame resulting from runGenphenRf, runGenphenSvm
  or runGenphenBayesian.}
}



\description{
This procedure plots the complementary T-test results obtained using
runGenphenRf, runGenphenSvm or runGenphenBayes.
}



\examples{
#Example 1:
data(genotype.snp)
#or data(genotype.snp.msa) in this case you cannot subset genotype.snp[, 1:5]
data(phenotype.snp)
genphen.rf <- runGenphenRf(genotype = genotype.snp[, 1:5],
                           phenotype = phenotype.snp,
                           cv.fold = 0.66,
                           cv.steps = 100,
                           hdi.level = 0.99,
                           ntree = 1000)
plotManhattan(genphen.results = genphen.rf)
}



\seealso{
runGenphenRf, runGenphenSvm, plotGenphenBayes, plotGenphenRfSvm,
plotSpecificGenotype
}



\details{
This procedure plots the results from a two-sample T-test analysis, conducted in
each of the three functions: runGenphenRf, runGenphenSvm or runGenphenBayesian.
The P-values resulting from the T-test are first corrected by the Benjamini
-Hochberg correction procedure (FDR), followed by a -log10 transformation. These
values are then plotted as points, where each point corresponds to a single
genotype-phenotype association as seen in a classical Manhattan plot.
}



\value{
\item{plot}{ggplot plot object.}
}



\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



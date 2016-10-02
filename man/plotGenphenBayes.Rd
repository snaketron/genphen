\name{plotGenphenBayes}



\alias{plotGenphenBayes}



\title{Visualizing the genphen results of runGenphenBayes}



\usage{
plotGenphenBayes(genphen.results, hdi)
}



\arguments{
  \item{genphen.results}{data.frame resulting from runGenphenBayes.}
  \item{hdi}{single HDI level (e.g. 0.95) which has previously been estimated
  with runGenphenBayes}
}



\description{
This procedure visualizes the results obtained using the function
runGenphenBayes.
}



\examples{
#Example:
data(genotype.saap)
# or data(genotype.saap.msa) in this case you cannot subset genotype.saap[, 1:5]
data(phenotype.saap)

# run analysis
genphen.bayes <- runGenphenBayes(genotype = genotype.saap[, 1:5],
                                 phenotype = phenotype.saap,
                                 mcmc.iter = 1000,
                                 chain.nr = 4,
                                 model = "tdist",
                                 hdi.levels = c(0.95, 0.99))
plotGenphenBayes(genphen.results = genphen.bayes)
}



\seealso{
runGenphenRf, runGenphenSvm, runGenphenBayes, plotGenphenRfSvm,
plotSpecificGenotype
}



\details{
This procedure plots the results of the function runGenphenBayes. Each
association is shown as a point with respect to its Bayesian effect size
(mu.effect), including the corresponding HDI.
}



\value{
\item{plot}{ggplot plot object.}
}



\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



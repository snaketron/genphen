\name{runGenphenBayes}



\alias{runGenphenBayes}




\title{Conducting genetic association analysis with Bayesian inference,
using hierarchical models.}



\description{
This procedure performs a Bayesian two-sample t-test, to quantify the
association between genotypes (SNPs or SAAPs) and phenotypes.
}



\usage{
runGenphenBayes(genotype, phenotype, chain.nr, mcmc.iter, model, hdi.levels)

}


\arguments{
\item{genotype}{Character matrix or data frame, containing SNPs/SAAPs as
columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{phenotype}{Numerical vector, where each element is a measured phenotype
corresponding to the observations of the genotype data.}
\item{chain.nr}{Number of MCMC chains to simulate (recommended > 1).}
\item{mcmc.iter}{Number of MCMC iterations used in the Bayesian inference
(recommended >= 1000).}
\item{model}{The hierarchical model used for the Bayesian inference (currently
only two models are available). In case of a normally distributed phenotype,
the model 'ndist' should be used, and 'tdist' if it is t-distributed (more
robust to outliers).}
\item{hdi.levels}{Highest density intervals (HDI) (default = c(0.95, 0.99)).}
}



\details{
This procedure takes two types of data as input: first a genotype data (SNPs
or SAAPs), each represented as a categorical vector; second a phenotype
represented by a vector of continuous numbers, whose elements correspond to the
elements of the genotype data. It then performs a Bayesian two-sample t-test
analysis, computing the contrast between the phenotype means in the two groups
(mu.1-mu.2). The high density intervals (HDIs) of the posteriors are also
provided to quantify the reliability of the contrast.

In the classical two-sample t-test two assumptions are made: first, the
phenotypes in each group should be normally distributed, and second they should
have equal variances. Here, we have designed hierarchial models such that these
assumptions are relaxed. Thus, we allow for the phenotype in either group to
either be normally (model = 'ndist') or t-distributed (model = 'tdist'), while
allowing heterogeneous variances.
}




\value{
Five classes of results are computed for each genotype-phenotype association,
stored as a row in a data frame:
\item{site, g.1, g.1, count.1, count.1}{General information about the genotype.}
\item{mu.1, mu.2, sd.1, sd.2}{Mean posterior estimates for the mean
phenotype associated with each genotype (mu.1 and mu.2), as well as the mean
standard deviation posteriors sd.1 and sd.2.}
\item{mu.effect, mu.effect.95.L, mu.effect.95.H}{The bayesian effect size
mu.effect = mu.1 - mu.2, including the high density interval (HDI).}
\item{ess.mu.1, ess.mu.2}{Effective sampling size estimated for mu.1 and mu.2}
\item{t.test.pvalue}{P-value score from a classical two-sample T-test.}
}



\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



\examples{
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
}



\seealso{
runGenpenRf, runGenpenSvm, plotGenphenBayes,
plotSpecificGenotype, plotManhattan
}



\references{
  \itemize{
    \item Kruschke, John. Doing Bayesian data analysis: A tutorial with R, JAGS,
    and Stan. Academic Press, 2014.
  }
}


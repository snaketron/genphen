\name{runDiagnostics}



\alias{runDiagnostics}




\title{
Data reduction procedure
}



\description{
Running genphen can be computationally costly for hundreds of thousands of
predictors. Motivated by the biological fact that a large fraction of the the
associations in a genetic association study are weak, genphen implements a
computationally light-weight diagnostic procedure to allow the user discard some
of the weakest SNPs before the main association analysis.
}



\usage{
runDiagnostics(genotype, phenotype, phenotype.type, mcmc.chains,
               mcmc.iterations, mcmc.warmup, mcmc.cores, hdi.level,
               cv.iterations, diagnostics.points, diagnostics.samples,
               diagnostics.rf.trees)
}

\arguments{
\item{genotype}{Character matrix/data frame or a vector, containing SNPs/SAAPs
as columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{phenotype}{Numerical vector for continuous-phenotype analysis, numerical
or character vector for dichotonous-phenotype analysis.}
\item{phenotype.type}{'continuous' or 'dichotomous' based on phenotype type.}
\item{mcmc.chains}{Number of MCMC chains used to test each association test.
We recomend mcmc.chains >= 2.}
\item{mcmc.iterations}{Length of MCMC chains (default = 1,000).}
\item{mcmc.warmup}{Length of adaptive MCMC chains (default = 500).}
\item{mcmc.cores}{Number of cores used for the MCMC (default = 1).}
\item{hdi.level}{Highest density interval (HDI) (default = 0.95).}
\item{cv.iterations}{cross-validation iterations (default = 1,000).}
\item{diagnostics.points}{Number of diagnostics points to analyze along the
genotype space (default = 10).}
\item{diagnostics.samples}{Number of diagnostics samples to take from each
diagnostics point (default = 10).}
\item{diagnostics.rf.trees}{Number of random forest trees to use for the
variable importance analysis (default = 50,0000).}
}



\details{
Procedure:
1) Using random forest on complete genotype-phenotype data and infer variable
importance for each genotype.
2) Sort genotypes by importance, and sample few genotypes at different points
along the importance space, performing for each genotype the procedure explained
in runGenphen.
3) Visualize results which can help the user to determine whether a sensible
data reduction can be done, i.e. to select X number of most important genotypes
for the main analysis.
}

\value{
\bold{General parameters:}
\item{site}{id of the site (e.g. position in the provided sequence alignment)}
\item{mutation}{type of polymorphism (e.g. T->A)}
\item{data}{number of data points for each allele (e.g. A:10, T:20)}

\bold{Association score parameters:}
\item{cohens.d or absolute.d}{Cohen's d effect size (continuous phenotype
analysis) or absolute effect size (dichotomous phenotype analysis)
point estimate}
\item{cohens.d.L/cohens.d.H or absolute.d.L/absolute.d.H}{The highest density
interval (HDI) of the estimated effect size}
\item{ca, ca.L, ca.H}{Classification accuracy (CA) estimate and its HDI}

\bold{Ranked variable importance scores:}
\item{site}{id of the site (e.g. position in the provided sequence alignment)}
\item{importance}{magnitude of importance of the site}
\item{diagnostics.point}{indicator of selected diagnostics.points}
}


\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



\examples{
# I: Continuous diagnostics
# genotype inputs:
data(genotype.saap)
# phenotype inputs:
data(phenotype.saap)

# run genphen
continuous.diagnostics <- runDiagnostics(genotype = genotype.saap,
                                         phenotype = phenotype.saap,
                                         phenotype.type = "continuous",
                                         mcmc.chains = 2,
                                         mcmc.iterations = 1500,
                                         mcmc.warmup = 500,
                                         mcmc.cores = 2,
                                         hdi.level = 0.95,
                                         cv.iterations = 200,
                                         diagnostics.points = 3,
                                         diagnostics.samples = 3,
                                         diagnostics.rf.trees = 50000)
}

\seealso{
runGenphen
}

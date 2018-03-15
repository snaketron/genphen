\name{runDiagnostics}



\alias{runDiagnostics}




\title{Data reduction procedure}



\description{
Running genphen for hundreds of thousands of predictors (e.g. SNPs) can be
computationally costly. Motivated by the biological fact that most SNPs have no
or weak associations with the phenotype, genphen allows the user to run a
light-weight diagnostic procedure which allows the user to discard large portion
of the non-informative SNPs before running the main association analysis. The
data reduction step proceeds as follows: 1) using random forest and their
measures of variable importance, we obtain one importance value for each SNP.
We can use the spectrum of importances as a 'rough' guide to determine the
importance level at which the non-informative SNPs are dominant; 2) the user
can then select different points along the importance spectrum (e.g. ranks 1:5,
100:105, 1000:1005, ...) for which genphen is run using its standard procedure.
We can then use the association scores produced by genphen to determine the
importance rank at which the SNPs are no longer informative and thereby achieve
data reduction.
}



\usage{
runDiagnostics(genotype, phenotype, phenotype.type, rf.importance.trees,
               with.anchor.points, mcmc.chains, mcmc.iterations, mcmc.warmup,
               mcmc.cores, hdi.level, anchor.points)
}

\arguments{
\item{genotype}{Character matrix/data frame or a vector, containing SNPs/SAAPs
as columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{phenotype}{Numerical vector for continuous-phenotype analysis, numerical
or character vector for dichotonous-phenotype analysis.}
\item{phenotype.type}{'continuous' or 'dichotomous' based on phenotype type.}
\item{rf.importance.trees}{Number of random forest trees to use for the variable
importance analysis (default = 50,000).}
\item{with.anchor.points}{Boolean whether to run the complete diagnostics
procedure (TRUE), or only the random forest based importance estimation (FALSE)}
\item{mcmc.chains}{Number of MCMC chains used to test each association test.
We recomend mcmc.chains >= 2.}
\item{mcmc.iterations}{Length of MCMC chains (default = 1,000).}
\item{mcmc.warmup}{Length of adaptive MCMC chains (default = 500).}
\item{mcmc.cores}{Number of cores used for the MCMC (default = 1). The same
parameter is for multicore execution of the statistical learning procedures.}
\item{hdi.level}{Highest density interval (HDI) (default = 0.95).}
\item{anchor.points}{Vector of ranks (based on the importance measure) at which
to select the genotypes, for which the diagnostics will be run.}
}



\details{
Procedure:
1) Run random forest on the complete genotype-phenotype data and infer variable
importance for each genotype.
2) Sort genotypes by importance, and sample few genotypes at different points
along the importance spectrum, performing for each genotype the procedure
explained in runGenphen.
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
\item{bc}{Bhattacharyya coefficient, degree of overlap between the posterior
predicted distributions of the phenotype in the two alleles of a SNP (or two
amino acid states of an SAAP.}
\item{anchor.point}{Indicator of selected anchor.point}

\bold{Ranked variable importance scores:}
\item{site}{id of the site (e.g. position in the provided sequence alignment)}
\item{importance}{magnitude of importance of the site}
\item{importance.rank}{rank based on the importance}
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
                                         rf.importance.trees = 50000,
                                         with.anchor.points = TRUE,
                                         mcmc.chains = 2,
                                         mcmc.iterations = 1500,
                                         mcmc.warmup = 500,
                                         mcmc.cores = 2,
                                         hdi.level = 0.95,
                                         anchor.points = c(1:10))
}

\seealso{
runGenphen, runPhyloBiasCheck
}

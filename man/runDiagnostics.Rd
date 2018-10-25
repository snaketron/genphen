\name{runDiagnostics}



\alias{runDiagnostics}




\title{Data reduction procedure}



\description{
The methods implemented in genphen are statistically superior to the ones 
implemented by most classical (frequentist) tools for GWAS. The major challenge, 
however, is the substantially increased computational cost when analyzing the 
effects of hundreeds of thousands of SNPs. Inspired by the biological assumption 
that the major fraction of the studied SNPs are non-informative (genetic noise) 
with respect to the selected phenotype, various data reduction techniques can be 
implemented to quickly scan the SNP and discard a substantial portion of the the 
SNPs deemed to be clearly non-informative.
}

\usage{
runDiagnostics(genotype, phenotype, phenotype.type, rf.trees,
               mcmc.chains, mcmc.iterations, mcmc.warmup,
               cores, hdi.level, diagnostic.points, ...)
}

\arguments{
\item{genotype}{Character matrix/data frame or a vector, containing SNPs/SAAPs
as columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{phenotype}{Numerical vector for continuous-phenotype analysis, numerical
or character vector for dichotonous-phenotype analysis.}
\item{phenotype.type}{'continuous' or 'dichotomous' based on phenotype type.}
\item{rf.trees}{Number of random forest trees to grow for the variable
importance analysis (default = 50,000).}
\item{mcmc.chains}{Number of MCMC chains used to test each association test.
We recomend mcmc.chains >= 2.}
\item{mcmc.iterations}{Length of MCMC chains (default = 1,000).}
\item{mcmc.warmup}{Length of adaptive MCMC chains (default = 500).}
\item{cores}{Number of cores used for the MCMC (default = 1). The same
parameter is for multicore execution of the statistical learning procedures.}
\item{hdi.level}{Highest density interval (HDI) (default = 0.95).}
\item{diagnostic.points}{Vector of ranks (based on the importance measure) at 
which to select the genotypes, for which the diagnostics will be run.}
\item{...}{Optional parameters include adapt_delta: STAN configuration 
(default = 0.9); max_treedepth: STAN configuration (default = 10).}
}



\details{
The data reduction procedure includes the following steps:
\enumerate{
\item The complete data (independent variables = genotype and response = 
phenotype) is used to train a random forest (RF) model, in order to quantify 
the association between each SNP and the phenotype using the variable importance 
measure implemented in the R package ranger.
\item Next, we rank the SNPs according to their their importance, and plot the
distribution of variable importances, to get an insight into the structure of 
the importances values and potentially detect outlying patterns.
\item We create a set of 'diagnostic points', i.e. SNPs having a specific rank
in terms of their importance value to be evaluated using a light-weight version 
of the standard genphen approach.
\item By plotting the estimated association scores as a function of the rank of 
importance,  can roughly determine the importance rank at which the SNPs no 
longer carry any signal of association, and discard that portian from the data
before conducting the main analysis using the standard genphen method.
}
}

\value{
\bold{General parameters:}
\item{site}{id of the site (e.g. position in the provided sequence alignment)}
\item{g1}{first genotype of the SNP}
\item{g0}{second genotype of the SNP}


\bold{Association score parameters:}
\item{beta.mean, beta.sd, beta.L/beta.H}{Mean, standard deviation and HDI of 
the estimated slope coefficient}
\item{alpha.mean, alpha.sd, alpha.L/alpha.H}{Mean, standard deviation and HDI 
of the estimated intercept coefficient}
\item{sd.mean, sd.L/sd.H, nu.mean, nu.L/nu.H}{Standard deviation and degrees of
freedom from the noise distribution used when analyzing continuous phenotypes,
including HDIs of the estimated parameters}
\item{ca, ca.L, ca.H}{CA estimate and corresponding HDI}
\item{kappa, kappa.L, kappa.H}{Cohen's kappa and corresponding HDI}
\item{bc}{Bhattacharyya coefficient}


\bold{Ranked variable importance scores:}
\item{site}{id of the site (e.g. position in the provided sequence alignment)}
\item{importance}{Magnitude of importance (impurity) of the site, estimated 
with random forest implemented in R package ranger}
\item{importance.rank}{rank of sites based on their importance}
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

# run diagnostics
continuous.diagnostics <- runDiagnostics(genotype = genotype.saap,
                                         phenotype = phenotype.saap,
                                         phenotype.type = "continuous",
                                         rf.trees = 5000,
                                         mcmc.chains = 2,
                                         mcmc.iterations = 1500,
                                         mcmc.warmup = 500,
                                         cores = 2,
                                         hdi.level = 0.95,
                                         diagnostic.points = seq(from = 1, 
                                                                 to = 100, 
                                                                 by = 15))
}

\seealso{
runGenphen, runPhyloBiasCheck
}

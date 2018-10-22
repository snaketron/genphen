\name{runGenphen}



\alias{runGenphen}




\title{
Genetic association analysis using Bayesian inference and statistical learning
methods
}



\description{
Given a set of genotypes such as single nucleotide polymorphisms (SNPs) or
single amino acid polymorphisms (SAAPs) for a set of N individuals, and the
corresponding N phenotypes, genphen computes severaly genotype-phenotypes
association scores using Bayesian inference and statistical learning.
}



\usage{
runGenphen(genotype, phenotype, phenotype.type, model.type,
           mcmc.chains, mcmc.iterations, mcmc.warmup,
           cores, hdi.level, stat.learn.method,
           cv.iterations, rpa.iterations, with.stan.obj)
}


\arguments{
\item{genotype}{Character matrix/data frame or a vector, containing SNPs/SAAPs
as columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{phenotype}{Numerical vector for continuous-phenotype analysis, numerical
or character vector for dichotonous-phenotype analysis.}
\item{phenotype.type}{'continuous' or 'dichotomous' depending on the type of
the phenotype.}
\item{model.type}{'univariate' or 'hierarchical'}
\item{mcmc.chains}{Number of MCMC chains used to test each association test.
We recomend mcmc.chains >= 2.}
\item{mcmc.iterations}{Length of MCMC chains (default = 1,000).}
\item{mcmc.warmup}{Length of adaptive MCMC chains (default = 500).}
\item{cores}{Number of cores to use (default = 1).}
\item{hdi.level}{Highest density interval (HDI) (default = 0.95).}
\item{stat.learn.method}{Character parameter used to specify the statistical
learning method used in the analysis. Currently two methods are available:
random forest ('rf') and support vector machine ('svm'). For no statistical
learning select 'none'.}
\item{cv.iterations}{cross-validation iterations (default = 1,000).}
\item{rpa.iterations}{Number of retrospective power analysis (RPA) iterations 
to be computed. If rpa.iterations = 0, no RPA will be performed.}
\item{with.stan.obj}{Logical indicator specifying whether the used stan objects
should be returrned (default = FALSE, the stan objects can be large).}
}



\details{
Input:
\itemize{
\item{genotype}
genotype data (e.g. set of 1,000 SNPs found along the aligned genomes of
10 individuals) - provided in one of three possible input types: 
\itemize{
  \item character vector of length N if only a single SNP/SAAP is 
  ment to be analyzed, containing the genotypes of N individuals.
  \item character matrix with dimensions NxM whereby the M columns represent 
  different SNPs or SAAPs, and the N rows represent different individuals for 
  which we have some measured phenotype.
  \item AAMultipleAlignment or DNAMultipleAlignment object; if the genotype 
  data is a multiple sequence alignment composed of N sequences.
}
\item{phenotype}
phenotypes of N individuals in the form of a N-sized vector. The type of the
phenotype can either be continuous or dichotomous. Therefore, genphen has two
analysis modes for each situation. The main difference between them is the
design of the Bayesian hierarchical model which are used.
}

Goal: To quantify the association between each genotype and phenotype. With
genphen, we provide several measures of association:
\itemize{
\item{Classification accuracy (CA): CA measures the degree of accuracy with 
which one can classify (predict) the alleles of a SNP from the phenotype. If 
there exists a strong association between a particular SNP and the phenotype, 
one should be able to train a statistical model (using RF or SVM) which 
accurately classifies the two alleles of that SNP solely from the phenotype 
data (CA close to 1). Otherwise, the model should perform poorly, with the 
classification accuracy of the model being approximately similar to that of 
simple guessing (CA close to 0.5)}

\item{Cohen's kappa statistic: There is one pitfall where the CA estimate 
can be misleading, and this is the case when the analyzed SNP is composed of 
unevenly represented genetic states (alleles). For instance, the allele A of 
a given SNP is found in 90\% of the individuals, while the other allele T in 
only 10\%. Such an uneven composition of the alleles can lead to misleading 
results, i.e. even without proper learning the algorithm can produce a high 
$CA close to 0.9 simply by always predicting the dominant label. The kappa 
statistics is a quality metric, which is to be used together with CA. Cohen 
defines the following meaningful kappa intervals: [kappa<0]: ``no agreement'', 
[0.0-0.2]: ``slight agreement'', [0.2-0.4]: ``fair agreement'', [0.4-0.6]: 
``moderate agreement'', [0.6-0.8]: ``substantial agreement'' and [0.8-1.0]: 
``almost perfect agreement''.}

\item{Slope coefficient (beta): for each SNP we compute beta (the main effect 
in the generalized linear model inferred with Bayesian inference) which 
indicates the strength of the association between the genetic states of the 
particular SNP and the phenotype. We also compute 95\% (for instance) highest 
density interval (HDI) of beta, which is defined as the interval that covers 
a 95\% of the posterior distribution, with every point inside the interval 
having a higher credibility than any point outside it. Thus we can define an 
association as significant if the null-effect, i.e. beta = 0 lies outside 
the 95\% HDI. We report for each parameter in the above models its posterior 
mean and the 95\% HDI.}

\item{Bhattacharyya coefficient (BC): the paramters inferred during the Bayesian
inferrence can be used to simulate new data (posterior prediction), i.e. the 
distributions of the phenotype for each allele of a given SNP. We then used the 
Bhattacharyya coefficient to quantify the overlap between these distributions, 
with BC close to 0 for low overlap (strong association) and BC close to 1 for a
complete overlap (weak association).}
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


\bold{MCMC convergence parameters:}
\item{site, g1, g0, n1, n0}{s=site, g1=fist genotype, g0=second genotype, 
n1=number g1 observations, n0=number of g0 observations}
\item{par}{Parameter name in the generalized linear model}
\item{n_ess}{Effective sampling size}
\item{Rhat}{Potential scale reduction factor}


\bold{Retrospective power analysis (RPA) parameters:}
\item{site, g1, g0, n1, n0}{s=site, g1=fist genotype, g0=second genotype, 
n1=number g1 observations, n0=number of g0 observations}
\item{rpa.power.error}{Ratio of RPA runs in which a statistically insignificant
effect, i.e. slope coefficient with HDI which contain the null effect. To
estimate power of SNP: rpa.power = 1 - rpa.power.error}
\item{rpa.sign.error}{Ratio of RPA runs in which the slope coefficient has a 
different sign compared to the originally estimated slope}
\item{rpa.beta.mean}{The mean of the slope coefficients computed during the RPA}
\item{rpa.beta.sd}{The standard deviation of the slope coefficients computed 
during the RPA}
\item{rpa.N}{Number of total RPA iterations}


\bold{Posterior predictive check parameters:}
\item{s, g, n}{s=site, g=genotype, n=number of observations}
\item{predicted.mu.i, predicted.mu.j, real.mu.i, real.mu.j}{Predicted and
real means in each group of a SNP}


\bold{Stan object:}
\item{stan.obj}{Stan object.}
}


\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


\examples{
# I: Continuous phenotype analysis

# genotype inputs:
data(genotype.saap)
# phenotype inputs:
data(phenotype.saap)

# run genphen
continuous.analysis <- runGenphen(genotype = genotype.saap[, 1:3],
                                  phenotype = phenotype.saap,
                                  phenotype.type = "continuous",
                                  model.type = "univariate",
                                  mcmc.chains = 2,
                                  mcmc.iterations = 1500,
                                  mcmc.warmup = 500,
                                  cores = 2,
                                  hdi.level = 0.95,
                                  stat.learn.method = "rf",
                                  cv.iterations = 200,
                                  rpa.iterations = 0,
                                  with.stan.obj = FALSE)

# II: Dichotomous phenotype analysis

# genotype inputs:
data(genotype.saap)
# phenotype inputs:
data(dichotomous.phenotype.saap)

# run genphen
dichotomous.analysis <- runGenphen(genotype = genotype.saap[, 1:3],
                                   phenotype = dichotomous.phenotype.saap,
                                   phenotype.type = "dichotomous",
                                   model.type = "univariate",
                                   mcmc.chains = 2,
                                   mcmc.iterations = 1500,
                                   mcmc.warmup = 500,
                                   cores = 2,
                                   hdi.level = 0.95,
                                   stat.learn.method = "rf",
                                   cv.iterations = 200,
                                   rpa.iterations = 0,
                                   with.stan.obj = FALSE)
}


\seealso{
runDiagnostics, runPhyloBiasCheck
}

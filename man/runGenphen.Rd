\name{runGenphen}



\alias{runGenphen}




\title{
Genetic association analysis using Bayesian inference and statistical learning
methods
}



\description{
Given a set of genotypes (single nucleotide polymorphisms - SNPs; or single 
amino acid polymorphisms - SAAPs) for a set of individuals, and a corresponding
set of phenotypes, genphen quantifies the association between each genotype and
phenotype using Bayesian inference and statistical learning.
}



\usage{
runGenphen(genotype, phenotype, phenotype.type, model.type,
           mcmc.chains, mcmc.steps, mcmc.warmup, cores, 
           hdi.level, stat.learn.method, cv.steps, ...)
}


\arguments{
\item{genotype}{Character matrix/data frame or a vector, containing SNPs/SAAPs
as columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{phenotype}{Numerical vector (for a single phenotype) or matrix with 
multiple phenotypes stored as columns.}
\item{phenotype.type}{Vector representing the type of each phenotype (of the 
phenotype input), with 'Q' identifier for quantitative, or 'D' for dichotomous 
phenotypes.}
\item{model.type}{Type of Bayesian model: 'univariate' or 'hierarchical'}
\item{mcmc.chains}{Number of MCMC chains (default = 2).}
\item{mcmc.steps}{Length of MCMC chains (default = 1,000).}
\item{mcmc.warmup}{Length of adaptive part of MCMC chains (default = 500).}
\item{cores}{Number of cores to use (default = 1).}
\item{hdi.level}{Highest density interval (HDI) (default = 0.95).}
\item{stat.learn.method}{Parameter used to specify the statistical learning 
method used in the analysis. Currently two methods are available: random 
forest ('rf') and support vector machine ('svm'). For no statistical learning
select 'none'.}
\item{cv.steps}{cross-validation steps (default = 1,000).}
\item{...}{Optional parameters include adapt_delta: STAN configuration 
(default = 0.9); max_treedepth: STAN configuration (default = 10); ntree: Number 
of random forest trees to grow, only in case stat.learn.method = 'rf' (default 
= 1000); cv.fold: Cross-validation fold (default = 0.66).}
}






\details{
Input:
\itemize{
\item{genotype}
genotype data (e.g. set of 1,000 SNPs found along the aligned genomes of
10 individuals) - provided in one of three possible input types: 
\itemize{
  \item character vector of length N (if only a single SNP/SAAP is provided), 
  containing the genotypes of N individuals.
  \item character matrix with dimensions NxS (N = individuals, S = SNPs/SAAPs).
  \item AAMultipleAlignment or DNAMultipleAlignment object; if the genotype 
  data is a multiple sequence alignment composed of N sequences.
}
\item{phenotype}
phenotype data (dichotomous or quantitative phenotypes allowed)
\itemize{
  \item numerical vector of length N if only a single phenotype is analyzed
  \item numerical matrix NxP, if P phenotypes are provided.
}
\item{phenotype.type}
Vector with identifiers specifying the type of the phenotypes with 'Q' (for 
quantitative) or 'D' (for dichotomous) for each column in the phenotype dataset.

\item{model.type}
Specifies the structure of Bayesian model used to estimate the effect size of 
each genotype. Options allow for either 'univariate' (each SNP/SAAP treated as 
completely independent) or 'hierarchical' (SNP/SAAP effects share information 
through partial pooling).
}

Metrics:
To quantify the association between each genotype and phenotype genphen computes 
multiple measures of association:
\itemize{

\item{Effect size  (beta): for each SNP we compute beta (effect) with Bayesian
inference). beta quantifies the strength of the association between the 
genotypes and the phenotype. We report for each beta its mean and 95\% (for 
instance) highest density interval (HDI) of beta, which is defined as the 
interval that covers a 95\% of the posterior distribution, with every point 
inside the interval having a higher credibility than any point outside it.}

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
}
}

\value{
\bold{General parameters:}
\item{site}{id of the site (e.g. position in the provided sequence alignment)}
\item{ref, alt}{reference and alternative genotypey}
\item{refN, altN}{count of ref and alt genotypes}
\item{phenotype.id}{Identifier of the studied phenotype}


\bold{Association scores:}
\item{beta.mean, beta.se, beta.sd, beta.hdi.low/beta.hdi.high}{
Estimates of the mean, standard error, standard deviation and HDI of the
slope coefficient}
\item{ca.mean, ca.hdi.low/ca.hdi.high}{CA estimate and HDI}
\item{kappa.mean, kappa.hdi.low/kappa.hdi.high}{Cohen's kappa and HDI}
\item{rank}{Pareto optimiazion based front (rank) of SNP/SAAP estimated by 
maximizing metrics beta.mean and kappa.mean}

\bold{MCMC convergence parameters:}
\item{Neff}{Effective sampling size}
\item{Rhat}{Potential scale reduction factor}


\bold{Posterior predictions:}
\item{ppc}{Posterior prediction check and real data summary for each genotype.}


\bold{Posterior summary:}
\item{complete.posterior}{Complete stan object containing the posterior of each 
parameter estimated during the Bayesian inference. The data can be used for 
model debugging, posterior predictive checks, etc.}
}


\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}


\examples{
# genotypes:
data(genotype.saap)
# quantitative phenotype:
data(phenotype.saap)
# dichotomous phenotype:
data(dichotomous.phenotype.saap)
# make phenotype matrix (column = phenotype)
phenotypes <- cbind(phenotype.saap, dichotomous.phenotype.saap)

# run genphen
out <- runGenphen(genotype = genotype.saap[, 80:82],
                  phenotype = phenotypes,
                  phenotype.type = c("Q", "D"),
                  model.type = "univariate",
                  mcmc.chains = 4,
                  mcmc.steps = 1500,
                  mcmc.warmup = 500,
                  cores = 2,
                  hdi.level = 0.95,
                  stat.learn.method = "rf",
                  cv.steps = 200)
}


\seealso{
runDiagnostics, runPhyloBiasCheck
}

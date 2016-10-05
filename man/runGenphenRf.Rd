\name{runGenphenRf}



\alias{runGenphenRf}




\title{Conducting genetic association analysis with random forest}



\description{
This procedure quantifies the accuracy with which one can predict a given
genotypes (SNPs or SAAPs) from the corresponding phenotypes using random
forest.
}



\usage{
runGenphenRf(genotype, phenotype, cv.fold, cv.steps, hdi.level, ntree)
}


\arguments{
\item{genotype}{Character matrix or data frame, containing SNPs/SAAPs as
columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{phenotype}{Numerical vector, where each element is a measured phenotype
corresponding to the observations of the genotype data.}
\item{cv.fold}{The cross-validation fraction (0, 1) of the data which is used
to train the classifier (default = 0.66). The remaining fraction (1-cv.fold) of
the data is used to test the classifier.}
\item{cv.steps}{Number of steps in the cross-validation to be performed to
estimate the classification accuracy and the corresponding highest density
intervals(recommended >= 100).}
\item{hdi.level}{Highest density level (default = 0.99).}
\item{ntree}{Number of trees to grow in the random forest model.}
}



\details{
This procedure takes two types of data as input: first a genotype data (SNPs
or SAAPs), each represented as a categorical vector; second a phenotype
represented by a vector of continuous numbers, whose elements correspond to the
elements of the genotype data. It then performs a random forest classification
analysis, computing the classification accuracy (CA) with which one can predict
the categories of a given genotype from the numerical values of the phenotype
(predictor). If there is a significant association between a given genotype and
phenotype, then one should be able to accurately predict the states of the
genotype from that phenotype (CA = 1 for a perfect classification ). To obtain
a robust CA, this procedure applies cross-validation, training the models on a
subset of the genotype-phenotype data, and then testing the models on the
remaining data. Thus, in addition to the CA point estimate, one also estimates
the corresponding confidence intervals.

In some instances the two states of the genotype are not evenly represented,
i.e. one allele of a SNP might be represented in 90\% of the individuals, while
the other allele might only be represented in 10\% of the individuals. Such
data compositions can lead to skewed CA estimates. Therefore, to verify the
CA, the procedure also computes the Cohen's kappa statistic (Cohen 1960), which
compares the observed classification accuracy (CA) with the classification
accuracy expected by chance (CAexp). If CA > CAexp, the procedure will yield
high kappas (kappa = 1 in an ideal case), otherwise low kappas are to be
expected (kappa ~ 0).

The function runGenphenRf also computes statistics such as Cohen's d (effect
size) and the P-value resulting from a two-sample T-test, allowing the user to
compare the random forest based results with those computed with simpler
techniques which are frequently used for genetic association studies.
}




\value{
Five classes of results are computed for each genotype-phenotype association,
stored as a row in a data frame:
\item{site, g.1, g.1, count.1, count.1}{General information about the genotype.}
\item{ca, ca.hdi.low, ca.hdi.high, ca.hdi.length}{Mean classification accuracy
and its HDI.}
\item{kappa, kappa.hdi.low, kappa.hdi.high, kappa.hdi.length}{Cohen's kappa
statistics and its HDI.}
\item{effect.size, effect.CI.low, effect.CI.high}{Cohen's effect size and CI.}
\item{t.test.pvalue}{P-value score from a classical two-sample T-test.}
}



\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



\examples{
data(genotype.saap)
#or data(genotype.saap.msa) in this case you cannot subset genotype.saap[, 1:5]
data(phenotype.saap)
genphen.rf <- runGenphenRf(genotype = genotype.saap[, 1:5],
                           phenotype = phenotype.saap,
                           cv.fold = 0.66,
                           cv.steps = 100,
                           hdi.level = 0.99,
                           ntree = 1000)
}



\seealso{
runGenpenSvm, runGenpenBayesian, plotGenphenResults, plotSpecificGenotype,
plotManhattan
}




\references{
  \itemize{
    \item
    Breiman, Leo. Random forests. Machine learning 45.1 (2001): 5-32.
    \item Cohen, Jacob. Statistical power analysis for the behavior science.
    Lawrance Eribaum Association (1988).
    \item Cohen, Jacob. A coefficient of agreement for nominal scales (1960).
  }
}



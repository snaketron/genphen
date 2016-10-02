\name{runGenphenSvm}



\alias{runGenphenSvm}




\title{Conducting genetic association analysis with linear support vector
machines (LSVM)}



\description{
This procedure quantifies the accuracy with which one can predict a given
genotypes (SNPs or SAAPs) from the corresponding phenotypes using linear
support vector machines (LSVM).
}



\usage{
runGenphenSvm(genotype, phenotype, cv.fold, cv.steps, hdi.level)
}


\arguments{
\item{genotype}{Character matrix or data frame, containing SNPs/SAAPs as
columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{phenotype}{Numerical vector, where each element is a measured phenotype
corresponding to the observations of the genotype data.}
\item{cv.fold}{The cross-validation fraction (0, 1) of the data which is used
to train the classifier (recommended = 0.66). The remaining fraction
(1-cv.fold) of the data is used to test the classifier.}
\item{cv.steps}{Number of steps in the cross-validation to be performed to
estimate the classification accuracy and the corresponding highest density
intervals(recommended >= 100).}
\item{hdi.level}{Highest density interval (default = 0.99).}
}



\details{
This procedure takes two types of data as input: first a genotype data
composed of a set of single nucleotide polymorphisms (SNPs) or alternatively
single amino acid polymorphisms (SAAPs), each of which is represented by a
column of character amino acids; second a numerical phenotype vector, where
the elements sorted to correspond to the rows of the genotype data. This method
quantifies the association between the polymorphic site (SNP or SAAP) and the
phenotype via a classification analysis using linear support vector machines.
The analysis results in a classification accuracy score between 0 and 1, where
1 indicates a perfect association between the genotype and the phenotype. To
validate the classification accuracy, the tool also computes the Cohen's kappa
statistics (Cohen 1960) which compares the observed classification accuracy
with the expected classification accuracy. If the expected and observed
classification accuracies are in concordance, the computed association can be
taken seriously, otherwise it can be discarded as noise.

The function runGenphenSvm also computes statistics such as Cohen's d (effect
size) and the P-value resulting from a two-sample T-test, allowing the user to
compare the linear support vector based results with those computed with
simpler techniques which are frequently used for genetic association studies.
}




\value{
Five classes of results are computed for each SAAP with respect
to the phenotype, resulting in a 18 element vector which is stored
as a row in the final data frame:
\item{effect.size, effect.CI.low, effect.CI.high}{Cohen's effect
size and CI.}
\item{ca, ca.hdi.low, ca.hdi.high, ca.hdi.length}{Mean classification
accuracy and its HDI.}
\item{kappa, kappa.hdi.low, kappa.hdi.high, kappa.hdi.length}{Cohen's
kappa statistics and its HDI.}
\item{site, g.1, g.2, count.1, count.2}{General information
about the genotype.}
\item{t.test.pvalue}{P-value score from an two-sample T-test.}
}



\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



\examples{
data(genotype.saap)
#or data(genotype.saap.msa) in this case you cannot subset genotype.saap[, 1:5]
data(phenotype.saap)
genphen.svm <- runGenpenSvm(genotype = genotype.saap[, 1:5],
                            phenotype = phenotype.saap,
                            cv.fold = 0.66,
                            cv.steps = 100,
                            hdi.level = 0.99)
}



\seealso{
runGenpenRf, runGenpenBayes, plotGenphenRfSvm, plotGenphenBayes,
plotSpecificGenotype, plotManhattan
}





\references{
  \itemize{
    \item
    Cortes, Corinna, and Vladimir Vapnik. Support-vector networks. Machine
    learning 20.3 (1995): 273-297.
    \item Cohen, Jacob. Statistical power analysis for the behavior science.
    Lawrance Eribaum Association (1988).
    \item Cohen, Jacob. A coefficient of agreement for nominal scales (1960).
  }
}

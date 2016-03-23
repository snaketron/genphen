\name{runGenphenSnp}




\alias{runGenphenSnp}




\title{Performing genetic association analysis between SNPs and phenotypes}




\description{
This procedure computes the association between single nucleotide polymorphisms
(SNPs) and phenotypes.
}




\usage{
runGenphenSnp(genotype, phenotype, technique, fold.cv, boots)
}




\arguments{
\item{genotype}{Character matrix or data frame, containing SNPs as columns
or alternatively a DNAMultipleAlignment Biostrings object}
\item{phenotype}{Numerical vector, where each element is a measured phenotype
corresponding to the observations of the genotype data.}
\item{technique}{Two techniques are provided: random forests (rf) or linear
support vector machines (svm) (recommended = svm).}
\item{fold.cv}{The cross-validation fraction (0, 1) of the data which is used
to train the classifier (recommended = 0.66). The ramaining fraction
(1-fold.cv) of the data is used to test the classifier.}
\item{boots}{Number of bootstraps to be performed to estimate the
classification accuracy and the corresponding confidence intervals
(recommended >= 100).}
}




\details{
This procedure takes as an input two types of data: first a genotype data
composed of single nucleotide polymorphism (SNP) sites, each of which is
represented by a column of alleles, whereby at most two types of alleles
should exist in each column; second a numerical phenotype vector, where the
elements sorted to correspond to the rows of the genotype data.

Using these two data types, it computes the association between each SNP and
the phenotype. For each SNP two metrics are computed, called "effect size" and
"classification accuracy".

The effect size of a given SNP is obtained by computing the Cohen's d
statistics (Cohen 1988). The 95\% confidence intervals are computed as well.

Classification accuracy is the second metric which is computed using
statistical learning techniques. This is the metric which is used to quantify
the strength of the association between a SNP and a phenotype. The idea is to
use either linear suppport vector machines or random forests to build a
classification model between the phenotype vector and the SNP vector. The
more accurate the model, the easier we can predict the two allele states of
the SNP from the phenotype and hence the stronger is the mutual association
between the two vectors. In order to obtain a robust classification accuracy
measure, the classification analysis is done in a bootstrapping fashion.
First a subset of the SNP-phenotype vectors is randomly selected to train a
classifier, while the remaining data is used to test the classifier. This
step is repeated multiple times after which the classification accuracies
of all the classifiers are averaged into a single classification accuracy
measure and the corresponding confidence intervals are computed.

In order to validate the classification accuracy, the tool also computes
the Cohen's kappa statistics (Cohen 1960) which compares the observed
classification accuracy with the expected classification accuracy. If the
expected and observed classification accuracies are in concordance, the
computed association can be taken seriously, otherwise it can be discarded
as noise.
}





\value{
Five classes of results are computed for each SNP with respect to the
phenotype, resulting in a 18 element vector which is stored as a row in
the final data frame:
\item{effect.size, effect.CI.low, effect.CI.high}{Cohen's effect size
and 95\% CI.}
\item{ca, ca.CI.low, ca.CI.high, ca.CI.length}{Mean classification
accuracy and its 95\% CI.}
\item{kappa, kappa.CI.low, kappa.CI.high, kappa.CI.length}{Cohen's
kappa statistics and its 95\% CI.}
\item{site, allele1, allele2, count.allele1, count.allele2}{General
information about the genotype.}
\item{anova.score}{P-value score from an ANOVA test.}
}




\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



\examples{
data(genotype.snp)
#or data(genotype.snp.msa) in this case you cannot subset genotype.snp[, 1:3]
data(phenotype.snp)
genphen.results <- runGenphenSnp(genotype = genotype.snp[, 1:3],
phenotype = phenotype.snp, technique = "svm", fold.cv = 0.66, boots = 100)
}


\seealso{
runGenphenSaap, plotGenphenResults, plotSpecificGenotype, plotManhattan
}



\references{
Cohen, J. (1988) Statistical power analysis for the behavioral sciences
(2nd ed.). New York:Academic Press.

Cohen, J. (1960) A coefficient of agreement for nominal scales.
}

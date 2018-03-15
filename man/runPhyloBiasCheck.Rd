\name{runPhyloBiasCheck}



\alias{runPhyloBiasCheck}




\title{
Check for phylogenetic bias
}



\description{
Given a set of genotypes such as single nucleotide polymorphisms (SNPs) or
single amino acid polymorphisms (SAAPs) for a set of N individuals, the
procedure can operate in two modes:

1) it computes a NxN kinship matrix (matrix populated with pairwise distances
(Hamming) between each two individuals computed using all the genotypes). Based
on the kinship matrix it then estimates the degree of phlyogenetic bias related
to each genotype as 1-mean.phylo.dist(allele)/mean.phylo.dist(all)

2) it uses a precomputed kinship matrix and then estimates the degree of
phlyogenetic bias related to each genotype using the same procedure.
}


\usage{
runPhyloBiasCheck(input.kinship.matrix, genotype)
}


\arguments{
\item{genotype}{Character matrix/data frame or a vector, containing SNPs/SAAPs
as columns or alternatively as DNAMultipleAlignment or AAMultipleAlignment
Biostrings object.}
\item{input.kinship.matrix}{precomputed kinship matrix provided by the user.}
}



\details{
Input:
\itemize{
\item{genotype}
P genotypes of N individuals in the form of NxP character matrix/data frame or
vector (if P = 1).
\item{input.kinship.matrix}
precomputed NxN matrix (row/column for each individual)
}
}

\value{
\bold{Genotype parameters:}
\item{site}{id of the site (e.g. position in the provided sequence alignment)}
\item{genotype}{allele of a SNP or amino acid of SAAP}
\item{bias}{number between 0 (no bias) or 1 (complete bias)}

\bold{Mutation bias:}
\item{site}{id of the site (e.g. position in the provided sequence alignment)}
\item{mutation}{allele of a SNP or amino acid of SAAP}
\item{bias}{number between 0 (no bias) or 1 (complete bias) for the mutation
computed as max(bias in genotype 1, bias in genotype 2)}

\bold{Kinship matrix:}
\item{kinship.matrix}{NxN matrix}
}


\author{
  Simo Kitanovski <simo.kitanovski@uni-due.de>
}



\examples{
# genotype inputs:
data(genotype.saap)
# phenotype inputs:
data(phenotype.saap)

# phylogenetic bias analysis
bias <- runPhyloBiasCheck(input.kinship.matrix = NULL,
                          genotype = genotype.saap)
}


\seealso{
runDiagnostics, runGenphen
}

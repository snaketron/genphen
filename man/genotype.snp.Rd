\name{genotype.snp}


\alias{genotype.snp}


\docType{data}


\title{
SNP genotype dataset
}


\description{
The genotype.snp data is a character matrix with dimensions 51x100.
It contains 100 SNPs across 51 mouse strains, taken from the publicly
available Mouse Hapmap data. We used it in combination with the
phenotype.snp data to compute the association between each SNP and
the phenotype data.
}

\usage{data(genotype.snp)}


\format{
A matrix with 51 observations (laboratory mouse strains) and 100
variables (SNPs).
}

\source{
http://mouse.cs.ucla.edu/mousehapmap/emma.html
}

\value{
Matrix with 51 rows and 100 columns, whereby each column is a SNP,
and the elements represent an alleles (nucleotides).
}



\examples{
data(genotype.snp)
}

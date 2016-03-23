\name{genotype.saap}


\alias{genotype.saap}


\docType{data}


\title{
SAAP genotype dataset
}


\description{
The genotype.saap data is a character matrix with dimensions 120x154.
It contains 154 amino acid protein sites across 120 organisms. The data
is used in combination with the phenotype.aa data to quantify the
association between each amino acid substitution pair and the phenotype
vector.
}

\usage{data(genotype.saap)}


\format{
A matrix with 120 observations and 154 columns (some of which qualify
as single amino acid polymorphisms).
}

\source{
http://www.ncbi.nlm.nih.gov/genbank/
}



\value{
Matrix with 120 rows and 154 columns, whereby each row is a protein
sequence and the elements represent an amino acids.
}



\examples{
data(genotype.saap)
}

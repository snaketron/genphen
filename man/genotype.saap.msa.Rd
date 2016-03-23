\name{genotype.saap.msa}



\alias{genotype.saap.msa}



\docType{data}



\title{
SAAP genotype dataset (msa)
}



\description{
The genotype.saap.msa data is a multiple sequence alignment in Biostrings
AAMultipleAlignment format. It contains 120 protein sequences, each with 154
sites (SAAPs). The data is used in combination with the phenotype.aa data to
quantify the association between each amino acid substitution pair and the
phenotype vector.
}



\usage{
data("genotype.saap.msa")
}



\format{
AAMultipleAlignment object with 120 sequences each made of 154 amino acid
sites (SNPs), some of which qualify as single amino acid polymorphisms.
}



\source{
http://www.ncbi.nlm.nih.gov/genbank/
}



\value{
AAMultipleAlignment object with 120 sequences each made of 154 amino acid
sites (SNPs), some of which qualify as single amino acid polymorphisms.
}



\examples{
data("genotype.saap.msa")
}

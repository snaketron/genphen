\name{genotype.snp.msa}


\alias{genotype.snp.msa}


\docType{data}


\title{
SNP genotype dataset (msa)
}


\description{
The genotype.snp.msa data is a multiple sequence alignment in Biostrings
DNAMultipleAlignment format. It contains 51 DNA sequences, each with 100
sites (SNPs), taken from the publicly available Mouse Hapmap data. We used
it in combination with the phenotype.snp data to compute the association
between each SNP and the phenotype data.
}

\usage{
data("genotype.snp.msa")
}


\format{
DNAMultipleAlignment object with 51 sequences each made of 100 nucleotides
(SNPs).
}

\source{
http://mouse.cs.ucla.edu/mousehapmap/emma.html
}

\value{
DNAMultipleAlignment object with 51 sequences each made of 100 nucleotides
(SNPs).
}


\examples{
data("genotype.snp.msa")
}

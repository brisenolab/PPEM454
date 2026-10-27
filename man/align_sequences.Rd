\name{align_sequences}
\alias{align_sequences}
\title{Align nucleotide sequences}
\usage{align_sequences(sequences, method = c("ClustalW", "ClustalOmega", "Muscle"))}
\arguments{\item{sequences}{Character vector with at least two nucleotide sequences.}\item{method}{Alignment method accepted by msa.}}
\value{A named character vector of aligned sequences.}
\description{Requires optional Biostrings and msa packages. Converts U to T and removes existing gaps on an in-memory copy before alignment.}

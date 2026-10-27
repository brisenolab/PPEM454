\name{get_sequences}
\alias{get_sequences}
\title{Load supplied sequences}
\usage{get_sequences(host, n = 8L, data_dir = NULL)}
\arguments{\item{host}{A host accepted by choose_host.}\item{n}{Maximum records, an integer from 2 to 8.}\item{data_dir}{Optional folder containing supplied FASTA files.}}
\value{A named character vector.}
\description{Read the first records from the supplied host FASTA. Original symbols are preserved. Host mapping follows supplied filenames and needs instructor review.}

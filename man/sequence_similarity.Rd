\name{sequence_similarity}
\alias{sequence_similarity}
\title{Summarize sequence identity}
\usage{sequence_similarity(alignment)}
\arguments{\item{alignment}{Character vector of at least two equal-length aligned sequences.}}
\value{A list containing pairwise percentage identity, overall mean, min, and max.}
\description{Only columns with canonical bases in both sequences are compared; U is treated as T. Gaps and ambiguities are excluded pairwise. Overall is the unweighted mean across unique pairs. Pairs without comparable sites are NA and excluded from summaries.}

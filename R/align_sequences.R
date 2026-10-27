#' Align nucleotide sequences with the optional msa package
#' @param sequences Named character vector of nucleotide sequences.
#' @param method An msa alignment method.
#' @return Named character vector of aligned sequences.
#' @export
align_sequences <- function(sequences, method = c("ClustalW", "ClustalOmega", "Muscle")) {
  method <- match.arg(method)
  if (!is.character(sequences) || length(sequences) < 2L ||
      anyNA(sequences) || any(!nzchar(sequences))) stop("Supply at least two nonempty sequences.")
  if (!requireNamespace("Biostrings", quietly = TRUE) ||
      !requireNamespace("msa", quietly = TRUE)) {
    stop("Install Biostrings and msa with BiocManager before aligning sequences.")
  }
  # Normalize RNA U to DNA T and remove supplied gaps for realignment.
  # Original FASTA files and get_sequences() output remain unchanged.
  dna <- chartr("U", "T", toupper(gsub("-", "", sequences, fixed = TRUE)))
  aligned <- msa::msa(Biostrings::DNAStringSet(dna), method = method,
                      order = "input")
  as.character(methods::as(aligned, "DNAStringSet"))
}

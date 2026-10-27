#' Read supplied host-specific FASTA sequences
#' @param host A host accepted by choose_host().
#' @param n Maximum sequences to return, between 2 and 8 (first records).
#' @param data_dir Optional folder containing the original FASTA files.
#' @return A named character vector, with original sequence symbols preserved.
#' @export
get_sequences <- function(host, n = 8L, data_dir = NULL) {
  host <- choose_host(host)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n < 2 || n > 8 || n != floor(n)) stop("n must be an integer from 2 to 8.")
  filenames <- c(human = "human", animal = "animal", fungus = "fungal",
                 bacteria = "phage", plant = "plant")
  filename <- paste0(filenames[[host]], "_virus_sequences.fasta")
  path <- if (is.null(data_dir)) system.file("extdata", filename, package = "ViralExamR") else file.path(data_dir, filename)
  if (!nzchar(path) || !file.exists(path)) stop("FASTA not found; install the package or supply data_dir.")
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, ";")]
  headers <- which(startsWith(lines, ">"))
  if (length(headers) < 2L || headers[1L] != 1L) stop("Expected at least two FASTA records.")
  ends <- c(headers[-1L] - 1L, length(lines))
  seqs <- vapply(seq_along(headers), function(i) {
    if (ends[i] <= headers[i]) stop("Empty FASTA record.")
    paste0(lines[seq.int(headers[i] + 1L, ends[i])], collapse = "")
  }, character(1))
  names(seqs) <- make.unique(substring(lines[headers], 2L))
  seqs[seq_len(min(n, length(seqs)))]
}

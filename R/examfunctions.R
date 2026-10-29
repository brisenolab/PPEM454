#' Summarize pairwise identity in an alignment
#' @param alignment Character vector of equal-length aligned sequences.
#' @return List with pairwise percentage identity matrix and overall, min,
#'   and max percentages. Overall is the unweighted mean across unique pairs.
#' @details Each pair uses only columns where both symbols are A, C, G, T,
#'   or U (U is treated as T). Gaps and ambiguity symbols are excluded.
#'   Pairs with no comparable sites are NA and excluded from summaries.
#' @export
sequence_similarity <- function(alignment) {
  if (!is.character(alignment) || length(alignment) < 2L || anyNA(alignment) ||
      any(!nzchar(alignment)) || length(unique(nchar(alignment))) != 1L) {
    stop("Supply at least two nonempty, equal-length aligned sequences.")
  }
  chars <- strsplit(chartr("U", "T", toupper(alignment)), "", fixed = TRUE)
  ids <- names(alignment)
  if (is.null(ids)) ids <- paste0("sequence_", seq_along(chars))
  result <- matrix(NA_real_, length(chars), length(chars), dimnames = list(ids, ids))
  for (i in seq_along(chars)) for (j in seq_len(i)) {
    valid <- chars[[i]] %in% c("A", "C", "G", "T") & chars[[j]] %in% c("A", "C", "G", "T")
    value <- if (any(valid)) 100 * mean(chars[[i]][valid] == chars[[j]][valid]) else NA_real_
    result[i, j] <- result[j, i] <- value
  }
  pairs <- result[upper.tri(result)]
  pairs <- pairs[!is.na(pairs)]
  list(pairwise = result, overall = if (length(pairs)) mean(pairs) else NA_real_,
       min = if (length(pairs)) min(pairs) else NA_real_,
       max = if (length(pairs)) max(pairs) else NA_real_)
}

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

#' k-mer sequence partition
#' @param k is to select the k-mer size
#' @export

# function to count k-mers
count_kmers <- function(sequences, k) {
  if (length(k) != 1L || !is.numeric(k) || is.na(k) ||
      !is.finite(k) || k < 1 || k != floor(k)) {
    stop("k must be a single positive integer")
  }
  
  # Also support Biostrings sequence sets
  if (inherits(sequences, "XStringSet")) {
    sequences <- as.character(sequences)
  }
  
  if (!is.character(sequences) && !is.list(sequences)) {
    stop("sequences must be a character vector or a list")
  }
  
  lapply(sequences, function(sequence) {
    # Support sequences stored as vectors of individual bases
    if (!is.character(sequence) || anyNA(sequence)) {
      stop("Each sequence must contain character data without NA")
    }
    
    sequence <- toupper(paste0(sequence, collapse = ""))
    n <- nchar(sequence)
    
    # Sequences shorter than k contain no k-mers
    if (n < k) {
      return(table(character(0)))
    }
    
    starts <- seq_len(n - k + 1L)
    kmers <- substring(sequence, starts, starts + k - 1L)
    
    #table(kmers, dnn = paste0(k, "-mers")) 
    return(table(kmers, dnn = paste0(k, "-mers"))) 
    })
}

#' Create an incidence matrix
#' @param set.seed to select a fix seed number for the incidence matrix
#' @param nrow to select the number of rows
#' @param ncol to select the number of columns
#' @return An incidence matrices with viruses and hosts
#' @export
#' 
create_incidence_matrix <- function(nrow = 10, ncol = 8, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  values <- rnorm(nrow * ncol, mean = 40, sd = 60) |>
    round(0)
  
  values[values < 0] <- 0
  
  mat <- matrix(values, nrow = nrow, ncol = ncol)
  rownames(mat) <- paste0("virus", seq_len(nrow))
  colnames(mat) <- paste0("sample", LETTERS[seq_len(ncol)])
  
  return(mat)
}

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

#' Choose a host for the sequence exercise
#' @param host One of human, animal, fungus, bacteria, or plant. If NULL,
#'   prompt in an interactive R session.
#' @return A host name as a character scalar.
#' @export
choose_host <- function(host = NULL) {
  hosts <- c("human", "animal", "fungus", "bacteria", "plant")
  if (is.null(host)) {
    if (!interactive()) stop("Supply a host in non-interactive sessions.")
    choice <- utils::menu(hosts, title = "Choose your host")
    if (choice == 0L) stop("Host selection cancelled.")
    return(hosts[[choice]])
  }
  if (!is.character(host) || length(host) != 1L || is.na(host) ||
      !tolower(host) %in% hosts) stop("Choose human, animal, fungus, bacteria, or plant.")
  tolower(host)
}



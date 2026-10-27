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

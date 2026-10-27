# Example viral nucleotide sequence
sequence <- "ATGCGATCGATCGATGCTAGCTAGCTAGCATCGATCGATCG"


# function to count k-mers
count_kmers <- function(sequence, k) {
  
  sequence <- toupper(sequence)
  
  if (k < 1 || k > nchar(sequence)) {
    stop("k must be between 1 and the sequence length")
  }
  
  # Extract overlapping k-mers
  kmers <- sapply(1:(nchar(sequence) - k + 1), function(i) {
    substr(sequence, i, i + k - 1)
  })
  
  # Count each unique k-mer
  frequencies <- table(kmers)
  
  return(frequencies)
}

# Select k-mer size
k <- 6

results <- count_kmers(sequence, k)

print(results)


barplot(
  results,
  main = paste("k-mer distribution (k =", k, ")"),
  xlab = "k-mer",
  ylab = "Frequency",
  las = 2,
  col = "steelblue",
  border = NA
)


par(mfrow = c(2, 2))

for (k in 2:5) {
  
  results <- count_kmers(sequence, k)
  
  barplot(
    results,
    main = paste("k =", k),
    xlab = "k-mer",
    ylab = "Frequency",
    las = 2,
    col = "steelblue",
    border = NA,
    cex.names = 0.6
  )
}

par(mfrow = c(1, 1))

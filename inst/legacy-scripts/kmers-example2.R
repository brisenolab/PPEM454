
# Viral sequences
virus_A <- "CGGAACTAAACTCGTGGTTCCTGTGGTTCACACCTGACCTCCTGAGCAGAAAAGAAAAAAGAAGGCGGCTCGGAGGAGCGCTTCAGGGATCCCCGGGGAAACCTGGAGCGAACTGGCAAAAAAGGACGGTGGGGAGTGCCCAGCGGCCGACAGGAGTAATTCCCGCCGAAACAGGGTTTTCACCCTTCCTTTCTTCGGGTGTCCTTCCTCGCGCCCGCAGGACCACCCCTCGCCCCCTTTGCGCTGTCGCTTCGGCTACTACCCGGTGGAAACAACTGAAGCTCCCGAGAACCGCTTTTTCTCTATCTTACTTGCTTCGGGGCGAGGGTGTTTAGCCCTTGGAACCGCAGTTGGTTCCT"
virus_B <- "TTTATTAGAACAAGAAGTGAGGATATGATTAAACTTTGTTTGACGAAACCAGGTCTGTTCCGACTTTCCGACTCTGAGTTTCGACTTGTGAGAGAAGGAGGAGTCGTGGTGAACTTTTATTAAAAAAATTAGTTCACTCGTCTTCAATCTCTTGATCACTTCGTCTCTTCAGGGAAAGATGGGAAGAACACTGATGAGTCTCGCAAGGTTTACTCCTCTATCTTCATTGTTTTTTTACAAAATCTTG"

# Select k-mer size
k <- 4

# Calculate frequencies
freq_A <- count_kmers(virus_A, k)
freq_B <- count_kmers(virus_B, k)

# Combine all observed k-mers
all_kmers <- sort(union(names(freq_A), names(freq_B)))

# Create frequency matrix
A <- setNames(numeric(length(all_kmers)), all_kmers)
B <- A

A[names(freq_A)] <- as.numeric(freq_A)
B[names(freq_B)] <- as.numeric(freq_B)

# Plot frequency distributions
barplot(
  rbind(A, B),
  beside = TRUE,
  col = c("steelblue", "orange"),
  names.arg = all_kmers,
  las = 2,
  cex.names = 0.8,
  xlab = "k-mer",
  ylab = "Frequency (counts)",
  main = paste("k-mer frequency comparison (k =", k, ")"),
  legend.text = c("Virus A", "Virus B"),
  args.legend = list(x = "topright", bty = "n")
)


# # Install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # Install with BiocManager
# BiocManager::install("msa")
# # Install with R base
# install.packages(c("ape", "phangorn"))

#---------------------------
library(Biostrings)

# Create the seqTences:
sequences <- DNAStringSet(c(
  "SL-CoVZC45" = "GGTGAATTTAAATTGTCTTCACATATGTACTGTTCTTTCTACCCTCCTGAAGATGAAGGGGAAGATGATTGTGAAGAAGGACAGTGTGAACCATCAACTCAATATGAGTATGGTACTGAGGATGACTACCAAGGTAAACCTTTGGAGTTTGGTGCTACTTCTTTTTCTTCTTCTTCACAGGAAGAAGAACAAGAAGAGGATTGGTTAGAATCTGATAGTCAGGAC",
  "ASino1-11" = "GGTGAAGAAAACTTTTCATCACGTATGTATTGTTCCTTTTACCCTCCAGATGAGGAAGAAGAGGACGATGCAGAGTGTGAGGAAGAAGAAATTGATGAAACCTGTGAACATGAGTACGGTACAGAGGATGATTATCAAGGTCTCCCTCTGGAATTTGGTGCCTCAGCTGAAACAGTTCGAGTTGAGGAAGAAGAAGAGGAAGACTGGCTGGATGATACTACTGAGCAA",
  "NS-1" = "GGTGAAGAAAACTTTTCATCACGTATGTATTGTTCCTTTTACCCTCCAGATGAGGAAGAAGAGGACGATGCAGAGTGTGAGGAAGAAGAAATTGATGAAACCTGTGAACATGAGTACGGTACAGAGGATGATTATCAAGGTCTCCCTCTGGAATTTGGTGCCTCAGCTGAAACAGTTCGAGTTGAGGAAGAAGAAGAGGAAGACTGGCTGGATGATACTACTGAGCAA",
  "GD01" = "GGTGAAGAAAACTTTTCATCACGTATGTATTGTTCCTTTTACCCTCCAGATGAGGAAGAAGAGGACGATGCAGAGTGTGAGGAAGAAGAAATTGATGAAACCTGTGAACATGAGTACGGTACAGAGGATGATTATCAAGGTCTCCCTCTGGAATTTGGTGCCTCAGCTGAAACAGTTCGAGTTGAGGAAGAAGAAGAGGAAGACTGGCTGGATGATACTACTGAGCAA",
  "AGZ02" = "GGTGAAGAAAACTTTTCATCACGTATGTATTGTTCCTTTTACCCTCCAGATGAGGAAGAAGAGGACGATGCAGAGTGTGAGGAAGAAGAAATTGATGAAACCTGTGAACATGAGTACGGTACAGAGGATGATTATCAAGGTCTCCCTCTGGAATTTGGTGCCTCAGCTGAAACAGTTCGAGTTGAGGAAGAAGAAGAGGAAGACTGGCTGGATGATACTACTGAGCAA",
  "BtRs-BetaCoV/YN2018B" = "GGTGAAGAAAAACTTTCTTCACGTATGTACTGTTCCTTTTATCCTCCCGATGATGAGGAGGATTGTGATGAGTATGAGGAAGAAGAGGAAGTCCTGGAAGAATCCTGTGCGCATGAATACGGTACAGAAGAAGACTACCAAGGTCTTCCACTGGAATTTGGTGCCTCAACTGAAATGCAAGTTGAGGAAGAAGAAGAAGAGGACTGGCTTGGTGATGCTACCGAATTA",
  "Rs7327" = "GGTGAAGAAAAACTTTCTTCACGTATGTACTGTTCCTTCTATCCTCCTGATGATGAGGAGGATTGTGATGAGTATGAGGAAGAAGAGGAAGTCCCGGAAGAATCCTGTGCGCATGAATACGGTACAGAAGAAGACTACCAAGGTCTTCCACTGGAATTTGGTGCCTCAACTGAAATGCAAGTTGAAGAAGAAGAAGAAGAGGACTGGCTTGGTGATGCTACCGAATTA",
  "Rs9401" = "GGTGAAGAAAAACTTTCTTCACGTATGTACTGTTCCTTTTATCCTCCTGATGATGAGGAGGATTGTGATGAGTATGAGGAAGAAGAGGAAGTCCCGGAAGAATCCTGTGCGCATGAATACGGTACAGAAGAAGACTACCGAGGTCTTCCACTGGAATTTGGTGCCTCAACTGAAATGCAAGTTGAAGAAGAAGAAGAAGAGGACTGGCTTGGTGATGCTACCGACTTA",
  "WIV16" = "GGTGAAGAAAAACTTTCTTCACGTATGTACTGTTCCTTCTATCCTCCTGATGATGAGGAGGATTGTGATGAGTATGAGGAAGAAGAGGAAGTCCCGGAAGAATCCTGTGCGCATGAATACGGTACAGAAGAAGACTACCAAGGTCTTCTACTGGAATTTGGTGCCTCAACTGAAATGCAAGTTGAAGAAGAAGAAGAAGAGGACTGGCTTGGTGATGCTACCGAATTA",
  "Rs4231" = "GGTGAAGAAAAACTTTCTTCACGTATGTACTGTTCCTTCTATCCTCCTGATGATGAGGAGGATTGTGATGAGTATGAGGAAGAAGAGGAAGTCCCGGAAGAATCCTGTGCGCATGAATACGGTACAGAAGAAGACTACCAAGGTCTTCCACTGGAATTTGGTGCCTCAACTGAAATGCAAGTTGAAGAAGAAGAAGAAGAGGACTGGCTTGGTGATGCTACTGAATTA"
))

# Inspect seqTences
sequences

# Function to introduce random substitutions
mutate_dna <- function(seq, n_mutations) {
  
  bases <- strsplit(as.character(seq), "")[[1]]
  
  positions <- sample(
    seq_along(bases),
    n_mutations,
    replace = FALSE
  )
  
  for (pos in positions) {
    bases[pos] <- sample(
      setdiff(c("A", "C", "G", "T"), bases[pos]),
      1
    )
  }
  
  paste0(bases, collapse = "")
}


#--------------------
# Add a little of artificial variation to the original sequences
# Group 1: SARS-related sequences
ancestor_1 <- mutate_dna(
  as.character(sequences["ASino1-11"]),
  12
)

# Group 2: Bat coronavirus sequences
ancestor_2 <- mutate_dna(
  as.character(sequences["Rs7327"]),
  12
)

# Generate modified sequences
sequences_mutated <- DNAStringSet(c(
  
  # More divergent sequence
  "SL-CoVZC45" = mutate_dna(
    as.character(sequences["SL-CoVZC45"]), 20
  ),
  
  # Group 1
  "ASino1-11" = mutate_dna(ancestor_1, 3),
  "NS-1"      = mutate_dna(ancestor_1, 5),
  "GD01"      = mutate_dna(ancestor_1, 8),
  "AGZ02"     = mutate_dna(ancestor_1, 10),
  
  # Group 2
  "BtRs-BetaCoV/YN2018B" = mutate_dna(
    ancestor_2, 10
  ),
  
  "Rs7327" = mutate_dna(ancestor_2, 3),
  "Rs9401" = mutate_dna(ancestor_2, 5),
  "WIV16"  = mutate_dna(ancestor_2, 8),
  "Rs4231" = mutate_dna(ancestor_2, 10)
  
))

sequences_mutated


#-----------------
# alignment 
library(msa)

alignment <- msa(
  sequences_mutated,
  method = "ClustalW"
)

print(alignment)

# generate Neighbor-Joining tree
library(ape)

dna <- as.DNAbin(alignment)
print(dna)


distances <- dist.dna(
  dna,
  model = "JC"
)

print(as.matrix(distances))

tree_nj <- nj(distances)

# Bootstrap Neighbor-Joining tree (100 replicates)
bsnj <- boot.phylo(
  phy = tree_nj,
  x = dna,
  FUN = function(x) {
    nj(dist.dna(x, model = "JC"))
  },
  B = 100,
  trees = TRUE
)
# Plot Neighbor-Joining tree with bootstrap support
plot(
  tree_nj,
  main = "Neighbor-Joining phylogeny",
  edge.width = 2,
  cex = 1.1
)
# Add bootstrap support
nodelabels(
  text = bsnj$BP,
  node = (Ntip(tree_nj) + 1):(Ntip(tree_nj) + tree_nj$Nnode),
  frame = "none",
  cex = 0.9,
  adj = c(-0.2, -0.2),
  col = "red"
)
add.scale.bar(0.5, 1.2)

# Create a Maximum-Likelihood tree
library(phangorn)

# Convert DNAbin to phyDat (phylip format)
phy_data <- phyDat(dna, type = "DNA")

# Initialize ML using NJ topology
fit <- pml(tree_nj, data = phy_data)

# Optimize likelihood
fit_ml <- optim.pml(
  fit,
  model = "GTR",
  optNni = TRUE,
  optBf = FALSE,
  optQ = FALSE,
  optGamma = FALSE,
  optInv = FALSE
)

# Extract tree
tree_ml <- fit_ml$tree


# Bootstrap ML tree (100 replicates)
bsml <- boot.phylo(
  phy = tree_ml,
  x = dna,
  FUN = function(x) {
    nj(dist.dna(x, model = "K80"))
  },
  B = 100,
  trees = TRUE
)
# Plot ML tree with bootstrap support
plot(
  tree_ml,
  main = "Maximum-Likelihood phylogeny",
  edge.width = 2,
  cex = 1.1
)
# Add bootstrap support
nodelabels(
  text = bsml$BP,
  node = (Ntip(tree_ml) + 1):(Ntip(tree_ml) + tree_ml$Nnode),
  frame = "none",
  cex = 0.9,
  adj = c(-0.2, -0.2),
  col = "red"
)
add.scale.bar(0.5, 1.2)

##-------------------- compare NJ and ML trees
par(mfrow = c(1, 2))
# Bootstrap Neighbor-Joining tree (100 replicates)
bs <- boot.phylo(
  phy = tree_nj,
  x = dna,
  FUN = function(x) {
    nj(dist.dna(x, model = "JC"))
  },
  B = 100,
  trees = TRUE
)
# Plot Neighbor-Joining tree with bootstrap support
plot(
  tree_nj,
  main = "Neighbor-Joining phylogeny",
  edge.width = 2,
  cex = 1.1
)
# Add bootstrap support
nodelabels(
  text = bsnj$BP,
  node = (Ntip(tree_nj) + 1):(Ntip(tree_nj) + tree_nj$Nnode),
  frame = "none",
  cex = 0.9,
  adj = c(-0.2, -0.2),
  col = "red"
)
##################
# Plot ML tree with bootstrap support
plot(
  tree_ml,
  main = "Maximum-Likelihood phylogeny",
  edge.width = 2,
  cex = 1.1
)
# Add bootstrap support
nodelabels(
  text = bsml$BP,
  node = (Ntip(tree_ml) + 1):(Ntip(tree_ml) + tree_ml$Nnode),
  frame = "none",
  cex = 0.9,
  adj = c(-0.2, -0.2),
  col = "red"
)
par(mfrow = c(1, 1))

# Viral ExamR : midterm 2 
#
# loading -> host -> sequences -> alignment -> similarity -> tree
#
# Open this script in RStudio and run section by section.
# It sources the editable R/ files, so reinstalling after each edit is unnecessary.

# 1. Load the current package functions ---------------------------------------

# Install pak (only needed once)
install.packages("pak")

# Install the package from GitHub
pak::pak("brisenolab/PPEM454")

# Load the package in each new R session
library(ViralExamR)

# package_dir <- "/Users/ria/git_local/PPEM454/"
# for (script in c("choose_host.R", "get_sequences.R", "align_sequences.R",
#                  "sequence_similarity.R", "kmers.R")) {
#   source(file.path(package_dir, "R", script))
# }

# Alternatively, after installing the package, use library(ViralExamR).

# Install dependencies once if needed (these lines are deliberately commented):
install.packages(c("BiocManager", "ape"))
BiocManager::install(c("Biostrings", "msa"))
install.packages("phangorn") # only for the optional ML section

# load libraries
library("Biostrings")
library("msa")
library("ape")
library("phangorn")

# 2. Choose a host and load up to eight existing sequences --------------------
host <- choose_host()      # interactive menu alternative

# load sequences from your selected host
sequences <- get_sequences(host)
# print the name of your sequences and lenght
print(data.frame(sequence = names(sequences), length = nchar(sequences)))

# 3. Align and inspect -------------------------------------------------------
# align your sequence, choose one of these three methods: 
## "ClustalW", "ClustalOmega", "Muscle"
alignment <- align_sequences(sequences, method = "ClustalW")
print(alignment)

# 4. Answer the similarity questions -----------------------------------------
# Run this section: from line 44 to 58
# Convert aligned sequences to ape's DNAbin format
alignment_matrix <- do.call(
  rbind,
  strsplit(chartr("U", "T", toupper(alignment)), "", fixed = TRUE)
)
dna <- ape::as.DNAbin(alignment_matrix)
# Raw distance is the proportion of differing bases.
raw_distances <- ape::dist.dna(
  dna,
  model = "raw",
  pairwise.deletion = TRUE
)
# Convert it to percentage identity; exclude gaps/ambiguities pairwise.
identity_matrix <- 100 * (1 - as.matrix(raw_distances))
diag(identity_matrix) <- 100

# Here are the pairwise sequence identity results
print(round(identity_matrix, 2))
# or 
View(identity_matrix)

##@@ QUESTION
# The overall is the mean of pairwise identities. 
# For each pair calculated, gaps and ambiguous symbols are excluded. 

# Q1. What is the overall pairwise identity? use the function mean()
# Write your result:_______________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
# Q2. Which pair is most and the least similar, and what percentage they share?
# Write your result:_______________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________


# 5. Build and plot a Neighbor-Joining tree -----------------------------------
# Select a model: "RAW" "JC69" "K80" "F81" "K81" "F84" "T92""TN93" "GG95"
distances <- ape::dist.dna(dna, model = "JC69", pairwise.deletion = FALSE)
# review the evolutionary distances
print(round(as.matrix(distances), 4))

# Calculate Neighbor-Joining tree
tree_nj <- ape::nj(distances)
# Plot Neighbor-Joining tree
plot(tree_nj, type = "phylogram", cex = 0.9,
     main = paste("Neighbor-Joining tree:", host, "virus"))
ape::add.scale.bar()

# 6. Maximum-Likelihood tree
phy_data <- phangorn::phyDat(dna, type = "DNA")
# Calculate Neighbor-Joining tree
fit <- phangorn::pml(tree_nj, data = phy_data)
fit_ml <- phangorn::optim.pml(fit, model = "GTR", optNni = TRUE,
                              optBf = FALSE, optQ = FALSE,
                              optGamma = FALSE, optInv = FALSE)
tree_ml <- fit_ml$tree
# Plot Maximum-Likelihood tree
plot(tree_ml, type = "phylogram", main = paste("Maximum-Likelihood tree:", host, "virus"))

##@@ QUESTIONS
# Q4. Do the most similar sequences cluster together? 
# Write your interpretation:_______________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
# __________________________________________________________________________
# Q5. Are the trees topologies from NJ and ML identical? Why?
# Write your interpretation:_______________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
# Q6. What happened if you change the evolutionary model in N-J and ML?
# Write your interpretation:_______________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
# Q7. Do you notice a change in the branch length? Why is that?
# Write your interpretation:_______________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

# 7. Bootstrapping
# We will perform a bootstrap test in the ML tree using 100 replicates
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
add.scale.bar(0.75, 1)

# We did not have a chance to review bootstrap values in a trees, check
# this manuscript: 
# Q8. What is bootstrap and what does it mean these values in your tree?
# Write your interpretation:_______________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
# Q9. What if your bootstrap values are low confidence?
# Write your interpretation:_______________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

# 8. k-mers. Describe the properties of k-mers.
k=3
# select the highest divergence sequences 
# go to line 63 and View() the identity_matrix
kmers_spp <- count_kmers(sequences[c(1,8)], k)

# 9. Comparing two different strains
# Plot frequency distributions
a <- kmers_spp[[1]]
b <- kmers_spp[[2]]

# Align k-mers and replace missing counts with zero
kmers <- sort(union(names(a), names(b)))
counts <- rbind(a[kmers], b[kmers])
counts[is.na(counts)] <- 0

# Plotting k-mer sizing between a pair of sequences
barplot(
  counts,
  beside = TRUE,
  names.arg = kmers,
  col = c("steelblue", "orange"),
  las = 2,
  cex.names = 0.6,
  ylab = "Count",
  legend.text = c(names(kmers_spp[1]), names(kmers_spp[2])),
  args.legend = list(cex = 0.7)
)

# calculating k-mer similarity using Jaccard index
similarity <- length(intersect(a, b)) / length(union(a, b))
similarity * 100

# Q10. Would kmer comparison between two viral isolates or more would reflect
# sequence differences?
# Write your interpretation:_______________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________




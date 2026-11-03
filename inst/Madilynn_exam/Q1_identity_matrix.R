# Q1 Matrix

identity_matrix <- matrix(c(
  100,    99.37, 99.10, 100,   98.83, 99.10, 97.84, 90.35,
  99.37, 100,   99.19, 99.37, 98.92, 99.19, 98.11, 90.35,
  99.10, 99.19, 100,   99.10, 99.73, 100,   98.02, 90.63,
  100,   99.37, 99.10, 100,   98.83, 99.10, 97.84, 90.35,
  98.83, 98.92, 99.73, 98.83, 100,   99.73, 97.75, 90.63,
  99.10, 99.19, 100,   99.10, 99.73, 100,   98.02, 90.63,
  97.84, 98.11, 98.02, 97.84, 97.75, 98.02, 100,   90.53,
  90.35, 90.35, 90.63, 90.35, 90.63, 90.63, 90.53, 100
), nrow = 8, byrow = TRUE)

accessions <- c(
  "KU687048", "KU687046", "KU687040", "KU687045",
  "KU687042", "KU687041", "KR350618", "AF416468"
)

rownames(identity_matrix) <- accessions
colnames(identity_matrix) <- accessions

print(identity_matrix)

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


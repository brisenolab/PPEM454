# Viral ExamR : midterm 2 part 2
#
# loading -> select data set -> crate incidence matrix -> diversity -> resultls -> plots 
#
# Open this script in RStudio and run section by section.
# It sources the editable R/ files, so reinstalling after each edit is unnecessary.

# Run once if needed:
install.packages(c("vegan", "ggplot2", "tidyr"))

library(ViralExamR)
library(vegan)
library(ggplot2)
library(tidyr)

# 1. Choose a seed number: more than 4 random numbers 
# Q1. Write down your seed number: _________________
seed =  999 # after equal, type those random numbers

# 2. Create incidence matrix for viruses (as rows), samples (as columns)
inc_matrix <- create_incidence_matrix(10, 8, seed = seed)
inc_matrix

# 3. Calculate the total abundance for each virus
rowSums(inc_matrix)
# Q2. What does total abundance by sample represent?
# Write your response:_____________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
# vegan expects rows = samples, columns = taxa.
community <- t(inc_matrix) # transpose matrix

# Let's calculate the alpha diversity of our sample sites
# 2. Calculate Richness and alpha diversities
Richness <- specnumber(community)
Richness
# Q3. What is the richness for each sample? What sample has the highest and 
# lowest richness? Overall what does this suggest?
# Write your response:_____________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

Shannon <- diversity(community, index = "shannon")
Shannon
# Q4. What is the Shannon diversity index for each sample? What sample has 
# the highest and lowest richness? Overall, what does this suggest?
# Write your response:_____________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

Simpson <-  diversity(community, index = "simpson")
Simpson
# Q5. What is the Simpson diversity index for each sample? What sample has 
# the highest and lowest richness? Overall, what does this suggest?
# Write your response:_____________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

Inverse_Simpson <-  diversity(community, index = "invsimpson")
Inverse_Simpson
# Q6. What is the Inverse Simpson diversity index for each sample? What sample 
# has the highest and lowest richness? Overall, what does this suggest?
# Write your response:_____________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

# Calculate Pielou's evenness
# Divide Shanon diversity by the Richness log
Evenness <- Shannon / log(Richness) # hint type the object for Shannon and the log for Richness
Evenness
# Q6. What is the Evenness index for each sample? What sample has the 
# highest and lowest evenness? Overall, what does this suggest?
# Write your response:_____________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

# Create a table with the Alpha diversity results
Alpha_diversity_results <- data.frame(
  Sample = rownames(community),
  Abundance = rowSums(community),
  # vegan package using diversity() function 
  Richness = Richness,
  Shannon = Shannon,
  Simpson = Simpson,
  Inverser_Simpson = Inverse_Simpson,
  Evenness = Evenness,
  row.names = NULL
)
print(Alpha_diversity_results)

# 3. Plot richness and diversity indices
boxplot(
  results$Richness,
  results$Shannon,
  results$Simpson,
  results$Inverse_Simpson,
  results$Evenness,
  names = c("Richness", "Shannon", "Simpson", "Inverse Simpson", "Evennes"),
  ylab = "Index value",
  main = "Viral alpha diversity",
  col = c("lightblue", "lightgreen", "salmon", "khaki", "gray")
)

# Q7. Write a figure caption suitable for a scientific paper.
# Describe what is shown.
# Write your response: _______________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________


# 4. Calculate Beta-Diversity using Bray-Curtis index
bray <- vegdist(community, method = "bray")
bray_matrix <- as.matrix(bray)

# Remove self-comparisons (distance to itself = 0)
diag(bray_matrix) <- NA
bray_matrix

# Q8. What is the Bray-Curtis dissimilarity index? What sample has the 
# highest and lowest values? Overall, what does this suggest?
# Write your response:_____________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

# Calculating the Principal Components Analysis
pcoa <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)

# Percentage of variation represented by each axis
percent <- 100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])

# Base R plot
{plot(pcoa$points,  pch = 19, 
     col = c("lightblue","lightgreen", "salmon", "khaki", "purple",  
             "orange", "blue4", "red3"),
     cex = 1.5, asp = 1,
     xlab = paste0("PCoA 1 (", round(percent[1], 1), "%)"),
     ylab = paste0("PCoA 2 (", round(percent[2], 1), "%)"),
     main = "PCoA of viral communities")
text(pcoa$points, labels = rownames(community), pos = 2,  cex = 0.9)
}

# Q9. Principal Coordinates Analysis (PCoA) is an ordination method that
# represents differences among samples along coordinate axes. 
# Describe the patterns in the figure. Which samples cluster together,
# and which appear most different? (click on Zoom to see figure)
# Write your response: ____________________________________________________
# ________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________


# 5. Plot relative abundance to explain the diversity patterns
relative_abundance <- decostand(community, method = "total")
# format data matrix 
composition <- data.frame(Sample = rownames(relative_abundance),
                          relative_abundance, check.names = FALSE )
# transform data matrix to a long table to be compatible with ggplot
composition_long <- pivot_longer(composition, cols = -Sample, 
                                 names_to = "Virus", values_to = "Proportion")
# plot the community compositio
ggplot(composition_long, aes(x = Sample, y = Proportion, fill = Virus))+
  geom_col(width = 0.7) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c("virus1" = "#4477AA", "virus2" = "#EEAA33",
                               "virus3" = "#228833","virus4" = "#AA3377",
                               "virus5"  = "#66CCEE", "virus6"  = "#CC6677",
                               "virus7"  = "#AA4499", "virus8"  = "#999933",
                               "virus9"  = "#44AA99", "virus10" = "#882255"))+
  labs(title="Viral community composition", x="Samples", y="Relative abundance", 
       fill="Virus")+
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

# Q10. Principal Coordinates Analysis (PCoA) is an ordination method that
# represents differences among samples along coordinate axes. 
# Describe the patterns in the figure. Which samples cluster together,
# and which appear most different? (click on Zoom to see figure)
# Write your response: ____________________________________________________
# ________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________


# Bonus Point. Calculate Beta-Diversity using a different index
# let's check the documentation, and select a different method and write it 
# down between the " " in line 193.
?vegdist
other_index <- vegdist(community, method = "chao")
other_index_matrix <- as.matrix(other_index)

# Remove self-comparisons (distance to itself = 0)
diag(other_index_matrix) <- NA
other_index_matrix

# Calculating the Principal Components Analysis
npcoa <- cmdscale(other_index, k = 2, eig = TRUE, add = TRUE)

# Percentage of variation represented by each axis
percent <- 100 * npcoa$eig / sum(npcoa$eig[npcoa$eig > 0])

# Base R plot
{plot(npcoa$points,  pch = 19, 
      col = c("lightblue","lightgreen", "salmon", "khaki", "purple",  
              "orange", "blue4", "red3"),
      cex = 1.5, asp = 1,
      xlab = paste0("PCoA 1 (", round(percent[1], 1), "%)"),
      ylab = paste0("PCoA 2 (", round(percent[2], 1), "%)"),
      main = "PCoA of viral communities")
  text(npcoa$points, labels = rownames(community), pos = 2,  cex = 0.9)
}

# Q10. Can you notice a difference with the new dissimilarity index?
# Yes or No?
# What did it change:_____________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________
#__________________________________________________________________________

# Note, you can click on the arrows to review previous plots.

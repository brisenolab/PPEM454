# Run once if needed:
install.packages(c("vegan", "ggplot2", "tidyr"))

library(vegan)
library(ggplot2)
library(tidyr)

# 1. Create incidence matrix for viruses (as rows), samples (as columns)
counts <- rbind(
  "Virus A" = c(40, 30,   0,  0,  0),
  "Virus B" = c( 0, 10,   0,  0,  0),
  "Virus C" = c( 0,  0, 100, 10, 10),
  "Virus D" = c( 0, 25,   0, 36, 32)
)
colnames(counts) <- paste0("S", 1:5)

# Check total abundance for each virus
rowSums(counts)

# vegan expects rows = samples, columns = taxa.
community <- t(counts) # transpose matrix

# 2. Calculate alpha diversity
results <- data.frame(
  Sample = rownames(community),
  Abundance = rowSums(community),
  # vegan package using diversity() function 
  Richness = specnumber(community),
  Shannon = diversity(community, index = "shannon"),
  Simpson = diversity(community, index = "simpson"),
  Inverse_Simpson = diversity(community, index = "invsimpson"),
  row.names = NULL
)
print(results)

# Calculate Pielou's evenness
# Divide Shanon diversity by Richness
results$Evenness <- NA
multiple_taxa <- results$Richness > 1
results$Evenness[multiple_taxa] <-
  results$Shannon[multiple_taxa] / log(results$Richness[multiple_taxa])

print(results, digits = 3)

# 3. Plot richness and diversity indices

boxplot(
  results$Richness,
  results$Shannon,
  results$Simpson,
  results$Inverse_Simpson,
  names = c("Richness", "Shannon", "Simpson", "Inverse Simpson"),
  ylab = "Index value",
  main = "Viral alpha diversity",
  col = c("lightblue", "lightgreen", "salmon", "khaki")
)

# 4. Bray-Curtis: differences in relative abundance
bray_matrix <- as.matrix(bray)

# Remove self-comparisons (distance to itself = 0)
diag(bray_matrix) <- NA

# Comparing between samples
boxplot(
  as.data.frame(bray_matrix),
  xlab = "Sample",
  ylab = "Bray-Curtis dissimilarity",
  main = "Beta diversity by sample",
  col = c("lightblue","lightgreen", "salmon", "khaki", "purple"),
  ylim = c(0, 1)
)

# Calculating the Principal Components Analysis
pcoa <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)

# Percentage of variation represented by each axis
percent <- 100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])

# Base R plot
plot(pcoa$points,  pch = 19, 
     col = c("lightblue","lightgreen", "salmon", "khaki", "purple"),
      cex = 1.5, asp = 1,
     xlab = paste0("PCoA 1 (", round(percent[1], 1), "%)"),
     ylab = paste0("PCoA 2 (", round(percent[2], 1), "%)"),
     main = "PCoA of viral communities")
text(pcoa$points, labels = rownames(community), pos = 2,  cex = 0.9)

# 6. Plot relative abundance to explain the diversity patterns
relative_abundance <- decostand(community, method = "total")

# format data matrix 
composition <- data.frame(Sample = rownames(relative_abundance),
                            relative_abundance, check.names = FALSE )
# transform data matrix
composition_long <- pivot_longer(composition, cols = -Sample, 
                                 names_to = "Virus", values_to = "Proportion")

ggplot(composition_long, aes(x = Sample, y = Proportion, fill = Virus))+
  geom_col(width = 0.7) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                      expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c("Virus A" = "#4477AA", "Virus B" = "#EEAA33",
                                "Virus C" = "#228833","Virus D" = "#AA3377"))+
  labs(title="Viral community composition", x="Sample", y="Relative abundance", 
       fill="Virus")+
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

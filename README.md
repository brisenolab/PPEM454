# Midterm II: Sequence diversity, phylogenies and diversity indices/

## This exam consists of two R scripts:

- **Part 1:** Sequence alignment and phylogenetic analysis.
- **Part 2:** Diversity indices.

Run both scripts and answer all **20 questions and the bonus question** based on your results.

You may write your answers directly in the R scripts as comments, after editing, save the script as Last_name_MidterII_date.R turn it into Canvas or send me an email to ria5282@psu.edu. 

If you opt to use a Word document, it must be **electronic Word with version history.** Use electronic word that records your editing history, and ensure that this history is accessible when you submit your work. Submissions without accessible version history in word will not be graded. Share the document to ria5282@psu.edu.

**Do not use AI tools to complete any part of this exam.**

```r
install.packages("devtools")
library("devtools")
# Installing hoboR through devtools
devtools::install_github("brisenolab/PPEM454")
library(PPEM454)

# After installing locally with R CMD INSTALL /path/to/ViralExamR:
library(PPEM454) # to load the commands for the exam
# PPEM454 functions:
choose_host.R()
get_sequences.R()
align_sequences.R()
sequence_similarity.R()
kmers.R()
incidence_matrix()

# Install these packages:
install.packages(c("BiocManager", "ape"))
BiocManager::install(c("Biostrings", "msa"))
install.packages(c("phangorn", "vegan", "ggplot2", "tidyr"))

# load libraries
library("Biostrings")
library("msa")
library("ape")
library(vegan)
library(ggplot2)
library(tidyr)
```



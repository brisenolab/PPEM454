# Midterm II: Sequence diversity, phylogenies and diversity indices/

## This exam consists of two R scripts:

- **Part 1:** Sequence alignment and phylogenetic analysis.
- **Part 2:** Diversity indices.

Run both scripts and answer all **20 questions and the bonus question** based on your results.

_You may submit your answers in either of the following formats:_

- **R scripts:** Write your answers directly in the scripts as comments. Save each script as `LastName_MidtermII_Part1_DATE.R` or `LastName_MidtermII_Part2_DATE.R`. Submit both files through Canvas or email them to ria5282@psu.edu.

- **Word document:** Copy the questions from the script in an online Microsoft Word document with version history enabled. Share the document with ria5282@psu.edu and ensure that its version history is accessible. 
Word submissions without accessible version history will not be graded.

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



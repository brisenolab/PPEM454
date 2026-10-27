# ViralExamR

Editable teaching package scaffold. Main functions are in `R/`.

```r
# After installing locally with R CMD INSTALL /path/to/ViralExamR:
library(ViralExamR)
host <- choose_host("human") # choose_host() prompts interactively
sequences <- get_sequences(host, n = 8)
# Optional: BiocManager::install(c("Biostrings", "msa"))
alignment <- align_sequences(sequences)
sequence_similarity(alignment)
```

Existing filenames identify the host mapping: human, animal, fungal -> fungus,
phage -> bacteria, and plant. No sequences were invented. This mapping follows
the supplied filenames; source provenance and biological suitability have not
been independently verified. Loading selects the first 2–8 available records.
Review accessions, homologous regions, lengths, and expected answers before use.
The animal file contains RNA U; some files have gaps and ambiguous bases.
Alignment converts U to T and strips pre-existing gaps on an in-memory copy.
Similarity excludes gaps and ambiguities pairwise; overall is the mean identity
of unique pairs, in percent. These definitions can be edited for your exam.

Original teaching scripts were moved unchanged to `inst/legacy-scripts/` because
they install dependencies and run analyses at top level. Adapt their functions
into `R/` as needed. Original FASTA files remain in `data/`, with installable
copies in `inst/extdata/`. Existing `outputs/` is preserved but excluded from builds.

Before distribution: replace the placeholder author email in DESCRIPTION,
confirm GPL-3 is appropriate for your material, and review sequence provenance.
The minimal NAMESPACE is maintained manually. Roxygen comments are supplied for
future documentation generation. Base-R checks live in `tests/smoke.R`.

## Complete editable example

Open `inst/examples/host_sequence_workflow.R` and run it section by section.
It sources the current R/ scripts, loads supplied sequences, aligns them,
summarizes identity, builds an unrooted Neighbor-Joining tree using JC distances,
and saves alignment, tables, Newick tree, and PDF into a new outputs/ folder.
An optional Maximum-Likelihood section uses phangorn with the K80 model.

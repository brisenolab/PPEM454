library(ViralExamR)
stopifnot(identical(choose_host("FUNGUS"), "fungus"))
for (host in c("human", "animal", "fungus", "bacteria", "plant")) {
  x <- get_sequences(host)
  stopifnot(length(x) >= 2L, length(x) <= 8L, all(nzchar(x)))
}
x <- sequence_similarity(c(a = "ACGT", b = "ACGA", c = "ACGG"))
stopifnot(x$overall == 75, x$min == 75, x$max == 75)
x <- sequence_similarity(c("AU-N", "ATGN"))
stopifnot(x$overall == 100)
x <- sequence_similarity(c("NN", "--"))
stopifnot(is.na(x$overall))
stopifnot(inherits(try(choose_host("invalid"), silent = TRUE), "try-error"))
stopifnot(inherits(try(sequence_similarity(c("A", "AA")), silent = TRUE), "try-error"))

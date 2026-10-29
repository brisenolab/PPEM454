#' Create an incidence matrix
#' @param set.seed to select a fix seed number for the incidence matrix
#' @param nrow to select the number of rows
#' @param ncol to select the number of columns
#' @return An incidence matrices with viruses and hosts
#' @export
#' 
create_incidence_matrix <- function(nrow = 10, ncol = 8, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  values <- rnorm(nrow * ncol, mean = 40, sd = 60) |>
    round(0)
  
  values[values < 0] <- 0
  
  mat <- matrix(values, nrow = nrow, ncol = ncol)
  rownames(mat) <- paste0("virus", seq_len(nrow))
  colnames(mat) <- paste0("sample", LETTERS[seq_len(ncol)])
  
  return(mat)
}
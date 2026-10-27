#' Choose a host for the sequence exercise
#' @param host One of human, animal, fungus, bacteria, or plant. If NULL,
#'   prompt in an interactive R session.
#' @return A host name as a character scalar.
#' @export
choose_host <- function(host = NULL) {
  hosts <- c("human", "animal", "fungus", "bacteria", "plant")
  if (is.null(host)) {
    if (!interactive()) stop("Supply a host in non-interactive sessions.")
    choice <- utils::menu(hosts, title = "Choose your host")
    if (choice == 0L) stop("Host selection cancelled.")
    return(hosts[[choice]])
  }
  if (!is.character(host) || length(host) != 1L || is.na(host) ||
      !tolower(host) %in% hosts) stop("Choose human, animal, fungus, bacteria, or plant.")
  tolower(host)
}

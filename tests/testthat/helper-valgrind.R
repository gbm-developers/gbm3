skip_on_valgrind <- function() {
  use_valgrind <- tolower(Sys.getenv("_R_CHECK_USE_VALGRIND_", "false"))
  valgrind_env <- nzchar(Sys.getenv("R_VALGRIND_OPTS")) ||
    nzchar(Sys.getenv("VALGRIND_OPTS"))

  if (use_valgrind %in% c("true", "yes", "1") || valgrind_env) {
    testthat::skip("Skipping long-running integration tests under valgrind")
  }
}

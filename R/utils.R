# Utility functions for the package

#' Get default cache directory
#' @keywords internal
get_default_cache_dir <- function() {
  tools::R_user_dir("genokits", which = "cache")
}

#' Search NCBI for accessions
#' @keywords internal
search_ncbi_accessions <- function(search_term, api_key = NULL) {
  # This is a simplified version - you might want to use rentrez for better search
  message("Searching NCBI for: ", search_term)

  # For now, return empty - you'd implement proper NCBI search here
  # This would use entrez_search from rentrez package
  character(0)
}

#' Get local accessions from cache
#' @keywords internal
get_local_accessions <- function(cache_dir) {
  metadata_dir <- file.path(cache_dir, "metadata")

  if (!dir.exists(metadata_dir)) {
    return(character(0))
  }

  files <- list.files(metadata_dir, pattern = "\\.csv$")
  gsub("\\.csv$", "", files)
}

#' Set NCBI API key
#' @param api_key Your NCBI API key
#' @export
set_ncbi_api_key <- function(api_key) {
  Sys.setenv(NCBI_API_KEY = api_key)
  message("NCBI API key set successfully")
}

#' Get package version
#' @export
genokits_version <- function() {
  utils::packageVersion("genokits")
}

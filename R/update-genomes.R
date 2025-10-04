#' Incremental Update Genomes
#'
#' Check for and download new genomes from NCBI. Default for ASFV
#'
#' @param search_term NCBI search term for genomes
#' @param cache_dir Directory to cache downloaded data
#'
#' @return An update system object
#' @importFrom utils head
#' @export
#'
#' @examples
#' \dontrun{
#' update_genomes <- update_genomes('"African swine fever virus"[Organism]')
#' update_info <- update_genomes$check_updates()
#' if (update_info$needs_update) {
#'   update_genomes$perform_update()
#' }
#' }
update_genomes <- function(search_term = NULL, cache_dir = NULL) {

  if (is.null(search_term)) {
    search_term <- '"African swine fever virus"[Organism] AND "complete genome"[Title] AND 160000:300000[SLEN]'
  }

  if (is.null(cache_dir)) {
    cache_dir <- get_default_cache_dir()
  }

  system <- list(
    search_term = search_term,
    cache_dir = cache_dir,

    check_updates = function(api_key = NULL) {
      message("Checking for updates...")

      # Get current NCBI accessions
      current_accessions <- search_ncbi_accessions(self$search_term, api_key)

      # Get local accessions
      local_accessions <- get_local_accessions(self$cache_dir)

      # Find new accessions
      new_accessions <- setdiff(current_accessions, local_accessions)

      update_info <- list(
        local_count = length(local_accessions),
        current_count = length(current_accessions),
        new_accessions = new_accessions,
        needs_update = length(new_accessions) > 0
      )

      class(update_info) <- "update_info"
      return(update_info)
    },

    perform_update = function(parallel = FALSE, api_key = NULL) {
      update_info <- self$check_updates(api_key)

      if (!update_info$needs_update) {
        message("Database is up to date!")
        return(update_info)
      }

      message("Found ", length(update_info$new_accessions), " new genomes")

      # Show examples
      if (length(update_info$new_accessions) > 0) {
        message("New examples: ",
                paste(head(update_info$new_accessions, 3), collapse = ", "))
        if (length(update_info$new_accessions) > 3) {
          message("   ... and ", length(update_info$new_accessions) - 3, " more")
        }
      }

      # Ask for confirmation in interactive sessions
      if (interactive()) {
        response <- readline(prompt = "Download new genomes? (y/n): ")
        if (!tolower(response) %in% c("y", "yes")) {
          message("! Update cancelled by user")
          return(update_info)
        }
      }

      # Perform download
      message("... Starting download...")
      processor <- batch_processor(update_info$new_accessions, self$cache_dir)
      results <- processor$process_all(parallel = parallel, api_key = api_key)

      # Update results
      update_info$download_results <- results
      update_info$success_count <- results$success_count
      update_info$failed_count <- results$failed_count

      message("Update completed!")
      message("   Success: ", results$success_count, " genomes")
      message("   Failed:  ", results$failed_count, " genomes")
      message("   Metadata:", results$total_metadata_rows, " rows")
      message("   Features:", results$total_features_rows, " rows")

      return(update_info)
    },

    get_stats = function() {
      metadata_files <- list.files(
        file.path(self$cache_dir, "metadata"),
        pattern = "\\.csv$",
        full.names = TRUE
      )

      feature_files <- list.files(
        file.path(self$cache_dir, "features"),
        pattern = "_features\\.csv$",
        full.names = TRUE
      )

      sequence_files <- list.files(
        file.path(self$cache_dir, "sequences"),
        pattern = "\\.fasta$",
        full.names = TRUE
      )

      list(
        total_genomes = length(metadata_files),
        complete_genomes = length(intersect(
          gsub("\\.csv$", "", basename(metadata_files)),
          gsub("_features\\.csv$", "", basename(feature_files))
        )),
        metadata_files = length(metadata_files),
        feature_files = length(feature_files),
        sequence_files = length(sequence_files),
        total_size = sum(
          file.info(c(metadata_files, feature_files, sequence_files))$size,
          na.rm = TRUE
        )
      )
    }
  )

  system$self <- system
  class(system) <- "update_genomes"

  return(system)
}

#' @export
print.update_info <- function(x, ...) {
  cat("Update Information\n")
  cat("=================\n")
  cat("Local genomes:   ", x$local_count, "\n")
  cat("Current on NCBI: ", x$current_count, "\n")
  cat("New genomes:     ", length(x$new_accessions), "\n")
  cat("Needs update:    ", ifelse(x$needs_update, "Yes", "No"), "\n")

  if (x$needs_update) {
    cat("\nNew accessions:\n")
    cat(paste(" ", head(x$new_accessions, 5), collapse = "\n"))
    if (length(x$new_accessions) > 5) {
      cat("\n  ... and", length(x$new_accessions) - 5, "more")
    }
  }
}

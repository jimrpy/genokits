#' Genome Data Management System
#'
#' Create a genome data manager for handling metadata, sequences, and feature tables.
#'
#' @param accession Genome accession number
#' @param cache_dir Directory to cache downloaded data
#'
#' @return A genome data manager object with methods for data access
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' \dontrun{
#' manager <- genome_manager("NC_044967.1")
#' data <- manager$get_all_data()
#' }
genome_manager <- function(accession, cache_dir = NULL) {

  # Set default cache directory
  if (is.null(cache_dir)) {
    cache_dir <- get_default_cache_dir()
  }

  # Create the manager object
  manager <- list(
    accession = accession,
    cache_dir = cache_dir,

    # File path methods
    get_metadata_path = function() {
      file.path(self$cache_dir, "metadata", paste0(self$accession, ".csv"))
    },

    get_features_path = function() {
      file.path(self$cache_dir, "features", paste0(self$accession, "_features.csv"))
    },

    get_sequence_path = function() {
      file.path(self$cache_dir, "sequences", paste0(self$accession, ".fasta"))
    },

    # Data access methods
    get_all_data = function(force_refresh = FALSE) {
      metadata <- self$get_metadata(force_refresh)
      features <- self$get_features(force_refresh)

      list(
        accession = self$accession,
        metadata = metadata,
        features = features,
        files = list(
          metadata = self$get_metadata_path(),
          features = self$get_features_path(),
          sequence = self$get_sequence_path()
        )
      )
    },

    get_metadata = function(force_refresh = FALSE) {
      metadata_path <- self$get_metadata_path()
      sequence_path <- self$get_sequence_path()

      # Check cache
      if (file.exists(metadata_path) && file.exists(sequence_path) && !force_refresh) {
        message("Loading metadata from cache: ", self$accession)
        return(read.csv(metadata_path, stringsAsFactors = FALSE))
      }

      # Download new data
      message("Downloading metadata and sequence: ", self$accession)
      data <- fgenome_download(self$accession)

      # Ensure directories exist
      dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(sequence_path), recursive = TRUE, showWarnings = FALSE)

      # Save metadata (without sequence)
      if ("sequence" %in% names(data)) {
        metadata <- data[, !names(data) %in% "sequence", drop = FALSE]
        utils::write.csv(metadata, metadata_path, row.names = FALSE)

        # Save sequence separately
        if (nchar(data$sequence[1]) > 0) {
          fasta_content <- c(
            paste0(">", self$accession, " ", data$gene_desc[1]),
            data$sequence[1]
          )
          writeLines(fasta_content, sequence_path)
        }
      } else {
        utils::write.csv(data, metadata_path, row.names = FALSE)
      }

      message("Metadata saved to: ", metadata_path)
      message("Sequence saved to: ", sequence_path)

      return(data)
    },

    get_features = function(force_refresh = FALSE) {
      features_path <- self$get_features_path()

      # Check cache
      if (file.exists(features_path) && !force_refresh) {
        message("Loading features from cache: ", self$accession)
        return(utils::read.csv(features_path, stringsAsFactors = FALSE))
      }

      # Download new features
      message("Downloading feature table: ", self$accession)
      features <- ft_download(self$accession)

      if (!is.null(features)) {
        dir.create(dirname(features_path), recursive = TRUE, showWarnings = FALSE)
        utils::write.csv(features, features_path, row.names = FALSE)
        message("Features saved to: ", features_path)
      } else {
        message("Warning: Could not download features for ", self$accession)
      }

      return(features)
    },

    # Utility methods
    is_complete = function() {
      list(
        metadata = file.exists(self$get_metadata_path()),
        sequence = file.exists(self$get_sequence_path()),
        features = file.exists(self$get_features_path())
      )
    },

    clear_cache = function() {
      paths <- c(
        self$get_metadata_path(),
        self$get_features_path(),
        self$get_sequence_path()
      )
      unlink(paths[file.exists(paths)])
      message("Cache cleared for: ", self$accession)
    }
  )

  # Add self-reference for method chaining
  manager$self <- manager
  class(manager) <- "genome_manager"

  return(manager)
}

#' @export
print.genome_manager <- function(x, ...) {
  cat("Genome Manager for:", x$accession, "\n")
  cat("Cache directory:", x$cache_dir, "\n")

  completeness <- x$is_complete()
  cat("Data completeness:\n")
  cat("  Metadata:", ifelse(completeness$metadata, "ok", "not completeness"), "\n")
  cat("  Sequence:", ifelse(completeness$sequence, "ok", "not completeness"), "\n")
  cat("  Features:", ifelse(completeness$features, "ok", "not completeness"), "\n")
}

#' Genome Data Management System
#'
#' Create a genome data manager for handling metadata, sequences, and feature tables.
#'
#' @param accession Genome accession number
#' @param db_path Directory to store genome data. If NULL, uses the current genome database.
#'
#' @return A genome data manager object with methods for data access
#' @importFrom utils read.csv write.csv
#' @export
#'
#' @examples
#' \dontrun{
#' manager <- genome_manager("NC_044967.1")
#' data <- manager$get_all_data()
#' }
genome_manager <- function(accession, db_path = NULL) {

  # use the given db_path or current genebank folder
  if (is.null(db_path)) {
    db_path <- get_genome_db()
    message("Debug: Using get_genome_db() result: ", db_path)
  } else {
    # make sure the given path is absolute path and exist
    db_path <- normalizePath(db_path, mustWork = FALSE)
    if (!dir.exists(db_path)) {
      dir.create(db_path, recursive = TRUE, showWarnings = FALSE)
    }
    message("Debug: Using provided db_path: ", db_path)
  }

  # use local variables, avoid recurrent reference
  .accession <- accession
  .db_path <- db_path

  # define the function path
  get_metadata_path <- function() {
    file.path(.db_path, "metadata", paste0(.accession, ".csv"))
  }

  get_features_path <- function() {
    file.path(.db_path, "features", paste0(.accession, "_features.csv"))
  }

  get_sequence_path <- function() {
    file.path(.db_path, "sequences", paste0(.accession, ".fasta"))
  }

  # data access method
  get_metadata <- function(force_refresh = FALSE) {
    metadata_path <- get_metadata_path()
    sequence_path <- get_sequence_path()

    # Check cache
    if (file.exists(metadata_path) && file.exists(sequence_path) && !force_refresh) {
      message("Loading metadata from database: ", .accession)
      return(utils::read.csv(metadata_path, stringsAsFactors = FALSE))
    }

    # Download new data
    message("Downloading metadata and sequence: ", .accession)
    data <- fgenome_download(.accession)

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
          paste0(">", .accession, " ", data$gene_desc[1]),
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
  }

  get_features <- function(force_refresh = FALSE) {
    features_path <- get_features_path()

    # Check cache
    if (file.exists(features_path) && !force_refresh) {
      message("Loading features from database: ", .accession)
      return(utils::read.csv(features_path, stringsAsFactors = FALSE))
    }

    # Download new features
    message("Downloading feature table: ", .accession)
    features <- ft_download(.accession)

    if (!is.null(features)) {
      dir.create(dirname(features_path), recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(features, features_path, row.names = FALSE)
      message("Features saved to: ", features_path)
    } else {
      message("Warning: Could not download features for ", .accession)
    }

    return(features)
  }

  get_all_data <- function(force_refresh = FALSE) {
    metadata <- get_metadata(force_refresh)
    features <- get_features(force_refresh)

    list(
      accession = .accession,
      metadata = metadata,
      features = features,
      files = list(
        metadata = get_metadata_path(),
        features = get_features_path(),
        sequence = get_sequence_path()
      )
    )
  }

  is_complete <- function() {
    list(
      metadata = file.exists(get_metadata_path()),
      sequence = file.exists(get_sequence_path()),
      features = file.exists(get_features_path())
    )
  }

  clear_cache <- function() {
    paths <- c(
      get_metadata_path(),
      get_features_path(),
      get_sequence_path()
    )
    unlink(paths[file.exists(paths)])
    message("Data removed for: ", .accession)
  }

  # Return public interface
  structure(
    list(
      accession = .accession,
      db_path = .db_path,
      get_all_data = get_all_data,
      get_metadata = get_metadata,
      get_features = get_features,
      is_complete = is_complete,
      clear_cache = clear_cache
    ),
    class = "genome_manager"
  )
}

#' @export
print.genome_manager <- function(x, ...) {
  cat("Genome Manager for:", x$accession, "\n")
  cat("Database path:", x$db_path, "\n")

  completeness <- x$is_complete()
  cat("Data completeness:\n")
  cat("  Metadata:", ifelse(completeness$metadata, "Ok", "Not completeness"), "\n")
  cat("  Sequence:", ifelse(completeness$sequence, "Ok", "Not completeness"), "\n")
  cat("  Features:", ifelse(completeness$features, "Ok", "Not completeness"), "\n")
}

#' Genome Data Import and Export Functions
#'
#' Functions for sharing genome data with colleagues and importing exported data.
#'
#' @name data-import-export
NULL

#' Export genome data for sharing
#'
#' Export genome data to a portable format that can be shared with colleagues.
#' Creates a standardized archive containing all data for the specified genome.
#'
#' @param accession Genome accession number
#' @param export_dir Directory where the export package will be created.
#'   If NULL, creates in current working directory.
#' @param db_path Optional path to genome database. If NULL, uses current database.
#' @param include_sequence Whether to include sequence data in the export
#' @param include_features Whether to include feature data in the export
#' @return Path to the created export directory
#' @export
#'
#' @examples
#' \dontrun{
#' # Export all data for a genome
#' export_genome("ON400500.1", "./exports")
#'
#' # Export only metadata
#' export_genome("ON400500.1", "./exports",
#'               include_sequence = FALSE, include_features = FALSE)
#' }
export_genome <- function(accession, export_dir = NULL, db_path = NULL,
                          include_sequence = TRUE, include_features = TRUE) {

  if (is.null(db_path)) {
    db_path <- get_genome_db()
  }

  if (is.null(export_dir)) {
    export_dir <- file.path(getwd(), "genome_exports")
  }

  # Create export directory structure
  export_path <- file.path(export_dir, accession)
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

  # Define source paths
  src_metadata <- file.path(db_path, "metadata", paste0(accession, ".csv"))
  src_sequence <- file.path(db_path, "sequences", paste0(accession, ".fasta"))
  src_features <- file.path(db_path, "features", paste0(accession, "_features.csv"))

  # Copy files
  files_copied <- c()

  # Always copy metadata (required)
  if (file.exists(src_metadata)) {
    dest_metadata <- file.path(export_path, paste0(accession, "_metadata.csv"))
    file.copy(src_metadata, dest_metadata, overwrite = TRUE)
    files_copied <- c(files_copied, "metadata")
  }

  # Copy sequence if requested and exists
  if (include_sequence && file.exists(src_sequence)) {
    dest_sequence <- file.path(export_path, paste0(accession, ".fasta"))
    file.copy(src_sequence, dest_sequence, overwrite = TRUE)
    files_copied <- c(files_copied, "sequence")
  }

  # Copy features if requested and exists
  if (include_features && file.exists(src_features)) {
    dest_features <- file.path(export_path, paste0(accession, "_features.csv"))
    file.copy(src_features, dest_features, overwrite = TRUE)
    files_copied <- c(files_copied, "features")
  }

  # Create README file with import instructions
  readme_content <- c(
    paste("# Genome Data Export:", accession),
    paste("Exported on:", Sys.time()),
    "",
    "## Contents:",
    if ("metadata" %in% files_copied) "- Metadata (CSV format)",
    if ("sequence" %in% files_copied) "- Sequence (FASTA format)",
    if ("features" %in% files_copied) "- Feature table (CSV format)",
    "",
    "## Import Instructions:",
    "1. Ensure you have the genokits package installed",
    "2. Use the import_genome() function:",
    paste0('   import_genome("', export_path, '")'),
    "",
    "## Data Source:",
    paste("Database:", db_path)
  )

  writeLines(readme_content, file.path(export_path, "README.md"))

  message("Genome '", accession, "' exported to: ", export_path)
  message("Files included: ", paste(files_copied, collapse = ", "))

  return(invisible(export_path))
}

#' Import genome data from export with validation
#'
#' Import genome data that was previously exported using export_genome().
#' Includes comprehensive validation and user prompts for incomplete data.
#'
#' @param import_path Path to the exported genome directory or file
#' @param db_path Optional path to genome database. If NULL, uses current database.
#' @param overwrite Whether to overwrite existing data (default: FALSE)
#' @param interactive Whether to prompt user for decisions (default: TRUE)
#' @return List with import results and details
#' @export
#'
#' @examples
#' \dontrun{
#' # Interactive import
#' result <- import_genome("./genome_exports/ON400500.1")
#'
#' # Non-interactive import (for scripts)
#' result <- import_genome("./genome_exports/ON400500.1", interactive = FALSE)
#' }
import_genome <- function(import_path, db_path = NULL, overwrite = FALSE,
                          interactive = TRUE) {

  if (is.null(db_path)) {
    db_path <- get_genome_db()
  }

  # show the target database path when start
  if (interactive) {
    cat("Target database:", db_path, "\n")
  }

  # Validate import path
  if (!dir.exists(import_path)) {
    stop("Import path does not exist: ", import_path)
  }

  # Extract accession from directory name
  accession <- basename(import_path)

  # Check for existing genome in database
  existing_check <- .check_existing_genome(accession, db_path, overwrite, interactive)
  if (!existing_check$proceed) {
    return(list(
      success = FALSE,
      message = "Import cancelled by user",
      accession = accession,
      db_path = db_path
    ))
  }

  # Validate exported data completeness
  validation <- .validate_export_data(import_path, accession)

  # Prompt user for incomplete data (if interactive)
  if (interactive && !validation$complete) {
    proceed <- .prompt_incomplete_data(validation, accession, db_path)
    if (!proceed) {
      return(list(
        success = FALSE,
        message = "Import cancelled due to incomplete data",
        accession = accession,
        validation = validation,
        db_path = db_path
      ))
    }
  }

  # Perform the actual import
  import_result <- .perform_import(import_path, db_path, accession, validation)

  # Return detailed results
  return(list(
    success = import_result$success,
    message = ifelse(import_result$success,
                     paste("Successfully imported genome:", accession, "to", db_path),
                     paste("Failed to import genome:", accession, "to", db_path)),
    accession = accession,
    files_imported = import_result$files_imported,
    validation = validation,
    import_path = import_path,
    db_path = db_path
  ))
}

#' Create shareable genome package
#'
#' Create a compressed archive of genome data for easy sharing.
#'
#' @param accession Genome accession number
#' @param output_file Path to the output archive file (.zip or .tar.gz)
#' @param db_path Optional path to genome database. If NULL, uses current database.
#' @return Path to the created archive file
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a zip archive
#' share_genome("ON400500.1", "./ON400500.1.zip")
#' }
share_genome <- function(accession, output_file = NULL, db_path = NULL) {

  if (is.null(output_file)) {
    output_file <- file.path(getwd(), paste0(accession, ".zip"))
  }

  # First export to temporary directory
  temp_dir <- tempfile("genome_export_")
  export_path <- export_genome(accession, temp_dir, db_path)

  # Create archive
  if (grepl("\\.zip$", output_file)) {
    if (!requireNamespace("zip", quietly = TRUE)) {
      stop("The 'zip' package is required for creating zip archives. Please install it.")
    }
    zip::zip(output_file, files = export_path, recurse = TRUE)
  } else {
    # Use tar for other formats
    utils::tar(output_file, files = export_path, compression = "gzip")
  }

  # Clean up temporary directory
  unlink(temp_dir, recursive = TRUE)

  message("Shareable package created: ", output_file)
  return(invisible(output_file))
}

#' Validate genome export before sharing
#'
#' Check if an exported genome directory contains valid and complete data.
#'
#' @param export_path Path to the exported genome directory
#' @return List with validation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Validate before sharing
#' validation <- validate_export("./genome_exports/ON400500.1")
#' print(validation)
#' }
validate_export <- function(export_path) {
  if (!dir.exists(export_path)) {
    stop("Export path does not exist: ", export_path)
  }

  accession <- basename(export_path)
  validation <- .validate_export_data(export_path, accession)

  # Additional checks for sharing readiness
  sharing_ready <- validation$has_metadata &&
    validation$has_sequence &&
    validation$has_features &&
    length(validation$validation_notes) == 0

  return(list(
    sharing_ready = sharing_ready,
    accession = accession,
    validation = validation,
    recommendations = if (!sharing_ready) {
      .generate_sharing_recommendations(validation)
    } else {
      "Export is ready for sharing"
    }
  ))
}

#' Show current genome database information
#'
#' Display information about the current genome database location and contents.
#'
#' @return Invisibly returns the database path
#' @export
#'
#' @examples
#' \dontrun{
#' show_database_info()
#' }
show_database_info <- function() {
  db_path <- get_genome_db()

  cat("=== Genome Database Information ===\n")
  cat("Database path:", db_path, "\n")

  # Check if database exists and has content
  if (!dir.exists(db_path)) {
    cat("Status: Database directory does not exist\n")
    return(invisible(db_path))
  }

  # Count files in each sub_directory
  metadata_count <- length(list.files(file.path(db_path, "metadata"), pattern = "\\.csv$"))
  sequence_count <- length(list.files(file.path(db_path, "sequences"), pattern = "\\.fasta$"))
  features_count <- length(list.files(file.path(db_path, "features"), pattern = "\\.csv$"))

  cat("Genomes in database:\n")
  cat("  Metadata files:", metadata_count, "\n")
  cat("  Sequence files:", sequence_count, "\n")
  cat("  Feature files: ", features_count, "\n")

  # Show disk usage
  if (requireNamespace("fs", quietly = TRUE)) {
    db_size <- fs::dir_info(db_path, recurse = TRUE)$size
    total_size <- sum(db_size, na.rm = TRUE) / (1024^2) # Convert to MB
    cat("Total disk usage:", round(total_size, 2), "MB\n")
  } else {
    # Fallback method
    total_size <- sum(file.info(list.files(db_path, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE) / (1024^2)
    cat("Total disk usage:", round(total_size, 2), "MB\n")
  }

  cat("To change database location, use: set_genome_db(\"your/path/here\")\n")

  return(invisible(db_path))
}

# ===== Internal Helper Functions =====

#' Check if genome already exists in database
#' @keywords internal
.check_existing_genome <- function(accession, db_path, overwrite, interactive) {

  existing_files <- c()
  if (file.exists(file.path(db_path, "metadata", paste0(accession, ".csv")))) {
    existing_files <- c(existing_files, "metadata")
  }
  if (file.exists(file.path(db_path, "sequences", paste0(accession, ".fasta")))) {
    existing_files <- c(existing_files, "sequence")
  }
  if (file.exists(file.path(db_path, "features", paste0(accession, "_features.csv")))) {
    existing_files <- c(existing_files, "features")
  }

  if (length(existing_files) > 0) {
    if (overwrite) {
      if (interactive) {
        cat("Warning: Genome '", accession, "' already exists in database: ", db_path, "\n")
        cat("Existing files: ", paste(existing_files, collapse = ", "), "\n")
        cat("Overwrite mode is enabled. These files will be replaced.\n")
      }
      return(list(proceed = TRUE, existing_files = existing_files))
    } else {
      if (interactive) {
        cat("Genome '", accession, "' already exists in database: ", db_path, "\n")
        cat("Existing files: ", paste(existing_files, collapse = ", "), "\n")
        response <- readline(prompt = "Overwrite? (y/N): ")
        if (tolower(response) %in% c("y", "yes")) {
          return(list(proceed = TRUE, existing_files = existing_files))
        } else {
          return(list(proceed = FALSE, existing_files = existing_files))
        }
      } else {
        # Non-interactive mode, don't overwrite
        return(list(proceed = FALSE, existing_files = existing_files))
      }
    }
  }

  return(list(proceed = TRUE, existing_files = character()))
}

#' Validate exported data completeness
#' @keywords internal
.validate_export_data <- function(import_path, accession) {

  expected_files <- c(
    metadata = file.path(import_path, paste0(accession, "_metadata.csv")),
    sequence = file.path(import_path, paste0(accession, ".fasta")),
    features = file.path(import_path, paste0(accession, "_features.csv"))
  )

  # Check which files exist
  file_exists <- sapply(expected_files, file.exists)

  # Check file sizes (avoid empty files)
  file_sizes <- sapply(expected_files, function(f) {
    if (file.exists(f)) file.info(f)$size else 0
  })

  # Determine completeness
  has_metadata <- file_exists["metadata"] && file_sizes["metadata"] > 0
  has_sequence <- file_exists["sequence"] && file_sizes["sequence"] > 0
  has_features <- file_exists["features"] && file_sizes["features"] > 0

  complete <- has_metadata && has_sequence && has_features

  # Validate file formats (basic checks)
  validation_notes <- c()

  if (has_metadata) {
    # Check if metadata file is readable CSV
    tryCatch({
      utils::read.csv(expected_files["metadata"], nrows = 1)
    }, error = function(e) {
      validation_notes <<- c(validation_notes,
                             "Metadata file appears to be invalid CSV")
      has_metadata <<- FALSE
    })
  }

  if (has_sequence) {
    # Check if sequence file looks like FASTA
    tryCatch({
      lines <- readLines(expected_files["sequence"], n = 2)
      if (length(lines) < 1 || !startsWith(lines[1], ">")) {
        validation_notes <<- c(validation_notes,
                               "Sequence file doesn't appear to be valid FASTA")
        has_sequence <<- FALSE
      }
    }, error = function(e) {
      validation_notes <<- c(validation_notes,
                             "Cannot read sequence file")
      has_sequence <<- FALSE
    })
  }

  return(list(
    complete = complete,
    has_metadata = has_metadata,
    has_sequence = has_sequence,
    has_features = has_features,
    file_exists = file_exists,
    file_sizes = file_sizes,
    validation_notes = validation_notes,
    expected_files = expected_files
  ))
}

#' Prompt user about incomplete data
#' @keywords internal
.prompt_incomplete_data <- function(validation, accession, db_path) {

  cat("\n=== Genome Import Validation ===\n")
  cat("Genome:", accession, "\n")
  cat("Target database:", db_path, "\n")
  cat("Data completeness check:\n")
  cat("  Metadata:", ifelse(validation$has_metadata, "Validated", "Invalidated"), "\n")
  cat("  Sequence:", ifelse(validation$has_sequence, "Validated", "Invalidated"), "\n")
  cat("  Features:", ifelse(validation$has_features, "Validated", "Invalidated"), "\n")

  if (length(validation$validation_notes) > 0) {
    cat("\nValidation notes:\n")
    for (note in validation$validation_notes) {
      cat("  -", note, "\n")
    }
  }

  if (!validation$has_metadata) {
    cat("\n WARNING: Metadata file is missing or invalid.\n")
    cat("Metadata is required for genome import.\n")
    return(FALSE)
  }

  missing_parts <- c()
  if (!validation$has_sequence) missing_parts <- c(missing_parts, "sequence")
  if (!validation$has_features) missing_parts <- c(missing_parts, "features")

  if (length(missing_parts) > 0) {
    cat("\n NOTE: The following data is missing:\n")
    cat("  -", paste(missing_parts, collapse = ", "), "\n")
    cat("You can still import the available data to", db_path, "\n")
    cat("and download missing parts later using genome_manager().\n")

    response <- readline(prompt = "Proceed with partial import? (y/N): ")
    return(tolower(response) %in% c("y", "yes"))
  }

  return(TRUE)
}

#' Perform the actual import operation
#' @keywords internal
.perform_import <- function(import_path, db_path, accession, validation) {

  # Ensure destination directories exist
  dir.create(file.path(db_path, "metadata"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(db_path, "sequences"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(db_path, "features"), recursive = TRUE, showWarnings = FALSE)

  files_imported <- c()

  # Import metadata (required)
  if (validation$has_metadata) {
    src <- validation$expected_files["metadata"]
    dest <- file.path(db_path, "metadata", paste0(accession, ".csv"))
    if (file.copy(src, dest, overwrite = TRUE)) {
      files_imported <- c(files_imported, "metadata")
    }
  }

  # Import sequence (optional)
  if (validation$has_sequence) {
    src <- validation$expected_files["sequence"]
    dest <- file.path(db_path, "sequences", paste0(accession, ".fasta"))
    if (file.copy(src, dest, overwrite = TRUE)) {
      files_imported <- c(files_imported, "sequence")
    }
  }

  # Import features (optional)
  if (validation$has_features) {
    src <- validation$expected_files["features"]
    dest <- file.path(db_path, "features", paste0(accession, "_features.csv"))
    if (file.copy(src, dest, overwrite = TRUE)) {
      files_imported <- c(files_imported, "features")
    }
  }

  success <- length(files_imported) > 0 && "metadata" %in% files_imported

  if (success) {
    cat("Importing to:", db_path, "\n")
  }

  return(list(
    success = success,
    files_imported = files_imported
  ))
}

#' Generate recommendations for improving export quality
#' @keywords internal
.generate_sharing_recommendations <- function(validation) {
  recommendations <- c()

  if (!validation$has_metadata) {
    recommendations <- c(recommendations,
                         "Metadata file is missing. This is required for import.")
  }

  if (!validation$has_sequence) {
    recommendations <- c(recommendations,
                         "Sequence file is missing. Consider including it for complete data.")
  }

  if (!validation$has_features) {
    recommendations <- c(recommendations,
                         "Features file is missing. Consider including it for complete data.")
  }

  if (length(validation$validation_notes) > 0) {
    recommendations <- c(recommendations,
                         "Fix validation issues before sharing:",
                         validation$validation_notes)
  }

  return(recommendations)
}



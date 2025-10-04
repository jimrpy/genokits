#' Genome Database Management
#'
#' Functions for managing genome database locations and contents.
#'
#' @name database-management
NULL

#' Set genome database directory
#'
#' Set the directory where genome data will be stored. This allows users
#' to specify a persistent location for their genome databases.
#'
#' @param path Path to the genome database directory. If NULL, uses the
#'   default location in the R user data directory.
#' @param create Logical. Whether to create the directory if it doesn't exist.
#' @return The path to the genome database directory (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' # Set to a custom directory
#' set_genome_db("~/my_genome_database")
#'
#' # Set to a project-specific directory
#' set_genome_db("./genome_data")
#'
#' # Reset to default location
#' set_genome_db()
#' }
set_genome_db <- function(path = NULL, create = TRUE) {
  if (is.null(path)) {
    # Reset to default
    options(genokits.genome_db = NULL)
    message("Genome database reset to default location")
    db_path <- get_default_genome_db()
  } else {
    # Validate and set custom path
    path <- normalizePath(path, mustWork = FALSE)
    if (create && !dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(path)) {
      stop("Directory does not exist: ", path)
    }
    options(genokits.genome_db = path)
    message("Genome database set to: ", path)
    db_path <- path
  }

  return(invisible(db_path))
}

#' Get current genome database path
#'
#' @return Path to the current genome database directory
#' @export
#'
#' @examples
#' \dontrun{
#' get_genome_db()
#' }
get_genome_db <- function() {
  custom_path <- getOption("genokits.genome_db")
  if (!is.null(custom_path)) {
    return(custom_path)
  }
  return(get_default_genome_db())
}

#' Get default genome database path
#'
#' @return Path to the default genome database directory
#' @keywords internal
get_default_genome_db <- function() {
  # Use R user data directory by default
  db_path <- tools::R_user_dir("genokits", which = "data")
  dir.create(db_path, recursive = TRUE, showWarnings = FALSE)
  return(db_path)
}

#' List available genomes in database
#'
#' @param db_path Optional path to genome database. If NULL, uses current database.
#' @return Data frame with information about available genomes
#' @export
#'
#' @examples
#' \dontrun{
#' list_genomes()
#' list_genomes("~/my_genomes")
#' }
list_genomes <- function(db_path = NULL) {
  if (is.null(db_path)) {
    db_path <- get_genome_db()
  }

  # confirm the directory exist
  dirs <- c("metadata", "sequences", "features")
  for (dir in dirs) {
    dir_path <- file.path(db_path, dir)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }

  # Look for metadata files to identify available genomes
  metadata_dir <- file.path(db_path, "metadata")
  if (!dir.exists(metadata_dir)) {
    message("No genomes found in database: ", db_path)
    return(data.frame(
      accession = character(),
      has_metadata = logical(),
      has_sequence = logical(),
      has_features = logical(),
      stringsAsFactors = FALSE
    ))
  }

  metadata_files <- list.files(metadata_dir, pattern = "\\.csv$", full.names = FALSE)
  accessions <- gsub("\\.csv$", "", metadata_files)

  if (length(accessions) == 0) {
    message("No genomes found in database: ", db_path)
    return(data.frame(
      accession = character(),
      has_metadata = logical(),
      has_sequence = logical(),
      has_features = logical(),
      stringsAsFactors = FALSE
    ))
  }

  # Check what data is available for each genome
  result <- lapply(accessions, function(acc) {
    list(
      accession = acc,
      has_metadata = file.exists(file.path(db_path, "metadata", paste0(acc, ".csv"))),
      has_sequence = file.exists(file.path(db_path, "sequences", paste0(acc, ".fasta"))),
      has_features = file.exists(file.path(db_path, "features", paste0(acc, "_features.csv")))
    )
  })

  df <- do.call(rbind, lapply(result, function(x) {
    data.frame(
      accession = x$accession,
      has_metadata = x$has_metadata,
      has_sequence = x$has_sequence,
      has_features = x$has_features,
      stringsAsFactors = FALSE
    )
  }))

  # sort by accession id
  df[order(df$accession), ]
}

#' Check if genome exists in database
#'
#' @param accession Genome accession number
#' @param db_path Optional path to genome database. If NULL, uses current database.
#' @return Logical indicating whether the genome exists in the database
#' @export
#'
#' @examples
#' \dontrun{
#' has_genome("ON400500.1")
#' }
has_genome <- function(accession, db_path = NULL) {
  if (is.null(db_path)) {
    db_path <- get_genome_db()
  }

  metadata_path <- file.path(db_path, "metadata", paste0(accession, ".csv"))
  return(file.exists(metadata_path))
}

#' Remove genome from database
#'
#' @param accession Genome accession number
#' @param db_path Optional path to genome database. If NULL, uses current database.
#' @return Logical indicating success
#' @export
#'
#' @examples
#' \dontrun{
#' remove_genome("ON400500.1")
#' }
remove_genome <- function(accession, db_path = NULL) {
  if (is.null(db_path)) {
    db_path <- get_genome_db()
  }

  paths <- c(
    file.path(db_path, "metadata", paste0(accession, ".csv")),
    file.path(db_path, "sequences", paste0(accession, ".fasta")),
    file.path(db_path, "features", paste0(accession, "_features.csv"))
  )

  existing_paths <- paths[file.exists(paths)]
  if (length(existing_paths) == 0) {
    message("No data found for genome: ", accession)
    return(FALSE)
  }

  # delete confirm
  cat("These files will be deleted:\n")
  for (path in existing_paths) {
    cat("  -", path, "\n")
  }

  response <- readline(prompt = "Delete Confirmed? (y/N): ")
  if (tolower(response) %in% c("y", "yes")) {
    unlink(existing_paths)
    message("Removed data for genome: ", accession)
    return(TRUE)
  } else {
    message("Delete Canceled")
    return(FALSE)
  }
}

#' Get database statistics
#'
#' @param db_path Optional path to genome database. If NULL, uses current database.
#' @return List with database statistics
#' @export
#'
#' @examples
#' \dontrun{
#' get_db_stats()
#' }
get_db_stats <- function(db_path = NULL) {
  if (is.null(db_path)) {
    db_path <- get_genome_db()
  }

  genomes <- list_genomes(db_path)

  list(
    total_genomes = nrow(genomes),
    genomes_with_metadata = sum(genomes$has_metadata),
    genomes_with_sequence = sum(genomes$has_sequence),
    genomes_with_features = sum(genomes$has_features),
    complete_genomes = sum(genomes$has_metadata & genomes$has_sequence & genomes$has_features),
    db_path = db_path
  )
}

#' Clean empty directories in database
#'
#' Remove empty subdirectories from the genome database
#'
#' @param db_path Optional path to genome database. If NULL, uses current database.
#' @return Number of directories removed
#' @export
#'
#' @examples
#' \dontrun{
#' clean_empty_dirs()
#' }
clean_empty_dirs <- function(db_path = NULL) {
  if (is.null(db_path)) {
    db_path <- get_genome_db()
  }

  dirs_to_check <- c("metadata", "sequences", "features")
  removed_count <- 0

  for (dir_name in dirs_to_check) {
    dir_path <- file.path(db_path, dir_name)
    if (dir.exists(dir_path)) {
      files <- list.files(dir_path, full.names = FALSE)
      if (length(files) == 0) {
        unlink(dir_path, recursive = TRUE)
        removed_count <- removed_count + 1
        message("Removed empty directory: ", dir_path)
      }
    }
  }

  # check the main fold is NULL
  files_in_db <- list.files(db_path, full.names = FALSE)
  if (length(files_in_db) == 0) {
    unlink(db_path, recursive = TRUE)
    message("Removed empty database directory: ", db_path)
    removed_count <- removed_count + 1
  }

  message("Removed ", removed_count, " empty directories")
  return(invisible(removed_count))
}

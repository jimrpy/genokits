#' Batch Genome Processor
#'
#' Process multiple genomes in batch with progress tracking and error handling.
#'
#' @param accessions Vector of genome accession numbers
#' @param cache_dir Directory to cache downloaded data
#'
#' @return A batch processor object
#' @export
#'
#' @examples
#' \dontrun{
#' processor <- batch_processor(c("NC_044967.1", "NC_044968.1"))
#' results <- processor$process_all()
#' }
batch_processor <- function(accessions, cache_dir = NULL) {

  if (is.null(cache_dir)) {
    cache_dir <- get_default_cache_dir()
  }

  processor <- list(
    accessions = accessions,
    cache_dir = cache_dir,

    process_single = function(accession, index, total, api_key = NULL) {
      message(sprintf("[%d/%d] Processing: %s", index, total, accession))

      start_time <- Sys.time()

      result <- tryCatch({
        manager <- genome_manager(accession, self$cache_dir)

        # Download both metadata and features
        metadata <- manager$get_metadata(api_key = api_key)
        features <- manager$get_features(api_key = api_key)

        end_time <- Sys.time()
        duration <- round(as.numeric(end_time - start_time), 1)

        list(
          accession = accession,
          status = "success",
          duration_seconds = duration,
          metadata_rows = nrow(metadata),
          features_rows = ifelse(!is.null(features), nrow(features), 0),
          sequence_length = ifelse(
            !is.null(metadata$sequence) && nchar(metadata$sequence[1]) > 0,
            nchar(metadata$sequence[1]),
            0
          ),
          error = NULL
        )
      }, error = function(e) {
        list(
          accession = accession,
          status = "failed",
          duration_seconds = 0,
          metadata_rows = 0,
          features_rows = 0,
          sequence_length = 0,
          error = e$message
        )
      })

      # Display result
      if (result$status == "success") {
        message("    Success (", result$duration_seconds, "s)")
        message("     Metadata: ", result$metadata_rows, " rows")
        message("     Features: ", result$features_rows, " rows")
        message("     Sequence: ", result$sequence_length, " bp")
      } else {
        message("   Failed: ", result$error)
      }

      # Small delay between requests
      if (index < total) {
        Sys.sleep(1)
      }

      return(result)
    },

    process_all = function(parallel = FALSE, max_workers = 4, api_key = NULL) {
      total <- length(self$accessions)
      message("Processing ", total, " genomes...")

      results <- list()

      if (parallel && requireNamespace("parallel", quietly = TRUE)) {
        # Parallel processing
        cl <- parallel::makeCluster(min(max_workers, length(self$accessions)))
        results <- parallel::parLapply(cl, seq_along(self$accessions), function(i) {
          accession <- self$accessions[i]
          self$process_single(accession, i, total, api_key)
        })
        parallel::stopCluster(cl)
      } else {
        # Sequential processing
        for (i in seq_along(self$accessions)) {
          results[[i]] <- self$process_single(self$accessions[i], i, total, api_key)
        }
      }

      return(self$summarize_results(results))
    },

    summarize_results = function(results) {
      success_count <- sum(sapply(results, function(x) x$status == "success"))
      failed_count <- length(results) - success_count

      total_metadata <- sum(sapply(results, function(x) x$metadata_rows))
      total_features <- sum(sapply(results, function(x) x$features_rows))
      total_sequence <- sum(sapply(results, function(x) x$sequence_length))
      total_duration <- sum(sapply(results, function(x) x$duration_seconds))

      failed_accessions <- sapply(
        results[sapply(results, function(x) x$status == "failed")],
        function(x) x$accession
      )

      summary <- list(
        total_processed = length(results),
        success_count = success_count,
        failed_count = failed_count,
        total_metadata_rows = total_metadata,
        total_features_rows = total_features,
        total_sequence_bp = total_sequence,
        total_duration_seconds = total_duration,
        failed_accessions = failed_accessions,
        details = results
      )

      class(summary) <- "batch_summary"
      return(summary)
    }
  )

  processor$self <- processor
  class(processor) <- "batch_processor"

  return(processor)
}

#' @export
print.batch_summary <- function(x, ...) {
  cat("Batch Processing Summary\n")
  cat("=======================\n")
  cat("Total processed: ", x$total_processed, "\n")
  cat("Successful:      ", x$success_count, "\n")
  cat("Failed:          ", x$failed_count, "\n")
  cat("Total metadata:  ", x$total_metadata_rows, " rows\n")
  cat("Total features:  ", x$total_features_rows, " rows\n")
  cat("Total sequence:  ", x$total_sequence_bp, " bp\n")
  cat("Total duration:  ", x$total_duration_seconds, " seconds\n")

  if (length(x$failed_accessions) > 0) {
    cat("\nFailed accessions:\n")
    cat(paste(" ", head(x$failed_accessions, 5), collapse = "\n"))
    if (length(x$failed_accessions) > 5) {
      cat("\n  ... and", length(x$failed_accessions) - 5, "more")
    }
  }
}

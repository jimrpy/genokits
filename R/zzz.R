# Package initialization

.onLoad <- function(libname, pkgname) {
  # Set default options if not already set
  if (is.null(getOption("genokits.genome_db"))) {
    options(genokits.genome_db = NULL)
  }

  # Add package to the search path
  # This ensures our functions are available
  invisible()
}

.onAttach <- function(libname, pkgname) {
  # Welcome message with database info
  db_path <- get_genome_db()
  packageStartupMessage(
    "genokits loaded. Genome database: ", db_path, "\n",
    "Use set_genome_db() to change the database location."
  )
}

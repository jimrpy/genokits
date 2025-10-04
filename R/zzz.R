# Package startup messages
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Welcome to genokits ",
    utils::packageVersion("genokits"),
    "!\n",
    "Use set_ncbi_api_key() to set your NCBI API key for better rate limits."
  )
}

# Clean up on package unload
.onUnload <- function(libpath) {
  # Clean up any temporary files if needed
}

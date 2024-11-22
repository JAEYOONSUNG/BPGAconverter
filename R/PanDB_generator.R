#' PanDB Generator
#'
#' This function generates a global list (\code{PanDB}) containing processed data frames for each GenBank file
#' in a specified directory. Additionally, it creates a flattened data frame (\code{PanDB_df}) for further analysis.
#'
#' @param directory A character string specifying the path to the directory containing GenBank files.
#'        If \code{NULL}, the current working directory (\code{getwd()}) is used. If no GenBank files are found
#'        in the current working directory, the function searches for a subdirectory named "genbank".
#' @param save_output Logical, whether to save outputs for each processed file. Default is \code{FALSE}.
#'
#' @details
#' The function processes GenBank files in the specified directory by:
#' \itemize{
#'   \item Identifying files with extensions \code{.gbk}, \code{.gb}, or \code{.gbff}.
#'   \item Processing each file using the \code{Genbank_organizer} function.
#'   \item Storing the processed data frames in a global list (\code{PanDB}) with filenames as keys.
#'   \item Flattening the list into a single data frame (\code{PanDB_df}) for convenience.
#' }
#'
#' The function checks if the specified directory exists and contains valid GenBank files. If not, it stops with an error.
#' Errors during file processing are logged, and the function continues processing the remaining files.
#'
#' @return
#' The function assigns the following objects to the global environment:
#' \itemize{
#'   \item \code{PanDB}: A list of processed GenBank data frames, indexed by filename (excluding extensions).
#'   \item \code{PanDB_df}: A flattened data frame containing all processed data, with an additional column for filenames.
#' }
#'
#' @export
#'
#' @examples
#' # Example usage:
#' PanDB_generator(directory = "path/to/genbank/files")
#'
#' # Automatically use the current working directory:
#' PanDB_generator()
#'
#' @seealso \code{\link{Genbank_organizer}}
#'

PanDB_generator <- function(directory = NULL, save_output = FALSE) {
  # If directory is NULL, set to 'getwd()' or search for 'genbank' folder
  if (is.null(directory)) {
    possible_directory <- list.files(getwd(), pattern = "genbank", full.names = TRUE, ignore.case = TRUE)
    if (length(possible_directory) > 0 && dir.exists(possible_directory[1])) {
      directory <- possible_directory[1]
      message("Directory set to: ", directory)
    } else {
      directory <- getwd()
      message("Directory set to current working directory: ", directory)
    }
  }

  # Check if the directory exists
  if (!dir.exists(directory)) {
    stop("The specified directory does not exist.")
  }

  # List all GenBank files in the directory
  genbank_files <- list.files(directory, pattern = "\\.gbk$|\\.gb$|\\.gbff$", full.names = TRUE)
  if (length(genbank_files) == 0) {
    stop("No GenBank files found in the specified directory.")
  }

  # Initialize PanDB as an empty list in the global environment
  assign("PanDB", list(), envir = .GlobalEnv)

  # Placeholder for error tracking
  errors <- list()

  # Loop through each GenBank file
  for (file in genbank_files) {
    message("Processing file: ", file)
    tryCatch({
      # Process the file using Genbank_organizer
      Genbank_organizer(gb_dir = file, save_output = save_output)

      # Extract filename without extension
      file_name <- tools::file_path_sans_ext(basename(file))

      # Add genbank_table to PanDB in the global environment
      PanDB <- get("PanDB", envir = .GlobalEnv)
      PanDB[[file_name]] <- get("genbank_table", envir = .GlobalEnv)
      assign("PanDB", PanDB, envir = .GlobalEnv)

      # Remove genbank_table from environment
      rm(list = "genbank_table", envir = .GlobalEnv)
    }, error = function(e) {
      # If an error occurs, store the error message
      errors[[basename(file)]] <- conditionMessage(e)
      warning("Error processing file: ", file, "\n", conditionMessage(e))
    })
  }

  # Summary of results and errors
  successful_files <- names(PanDB)
  failed_files <- names(errors)

  message("\nProcessing complete.")
  message("Successfully processed files: ", length(successful_files))
  if (length(failed_files) > 0) {
    message("Files with errors: ", length(failed_files))
    for (file in failed_files) {
      message("- ", file, ": ", errors[[file]])
    }
  }
  # data.frame format (flatten)
  PanDB_df <- dplyr::bind_rows(PanDB, .id = "File")
  assign("PanDB_df", PanDB_df, envir = .GlobalEnv)
}

#' Combine GenBank Genome Information
#'
#' This function processes multiple GenBank files from a specified directory, extracts genome information
#' using the `gb_info` function, and compiles the results into a summary data frame. The processed data
#' are saved in the global environment as `Genome_summary_list` and `Genome_summary_df`.
#'
#' @param directory A character string specifying the path to the directory containing GenBank files.
#'        If \code{NULL}, the current working directory (\code{getwd()}) or a subdirectory named
#'        "genbank" will be used.
#' @param save_output Logical, whether to save outputs for each processed file. Default is \code{FALSE}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Searches the specified directory for GenBank files with extensions \code{.gbk}, \code{.gb}, or \code{.gbff}.
#'   \item Calls the \code{gb_info} function to process each GenBank file and extract genome information.
#'   \item Compiles the results into a list (\code{Genome_summary_list}) and a flattened data frame (\code{Genome_summary_df}).
#'   \item Saves both the list and data frame in the global environment.
#' }
#'
#' If the specified directory does not contain GenBank files, or if the directory is invalid, the function will stop with an error.
#'
#' @return A data frame \code{Genome_summary_df} containing combined genome information from all processed GenBank files.
#'
#' @export
#'
#' @examples
#' # Example usage
#' combine_gb_info(directory = "path/to/genbank/files")
#'
#' # Automatically detect the 'genbank' folder or use the current working directory
#' combine_gb_info()
#'
#' @seealso \code{\link{gb_info}}, \code{\link{Genbank_organizer_batch}}
#'

Genbank_info_combiner <- function(directory = NULL, save_output = FALSE) {
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

  Genome_summary_list <- list() # Initialize Genome_summary_list as an empty list
  errors <- list() # Placeholder for error tracking

  # Loop through each GenBank file
  for (file in genbank_files) {
    message("Processing file: ", file)
    tryCatch({
      gb_info(gb_dir = file) # Call gb_info with the file path
      file_name <- tools::file_path_sans_ext(basename(file)) # Extract filename without extension
      Genome_summary <- get("Genome_summary", envir = .GlobalEnv) # Retrieve Genome_summary from the global environment
      Genome_summary_list[[file_name]] <- Genome_summary # Append the result for the current file
      rm(list = "Genome_summary", envir = .GlobalEnv) # Remove the single-file Genome_summary from the environment
    }, error = function(e) {
      errors[[basename(file)]] <- conditionMessage(e) # If an error occurs, store the error message
      warning("Error processing file: ", file, "\n", conditionMessage(e))
    })
  }

  # Summary of results and errors
  successful_files <- names(Genome_summary_list)
  failed_files <- names(errors)

  # Flatten the list into a data.frame
  Genome_summary_df <- dplyr::bind_rows(Genome_summary_list)

  message("\nProcessing complete.")
  message("Successfully processed files: ", length(successful_files))
  if (length(failed_files) > 0) {
    message("Files with errors: ", length(failed_files))
    for (file in failed_files) {
      message("- ", file, ": ", errors[[file]])
    }
  }

  # Assign Genome_summary_list to the global environment
  assign("Genome_summary_list", Genome_summary_list, envir = .GlobalEnv)
  assign("Genome_summary_df", Genome_summary_df, envir = .GlobalEnv)

  return(Genome_summary_df) # Return the data frame for immediate use
}

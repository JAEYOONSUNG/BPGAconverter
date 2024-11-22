#' Extract Genome Information from GenBank Files
#'
#' This function processes GenBank files to extract genome metadata and annotation information,
#' generating a summarized genome table. The function supports `.gb`, `.gbk`, and `.gbff` file formats.
#'
#' @param gb_dir A character string specifying the directory containing GenBank files, or the path to a single GenBank file. Defaults to the current working directory if `NULL`.
#' @return A data frame (`Genome_summary`) containing the extracted genome information. The data frame is also assigned to the global environment.
#' @export
#'
#' @examples
#' # Extract genome information from GenBank files in the current directory
#' gb_info()
#'
#' # Extract genome information from a specific directory
#' gb_info(gb_dir = "path/to/genbank/files")
#'
#' # Extract genome information from a single GenBank file
#' gb_info(gb_dir = "path/to/single/file.gb")

gb_info <- function(gb_dir = NULL) {
  # Load necessary libraries
  library(dplyr)
  library(readr)
  library(tidyr)
  library(openxlsx)
  library(plyr)

  # Initialize empty data frame for storing genome table
  Genome_summary <- data.frame()

  # Read the specified GenBank file
  if (is.null(gb_dir)) {
    gb_dir <- getwd()
    # Get the list of GenBank files in the directory
    gb_files <- list.files(
      path = gb_dir,
      pattern = "\\.gbk$|\\.gb$|\\.gbff$",
      full.names = TRUE
    )
  } else {
    # If gb_dir is provided, use it as the file directly
    gb_files <- gb_dir
  }
  # Exclude temporary files starting with '~$'
  gb_files <- gb_files[!grepl("^~\\$", gb_files)]

  # Check if InterPro search files exist
  if (length(gb_files) == 0) {
    stop("No genbank files found.")
  }

  #origin <- readr::read_fwf(gb_files, col_positions = fwf_empty(gb_files), show_col_types = FALSE)
  # Step 1: Attempt to read the GenBank file
  origin <- tryCatch({
    readr::read_fwf(
      gb_files,
      col_positions = fwf_empty(gb_files),
      show_col_types = FALSE
    )
  }, error = function(e) {
    stop("Error reading GenBank file: ", gb_files, "\n", e)
  })

  # Step 2: Check if the data was split into at least two columns
  if (ncol(origin) < 2) {
    warning("File ", gb_files, " has fewer than 2 columns. Attempting to split manually.")

    # Step 3: Manually split using fixed positions
    origin <- readr::read_fwf(
      gb_files,
      col_positions = fwf_positions(
        start = c(1, 13),  # Adjust based on expected GenBank structure
        end = c(12, NA),
        col_names = c("X1", "X2")
      ),
      show_col_types = FALSE
    )
  }

  # Extract portion before "FEATURES"
  origin <- origin %>%
    dplyr::slice(., 1:as.numeric(grep("FEATURES", origin$X1, ignore.case = FALSE)[1]) - 1) %>%
    sapply(function(x) stringr::str_squish(x)) %>%
    as.data.frame()

  # Process genome info table
  Genome_info_table <- origin %>%
    dplyr::slice(., 1:as.numeric(grep("##", origin$X2, ignore.case = FALSE)[1]) - 1) %>%
    dplyr::mutate_all(~replace(., is.na(.), NA)) %>%
    tidyr::fill(X1) %>%
    dplyr::group_by(X1) %>%
    dplyr::summarise(X2 = paste(X2, collapse = " ")) %>%
    dplyr::mutate(X2 = ifelse(X2 == "NA", "", X2))

  # Process genome annotation table
  Genome_annotation_table <- origin %>%
    dplyr::slice(as.numeric(grep("^##Genome-Annotation-Data-START##", .$X2, ignore.case = FALSE)[1] + 1):
            as.numeric(grep("^##Genome-Annotation-Data-END##", .$X2, ignore.case = FALSE)[1] - 1)) %>%
    tidyr::separate(col = "X2", into = c("X1", "X2"), sep = " :: ", fill = "left", extra = "merge") %>%
    dplyr::mutate_all(~replace(., is.na(.), NA)) %>%
    tidyr::fill(X1)   %>%
    dplyr::group_by(X1) %>%
    dplyr::summarise(X2 = paste(X2, collapse = " "))

  # Combine info and annotation tables
  Genome_summary <- rbind(Genome_info_table, Genome_annotation_table) %>% as.data.frame()

  # Remove empty first row if necessary
  Genome_summary <- Genome_summary %>%
    dplyr::slice(if (is.na(Genome_summary[1, 1]) | Genome_summary[1, 1] == "") 2:n() else 1:n())

  # Final processing on the genome table
  rownames(Genome_summary) <- Genome_summary$X1
  Genome_summary <- Genome_summary %>% dplyr::select(-c(X1))
  Genome_summary <- rbind("File" = gsub("\\.gb(k|ff)?$", "", basename(gb_files)), Genome_summary)
  Genome_summary <- t(Genome_summary) %>% as.data.frame()
  rownames(Genome_summary) <- NULL

  # Assign the table to the global environment (optional)
  assign("Genome_summary", Genome_summary, envir = .GlobalEnv)
}

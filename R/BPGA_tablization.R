#' BPGA tablization
#'
#' @param supporting_files_dir Path to the directory containing Supporting_files. Defaults to the current working directory.
#' @return A reconstructed data frame.
#' @export
#'

BPGA_tablization <- function(supporting_files_dir = NULL) {
  # Set the directory for Supporting_files
  if (is.null(supporting_files_dir)) {
    supporting_files_dir <- file.path(getwd(), "Supporting_files")
  }

  # Check if directory exists
  if (!dir.exists(supporting_files_dir)) {
    stop("Supporting_files directory not found at: ", supporting_files_dir)
  }

  # Load strain list
  list_file <- file.path(supporting_files_dir, "list")
  if (!file.exists(list_file)) {
    stop("list file not found in Supporting_files directory.")
  }

  list <- read.delim(list_file, sep = "\t", header = FALSE, col.names = c("Strain_number", "Species"))

  # Load u_clusters data
  u_clusters_file <- file.path(supporting_files_dir, "u_clusters.txt")
  if (!file.exists(u_clusters_file)) {
    stop("u_clusters.txt file not found in Supporting_files directory.")
  }

  raw <- read.csv(u_clusters_file, sep = "\t", header = FALSE)

  # Check uclust code summary
  uclust_summary <- raw %>%
    dplyr::group_by(uclust_code = V1) %>%
    dplyr::summarise(count = n())
  print(uclust_summary)

  # Extract and process relevant columns
  trim <- raw[, c(1, 2, 3, 4, 7, 8, 9, 10)]
  colnames(trim) <- c("Hit", "Gene_number", "AA_length", "Sequence_id", "hitC", "Variable_region", "Multi", "Ref_seq")

  # Split the "Multi" column into separate components
  trim <- tidyr::separate(trim, col = "Multi", into = c("Strain_number", "Accession_number", "Species"), sep = "\\|")

  # Split the "Species" column to extract GC content
  trim <- tidyr::separate(trim, col = "Species", into = c("Species", "GC"), sep = "_GC=")

  # Clean up GC column
  trim$GC <- sapply(trim$GC, function(x) as.integer(gsub("[#*]", "", x)))

  # Filter rows where Hit contains "S" or "H"
  trim <- trim[grep("S|H", trim$Hit, ignore.case = TRUE), ]

  # Sort and remove duplicates
  trim <- trim %>%
    dplyr::arrange(dplyr::desc(Sequence_id)) %>%
    tidyr::unite(hitC, Gene_number, Species, sep = "_", remove = FALSE) %>%
    dplyr::distinct(hitC, .keep_all = TRUE)

  # Sort the final table
  trim <- dplyr::arrange(trim, trim$Gene_number)

  # Convert to data frame
  trim <- as.data.frame(trim)

  # Assign contig_length to the global environment
  assign("u_cluster", trim, envir = .GlobalEnv)
  message("'u_cluster' has been generated and assigned to the global environment.")

  # Return the processed table
  return(dplyr::glimpse(trim))
}

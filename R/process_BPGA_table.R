#' Process Table
#'
#' This function processes and filters a final data table by merging input data frames,
#' performing filtering and validation steps, and saving the result in the global environment as \code{PanDB_df}.
#'
#' @param u_cluster A data frame containing cluster information. If \code{NULL}, the function will try
#'        to retrieve \code{u_cluster} from the global environment.
#' @param PanDB_df A data frame containing processed PanDB information. If \code{NULL}, the function will
#'        attempt to retrieve \code{PanDB_df} from the global environment.
#' @param Genome_summary_df A data frame containing genome summary information. If \code{NULL}, the function
#'        will attempt to retrieve \code{Genome_summary_df} from the global environment.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Merges \code{u_cluster} with \code{PanDB_df}.
#'   \item Appends genome summary data from \code{Genome_summary_df}.
#'   \item Filters and removes duplicated or unmatched sequences based on cluster and genome data.
#'   \item Performs validation to ensure all processed files are included.
#'   \item Saves the final processed table as \code{PanDB_df} in the global environment.
#' }
#'
#' @return The function does not return a value directly but assigns the processed \code{PanDB_df}
#'         to the global environment.
#'
#' @export
#'
#' @examples
#' # Example usage
#' process_final_table(u_cluster = my_cluster, PanDB_df = my_PanDB_df, Genome_summary_df = my_genome_summary_df)
#'
#' # Automatically fetch objects from the global environment
#' process_final_table()
#'
#' @seealso \code{\link{Genbank_organizer_batch}}, \code{\link{combine_gb_info}}
#'

process_BPGA_table <- function(u_cluster = NULL, PanDB_df = NULL, Genome_summary_df = NULL) {
  # Load objects from the global environment if not provided
  if (is.null(u_cluster)) {
    if (exists("u_cluster", envir = .GlobalEnv)) {
      u_cluster <- get("u_cluster", envir = .GlobalEnv)
    } else {
      stop("u_cluster is not provided and does not exist in the global environment.")
    }
  }

  if (is.null(PanDB_df)) {
    if (exists("PanDB_df", envir = .GlobalEnv)) {
      PanDB_df <- get("PanDB_df", envir = .GlobalEnv)
    } else {
      stop("PanDB_df is not provided and does not exist in the global environment.")
    }
  }

  if (is.null(Genome_summary_df)) {
    if (exists("Genome_summary_df", envir = .GlobalEnv)) {
      Genome_summary_df <- get("Genome_summary_df", envir = .GlobalEnv)
    } else {
      stop("Genome_summary_df is not provided and does not exist in the global environment.")
    }
  }

  # Step 1: Merge clusters with PanDB_df
  Final <- merge(u_cluster, PanDB_df, by.x = "Accession_number", by.y = "protein_id", all.x = TRUE) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate("Species" = gsub("\\_", " ", .$Species))
  Final <- merge(Final, Genome_summary_df %>% dplyr::select(File, SOURCE), by = "File", all.x = TRUE)

  # Step 2: Create and merge locus_prefix
  locus_prefix <- PanDB_df %>%
    dplyr::select(locus_tag, File) %>%
    dplyr::mutate(locus_prefix = qdap::beg2char(locus_tag, "_")) %>%
    dplyr::distinct(locus_prefix, .keep_all = TRUE) %>%
    dplyr::select(-locus_tag)

  Final <- merge(Final, locus_prefix, by = "File", all.x = TRUE)

  # Step 3: Filter duplicated matching sequences
  Final <- Final %>%
    dplyr::mutate(SOURCE_temp = stringr::str_replace_all(SOURCE, "[^a-zA-Z0-9\\.]", " ")) %>%
    dplyr::filter(SOURCE_temp == Species) %>%
    dplyr::select(-SOURCE_temp)

  # Step 4: Check for matches
  identifier <- Final %>%
    dplyr::select(File, Species, SOURCE, contig) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::mutate(isMatch = sapply(
      strsplit(File, " "),
      function(words) all(grepl(paste0("\\b", words, "\\b", collapse = "|"), contig))
    ))

  # Step 5: Validate processed files
  unmatched_files <- setdiff(
    Final %>% dplyr::distinct(File) %>% dplyr::pull(),
    PanDB_df %>% dplyr::distinct(File) %>% dplyr::pull()
  )

  if (length(unmatched_files) == 0) {
    message("Go to next step")
    PanDB_df <- Final
  } else {
    message("Check the group_by factors and follow strains")
    print(list[unmatched_files, ])
  }

  # Step 6: Display results
  message("No. of processed strains: ", PanDB_df %>% dplyr::distinct(File) %>% dplyr::n_distinct())
  message("Processed strains:\n", paste(PanDB_df %>% distinct(SOURCE) %>% pull(), collapse = "\n"))

  # Save the final PanDB_df to the global environment
  assign("PanDB_df", PanDB_df, envir = .GlobalEnv)
}

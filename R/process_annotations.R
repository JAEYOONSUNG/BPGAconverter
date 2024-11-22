#' Process Gene and COG Annotations
#'
#' This function processes `.annotations.` files found in a specified directory, extracts gene and COG category annotations,
#' and concatenates the results into structured outputs. The outputs include gene and COG annotations grouped by clusters.
#'
#' @param directory A character string specifying the base directory for processing. If \code{NULL}, the current working directory 
#'        (\code{getwd()}) will be used, and the function will validate the presence of required subdirectories 
#'        (\code{"Results"} and \code{"Supporting_files"}).
#' @param PanDB A global environment list containing the processed GenBank files. If \code{NULL}, the function assumes the object 
#'        exists in the global environment.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Checks for the presence of `.annotations.` files within the specified directory or its subdirectories.
#'   \item Extracts gene and COG annotations by concatenating relevant information from the provided `PanDB`.
#'   \item Outputs cleaned and concatenated gene (`Gene_conc.`) and COG category (`COG_conc.`) annotations.
#'   \item Skips processing if no `.annotations.` files are detected.
#' }
#'
#' The function validates the directory structure to ensure the required subdirectories are present. If they are missing, the 
#' function stops with an error.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Gene_conc}: A data frame of concatenated gene annotations.
#'   \item \code{COG_conc}: A data frame of concatenated COG category annotations.
#' }
#' If no `.annotations.` files are found, the function returns \code{NULL}.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' process_annotations(directory = "path/to/directory", PanDB = PanDB)
#'
#' # Automatically detect and process using the current working directory:
#' process_annotations()
#'
#' @seealso \code{\link{Genbank_organizer_batch}}, \code{\link{combine_gb_info}}
#'

# add gene and COG annotations per cluster
process_annotations <- function(directory = NULL, PanDB = NULL) {
  # If directory is NULL, set to 'getwd()' and validate subdirectories
  if (is.null(directory)) {
    base_dir <- getwd()
    message("Directory set to current working directory: ", base_dir)
    
    # Check if required subdirectories exist
    required_dirs <- c("Results", "Supporting_files")
    missing_dirs <- required_dirs[!dir.exists(file.path(base_dir, required_dirs))]
    
    if (length(missing_dirs) > 0) {
      stop("The following required subdirectories are missing: ", paste(missing_dirs, collapse = ", "))
    }
    
    directory <- base_dir
  }
  
  # Check for .annotations. files in the directory
  annotation_files <- list.files(directory, pattern = "\\.annotations\\.", full.names = TRUE, recursive = TRUE)
  
  if (length(annotation_files) > 0) {
    message("Detected .annotations. files. Proceeding with processing.")
    
    # Concatenate detected gene names in set
    Gene_conc. <- tryCatch({
      purrr::map2_dfc(
        casting_by_gene[, 2:(length(PanDB) + 1)],
        casting_by_Preferred_name[, 2:(length(PanDB) + 1)],
        ~ str_c(.x, .y, sep = ",")
      ) %>%
        dplyr::mutate(
          "Gene_number" = casting$Gene_number,
          tidyr::unite(., "Gene_conc.", 1:length(PanDB), sep = ",")
        ) %>%
        dplyr::mutate("Gene_conc." = gsub("Not assigned|\\-", "", Gene_conc.)) %>%
        dplyr::mutate("Gene_conc." = gsub("^,{1,}|,{1,}$", "", Gene_conc.)) %>%
        dplyr::mutate("Gene_conc." = gsub("(?<=[\\,])\\,*|^\\,+|\\,+$", "", Gene_conc., perl = TRUE)) %>%
        dplyr::select(Gene_number, Gene_conc.) %>%
        splitstackshape::cSplit(., 'Gene_conc.', ',') %>%
        apply(., 1, unique) %>%
        data.frame(t(sapply(., `length<-`, max(lengths(.))))) %>%
        unite(., "Gene_conc.", 2:length(.), sep = ",", na.rm = TRUE) %>%
        setNames(c("Gene_number", "Gene_conc."))
    }, error = function(e) {
      warning("Error in processing Gene_conc.: ", conditionMessage(e))
      NULL
    })
    
    # Concatenate detected COG_category in set
    COG_conc. <- tryCatch({
      as.data.frame(casting_by_COG_category[, 2:(length(PanDB) + 1)]) %>%
        dplyr::mutate(
          "Gene_number" = casting$Gene_number,
          tidyr::unite(., "COG_conc.", 1:length(PanDB), sep = ",")
        ) %>%
        dplyr::mutate("COG_conc." = gsub("Not assigned|\\-", "", COG_conc.)) %>%
        dplyr::mutate("COG_conc." = gsub("^,{1,}|,{1,}$", "", COG_conc.)) %>%
        dplyr::mutate("COG_conc." = gsub("(?<=[\\,])\\,*|^\\,+|\\,+$", "", COG_conc., perl = TRUE)) %>%
        dplyr::select(Gene_number, COG_conc.) %>%
        splitstackshape::cSplit(., 'COG_conc.', ',') %>%
        apply(., 1, unique) %>%
        data.frame(t(sapply(., `length<-`, max(lengths(.))))) %>%
        tidyr::unite(., "COG_conc.", 2:length(.), sep = ",", na.rm = TRUE) %>%
        setNames(c("Gene_number", "COG_conc."))
    }, error = function(e) {
      warning("Error in processing COG_conc.: ", conditionMessage(e))
      NULL
    })
    
    # Message to confirm successful processing
    if (!is.null(Gene_conc.) && !is.null(COG_conc.)) {
      message("Successfully processed Gene_conc. and COG_conc.")
    }
    
    return(list(Gene_conc = Gene_conc., COG_conc = COG_conc.))
  } else {
    message("No .annotations. files detected. Skipping processing.")
    return(NULL)
  }
}

# To be updated
# InterPro <- reshape2::melt(InterPro_search, id.vars=c("Protein accession","Analysis"),measure.vars=c("Signature accession","Signature description","Score"))
# InterPro <- reshape2::dcast(InterPro, formula = `Protein accession`  ~ variable + Analysis, fun.aggregate = toString)

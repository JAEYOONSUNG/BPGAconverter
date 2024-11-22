#' BPGA Conversion and Table Generation
#'
#' This function orchestrates the generation of comparative genomics tables using BPGA pipeline outputs.
#' It processes PanDB data, combines genome summary information, and generates factor-specific tables.
#'
#' @param directory A character string specifying the working directory where the required files and outputs are located. Defaults to the current working directory.
#' @param save_output Logical. If `TRUE` (default is `FALSE`), intermediate and final outputs will be saved to files.
#' @param Table_factor A character vector of column names from `PanDB_df` used for generating factor-specific tables. Defaults to `c("locus_tag", "product", "Sequence_id")`.
#' @return A combined table (`casting_merge`) of comparative genomics data, saved as an Excel file if `save_output` is `TRUE`.
#' @export
#'
#' @examples
#' # Run the function with default settings
#' run_BPGAconvert()
#'
#' # Run with a specified directory and save the outputs
#' run_BPGAconvert(directory = "path/to/working/directory", save_output = TRUE)
#'
#' # Generate tables for specific factors
#' run_BPGAconvert(Table_factor = c("locus_tag", "Preferred_name"))

run_BPGAconvert <- function(
    directory = NULL,
    save_output = FALSE,
    Table_factor = c("locus_tag", "product", "Sequence_id")
) {
  # Step 1: Generate PanDB
  PanDB_generator(directory = directory, save_output = save_output)

  # Step 2: Combine Genome Summary Info
  Genbank_info_combiner(directory = directory, save_output = save_output)

  # Step 3: Uclust output file
  BPGA_tablization()

  # Step 4: Process Final Table
  process_BPGA_table()

  # Step 5: Table generation by factors
  factor <- c("locus_tag", "Accession_number", "product", "Sequence_id", "gene", "Preferred_name", "COG_category")
  for (i in 1:length(factor)) {
    if (factor[i] %in% colnames(PanDB_df)) {
      tryCatch({
        casting <- reshape2::dcast(PanDB_df, formula = Gene_number ~ Strain_number,
                                   value.var = factor[i],
                                   fun.aggregate = function(x) paste(x, collapse = ", "),
                                   drop = TRUE)
        casting <- casting %>% dplyr::mutate_at(vars(Gene_number), as.integer)
        casting <- dplyr::arrange(casting, casting$Gene_number)
        casting <- casting[ ,c("Gene_number", c(1:length(PanDB)))]
        names(casting) <- Genome_summary_df$File[match(names(casting), rownames(Genome_summary_df))] # name matching
        names(casting)[1] <- "Gene_number"
        casting <- casting %>% dplyr::mutate(across(everything(), ~ stringr::str_replace(., "^NA$", "Not assigned")))

        assign(paste("casting", "by", factor[i], sep = "_"), casting)
        message("Successfully created casting for: ", factor[i])
      }, error = function(e) {
        warning("Error processing casting for: ", factor[i], "\n", conditionMessage(e))
      })
    } else {
      message("Skipping ", factor[i], ": Column not found in PanDB_df.")
    }
  }

  # Step 6: add on EggNOG annotations [Optional]
  process_annotations()


  # Step 7: Factor-specific casting
  factor <- c("locus_tag", "Accession_number", "product", "Sequence_id", "gene", "Preferred_name", "COG_category")
  for (f in factor) {
    # Create the object name dynamically
    casting_name <- paste("casting_by", f, sep = "_")

    # Check if the object exists in the environment
    if (exists(casting_name)) {
      tryCatch({
        # Get the object and modify its column names
        casting_object <- get(casting_name)
        names(casting_object)[2:ncol(casting_object)] <- paste(names(PanDB), f, sep = "_")

        # Reassign the modified object back to its original name
        assign(casting_name, casting_object, envir = .GlobalEnv)
        message("Successfully renamed columns for: ", casting_name)
      }, error = function(e) {
        warning("Error renaming columns for: ", casting_name, "\n", conditionMessage(e))
      })
    } else {
      message("Skipping ", casting_name, ": Object not found.")
    }
  }

  casting_merge <- u_cluster %>% dplyr::select("Gene_number") %>% dplyr::distinct()

  # Initialize casting_merge with the first table
  casting_merge <- get(paste("casting", "by", Table_factor[1], sep = "_"))

  # Sequentially merge tables
  for (i in 2:length(Table_factor)) {
    current_table <- get(paste("casting", "by", Table_factor[i], sep = "_"))

    casting_merge <- merge(
      x = casting_merge,
      y = current_table,
      by = "Gene_number",
      suffixes = c("", paste0("_", Table_factor[i])) # Explicitly add suffixes
    )
  }

  casting_merge <- casting_merge[, c("Gene_number", sort(setdiff(names(casting_merge), "Gene_number")))]
  casting_merge <- casting_merge %>% dplyr::mutate_at(vars(Gene_number), as.integer)
  casting_merge <- dplyr::arrange(casting_merge, casting_merge$Gene_number)
  casting_merge <- data.table::setcolorder(casting_merge, order(sub('_.*', '', names(casting_merge))))
  assign("casting_merge", casting_merge, envir = .GlobalEnv)

  # Step 8: Conditional merge with Gene_conc. and COG_conc.
  if (exists("Gene_conc.", envir = .GlobalEnv)) {
    casting_merge <- merge(x = casting_merge, y = get("Gene_conc.", envir = .GlobalEnv), by = "Gene_number")
    message("Gene_conc. merged successfully.")
  } else {
    message("Gene_conc. not found in the environment. Skipping merge.")
  }

  if (exists("COG_conc.", envir = .GlobalEnv)) {
    casting_merge <- merge(x = casting_merge, y = get("COG_conc.", envir = .GlobalEnv), by = "Gene_number")
    message("COG_conc. merged successfully.")
  } else {
    message("COG_conc. not found in the environment. Skipping merge.")
  }

  # Step 9: Save the final table (optional)
  output_filename <- paste0("Table. comparative genomics ", length(PanDB), ".xlsx")
  xlsx::write.xlsx(casting_merge, output_filename, row.names = FALSE, col.names = TRUE, sheetName = "BPGA")
  message("Final table saved as: ", output_filename)
}

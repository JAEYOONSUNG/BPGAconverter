#' Generate Combined Plot for Flower, Pan-Core, and Phylogeny
#'
#' This function generates a combined plot of the flower plot, pan-core plot, and phylogeny plot.
#' The plots are labeled "a", "b", and "c" respectively, and are arranged in a grid format.
#'
#' @param directory A character string specifying the working directory where the required files and outputs are located. Defaults to the current working directory.
#' @param save_combined Logical. If `TRUE` (default), the combined plot will be saved as a PDF.
#' @param combined_filename A character string specifying the filename for the combined plot PDF. Defaults to `"combined_plot.pdf"`.
#' @return A combined plot object.
#' @export
#'
#' @examples
#' # Generate and save the combined plot in the current directory
#' combined_plot()
#'
#' # Generate the combined plot without saving
#' combined_plot(save_combined = FALSE)
#'
#' # Save the combined plot to a specific directory
#' combined_plot(directory = "path/to/output", combined_filename = "custom_combined_plot.pdf")

run_BPGAplot <- function(directory = NULL, save_combined = TRUE, combined_filename = "combined_plot.pdf") {
  # Step 1: Set default directory
  if (is.null(directory)) {
    directory <- getwd()
    message("Directory set to current working directory: ", directory)
  }

  # Step 2: Generate individual plots
  message("Generating flower plot...")
  flower <- flower_plot(save_plot = FALSE)  # Ensure flower_plot() returns a ggplot object

  message("Generating pan-core plot...")
  pan_core <- pan_core_plot(save_plot = FALSE)  # Do not save individual plot

  message("Generating phylogeny plot...")
  phylogeny <- phylogeny_plot(save_plot = FALSE, save_palette = FALSE)  # Do not save individual plot

  # Step 3: Resize flower_plot and arrange plots
  # Add labels a, b, c
  labeled_flower <- gridExtra::arrangeGrob(flower,
                                           top = grid::textGrob("a",
                                                                x = unit(0, "npc"),
                                                                just = "left",
                                                                gp = grid::gpar(fontsize = 16, fontface = "bold")))
  labeled_pan_core <- gridExtra::arrangeGrob(pan_core,
                                             top = grid::textGrob("b",
                                                                  x = unit(0, "npc"),
                                                                  just = "left",
                                                                  gp = grid::gpar(fontsize = 16, fontface = "bold")))
  labeled_phylogeny <- gridExtra::arrangeGrob(phylogeny,
                                              top = grid::textGrob("c",
                                                                   x = unit(0, "npc"),
                                                                   just = "left",
                                                                   gp = grid::gpar(fontsize = 16, fontface = "bold")))

  # Combine plots
  combined <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(
      labeled_flower,
      labeled_pan_core,
      ncol = 2,  # Arrange flower and pan-core plots in one row
      widths = c(0.75, 1)  # Adjust relative widths of the plots
    ),
    labeled_phylogeny,
    ncol = 1,  # Phylogeny plot occupies the next row
    heights = c(3, 2)  # Adjust the relative height of rows
  )

  # Step 4: Save combined plot if required
  if (save_combined) {
    combined_filepath <- file.path(directory, combined_filename)
    grDevices::pdf(combined_filepath, width = 16, height = 10)  # Adjust dimensions as needed
    grid::grid.draw(combined)
    grDevices::dev.off()
    message("Combined plot saved as: ", combined_filepath)
  }

  # Step 5: Return the combined plot
  return(combined)
}

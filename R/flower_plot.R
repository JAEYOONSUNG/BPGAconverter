#' Generate Flower Plot for Genome Analysis
#'
#' This function generates a flower plot to visualize core, unique, accessory, and absent genes for various organisms.
#'
#' @param file_path The file path to the input data (e.g., a `.xls` or `.csv` file). If `NULL`, the function will look for a default file in the current working directory under `"Results/stats.xls"`.
#' @param data An optional data frame containing the required columns: `'Organism name'`, `'No. of core genes'`, `'No. of unique genes'`, `'No. of accessory genes'`, and `'No. of exclusively absent genes'`. If provided, this will override the file input.
#' @param output_dir The directory where the resulting flower plot PDF will be saved. If `NULL`, the plot will be saved in the current working directory.
#' @param label_outer Logical. If `TRUE`, aligns labels for organisms on the outer side of the petals. Defaults to `TRUE`.
#' @param save_plot Logical. If `TRUE` (default), saves the plot as a PDF. If `FALSE`, the plot will not be saved.
#' @return A ggplot object representing the flower plot.
#' @export
#'
#' @examples
#' # Generate and save the flower plot in the current directory
#' flower_plot()
#'
#' # Use custom input data and save to a specific directory
#' flower_plot(data = custom_data, output_dir = "path/to/output")
#'
#' # Generate plot with inner-aligned labels without saving
#' flower_plot(file_path = "path/to/stats.xls", label_outer = FALSE, save_plot = FALSE)

flower_plot <- function(file_path = NULL, data = NULL, output_dir = NULL, label_outer = TRUE, save_plot = TRUE){

  # Check if file_path is NULL and read from default stats.xlsx
  if (is.null(file_path)) {
    file_path <- paste(getwd(),"Results/stats.xls", sep="/")  # Default file location
    file_path <- normalizePath(file_path)
  }

  # Read and preprocess data
  if (is.null(data)) {
    if (!file.exists(file_path)) {
      stop("File not found: Please provide a valid file path or data frame.")
    }
    data <- readr::read_delim(file_path) %>%
      mutate(`Organism name` = gsub("_", " ", `Organism name`))  # Replace underscores with spaces
  }

  # Check for required columns
  if (!all(c("Organism name", "No. of core genes", "No. of unique genes") %in% colnames(data))) {
    stop("Data must contain 'Organism name', 'No. of core genes', and 'No. of unique genes' columns.")
  }

  # Extract values for plotting
  categories <- data$`Organism name`
  core_genes <- data$`No. of core genes`  # Assume the first row contains the core gene count
  unique_genes <- data$`No. of unique genes`
  accessory_genes <- data$`No. of accessory genes`
  absent_genes <- data$`No. of exclusively absent genes`

  # Number of petals
  num_petals <- length(categories)

  # Fixed petal size
  petal_width <- num_petals / 4   # Width of the petal
  petal_height <- petal_width * 5  # Height of the petal

  # Create a data frame for angles and categories
  flower_data <- data.frame(
    category = categories,
    unique = unique_genes,
    accessory = accessory_genes,
    absent = absent_genes,
    angle = seq(0, 360, length.out = num_petals + 1)[-1]  # Angles for petals
  )

  # Function to generate petal coordinates (half ellipse)
  create_petal <- function(center_x, center_y, width, height, angle, label) {
    t <- seq(0, pi, length.out = 100)  # Half ellipse (0 to π)
    x <- width * cos(t)  # Ellipse x-coordinates
    y <- height * sin(t) # Ellipse y-coordinates

    # Rotate the petal to the desired angle
    rotated_x <- x * cos(angle) - y * sin(angle) + center_x
    rotated_y <- x * sin(angle) + y * cos(angle) + center_y

    data.frame(x = rotated_x, y = rotated_y, label = label)
  }

  # Generate petals at adjusted positions
  petals <- do.call(rbind, lapply(1:nrow(flower_data), function(i) {
    angle_rad <- flower_data$angle[i] * pi / 180  # Convert angle to radians

    # Shift the petal so its base touches the center
    shift_x <- -petal_height / 8 * sin(angle_rad)  # Adjust shift_x
    shift_y <- petal_height / 8 * cos(angle_rad)   # Adjust shift_y

    # Create the petal at the rotated position
    create_petal(
      center_x = shift_x,
      center_y = shift_y,
      width = petal_width,
      height = petal_height,
      angle = angle_rad,
      label = flower_data$category[i]
    )
  }))

  # Set radii for the circles
  core_circle_radius <- petal_width * 1.2  # Radius for the core circle
  accessory_circle_radius <- core_circle_radius * 2.0  # Radius for the accessory circle

  # Create data for the center circles
  core_circle_data <- data.frame(
    x = core_circle_radius * cos(seq(0, 2 * pi, length.out = 100)),
    y = core_circle_radius * sin(seq(0, 2 * pi, length.out = 100))
  )

  accessory_circle_data <- data.frame(
    x = accessory_circle_radius * cos(seq(0, 2 * pi, length.out = 100)),
    y = accessory_circle_radius * sin(seq(0, 2 * pi, length.out = 100))
  )

  flower_organism <- flower_data %>%
    dplyr::mutate(
      text_x = (petal_height * 1.15) * -sin(angle * pi / 180),  # Adjust text_x to align with petal base
      text_y = (petal_height * 1.15) * cos(angle * pi / 180),   # Adjust text_y similarly
      text_angle = ifelse(angle <= 180, angle - 90, angle + 90),  # Adjust text rotation to match petals
      hjust = if (label_outer) {  # Dynamically set hjust based on label_outer
        ifelse(angle <= 180, 1, 0)  # Outer-petal align
      } else {
        ifelse(angle <= 180, 0, 1)  # Inner-petal align
      }
    )


  flower_unique <- flower_data %>%
    dplyr::mutate(
      text_x = (petal_height * 0.9) * -sin(angle * pi / 180),  # Adjust text_x to align with petal base
      text_y = (petal_height * 0.9) * cos(angle * pi / 180),   # Adjust text_y similarly
      text_angle = ifelse(angle <= 180, angle - 90, angle + 90)  # Adjust text rotation to match petals
    )

  flower_accessory <- flower_data %>%
    dplyr::mutate(
      text_x = (petal_height * 0.36) * -sin(angle * pi / 180),  # Adjust text_x to align with petal base
      text_y = (petal_height * 0.36) * cos(angle * pi / 180),   # Adjust text_y similarly
      text_angle = ifelse(angle <= 180, angle - 90, angle + 90)  # Adjust text rotation to match petals
    )

  # Create flower plot
  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = petals, aes(x = x, y = y, fill = label, group = label), color = "white", alpha = 0.25) +
    ggplot2::geom_polygon(data = accessory_circle_data, aes(x = x, y = y), fill = "yellow", color = "white", alpha = 0.5) +  # Add accessory circle
    ggplot2::geom_polygon(data = core_circle_data, aes(x = x, y = y), fill = "white", color = "black") +  # Add core circle
    ggplot2::geom_text(data = flower_organism, aes(x = text_x, y = text_y, label = category, angle = text_angle, hjust = hjust), size = 5) +  # Corrected text rotation
    ggplot2::geom_text(data = flower_unique, aes(x = text_x, y = text_y, label = paste(unique, paste0("(", formatC(as.numeric(absent), format = "f", digits = 0, big.mark = ","), ")"), sep = " "), angle = text_angle), size = 5, fontface = "bold") +  # Corrected text rotation
    ggplot2::geom_text(data = flower_accessory, aes(x = text_x, y = text_y, label = formatC(as.numeric(accessory), format = "f", digits = 0, big.mark = ","), angle = text_angle), size = 5, fontface = "bold") +  # Corrected text rotation
    ggplot2::geom_text(data = data.frame(x = 0, y = 0), aes(x = x, y = y, label = paste("Core gene", formatC(sum(as.numeric(core_genes)) / num_petals, format = "f", digits = 0, big.mark = ","), sep = "\n")), size = 6, fontface = "bold") +  # Core gene count in the center
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  # Determine output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()  # Default to current directory
  }

  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save the plot to a PDF
  # Save plot if save_plot is TRUE
  if (save_plot) {
    if (is.null(output_dir)) {
      output_dir <- getwd()
    }
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    output_file <- file.path(output_dir, "flower_plot.pdf")
    ggplot2::ggsave(output_file, plot, width = 10, height = 10)
    message("Plot saved to: ", output_file)
  }

  # Return the plot
  return(plot)
}

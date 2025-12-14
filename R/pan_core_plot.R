#' Generate Pan-Core Genome Plot
#'
#' This function generates a pan-core genome plot based on provided genomic data files. The plot includes both
#' pan-genome and core-genome trends, along with associated boxplots and jitter plots for visualization.
#'
#' @param directory The directory containing the required files: `curve.xls`, `list`, `pan_genome.txt`,
#'   and `core_genome.txt`. Defaults to the current working directory if `NULL`.
#' @param save_plot Logical value indicating whether to save the plot as a PDF. Defaults to `TRUE`.
#' @param plot_filename The name of the output PDF file. If `NULL`, defaults to `"pan_core_plot.pdf"`.
#' @return A `ggplot` object representing the pan-core genome plot.
#' @export
#'
#' @examples
#' # Generate and save the plot in the current directory
#' pan_core_plot()
#'
#' # Generate the plot without saving
#' pan_core_plot(save_plot = FALSE)
#'
#' # Generate the plot and save to a specific directory
#' pan_core_plot(directory = "path/to/data", plot_filename = "custom_plot.pdf")

pan_core_plot <- function(directory = NULL, save_plot = TRUE, plot_filename = NULL) {
  # Step 1: Set default directory and filename
  if (is.null(directory)) {
    directory <- getwd()
    message("Directory set to current working directory: ", directory)
  }
  if (is.null(plot_filename)) {
    plot_filename <- "pan_core_plot.pdf"
    message("Plot filename set to default: ", plot_filename)
  }

  results_file <- file.path(directory, "Results/curve.xls")
  supporting_file <- file.path(directory, "Supporting_files/list")
  pan_genome_file <- file.path(directory, "Supporting_files/pan_genome.txt")
  core_genome_file <- file.path(directory, "Supporting_files/core_genome.txt")

  # Step 2: Read and preprocess curve data
  if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file)
  }
  data <- readr::read_tsv(results_file, show_col_types = FALSE)
  data <- tibble::column_to_rownames(data, var = "...1")
  assign("data", data, envir = .GlobalEnv)

  # Step 3: Extract equations
  equations <- data[rownames(data) == "Equation", ]
  for (col in colnames(equations)) {
    equation <- equations[[col]] %>%
      stringr::str_extract("(?<=\\=).*") %>%
      stringr::str_squish()
    assign(paste0("equation_", col), equation, envir = .GlobalEnv)
    message("Extracted for ", col, ": ", equation)
  }

  # Step 4: Preprocess parameters
  df <- data %>%
    tibble::rownames_to_column(var = "rownames") %>%
    dplyr::mutate(across(everything(), ~ stringr::str_squish(gsub("_", " ", .)))) %>%
    dplyr::mutate(across(everything(), ~ ifelse(. == "", NA, .))) %>%
    tidyr::fill(everything(), .direction = "down")

  parameters_row <- df %>% dplyr::filter(rownames == "Parameters")
  for (col in colnames(parameters_row)[-1]) {
    key_value_pairs <- parameters_row[[col]] %>%
      stringr::str_squish() %>%
      stringr::str_extract_all("[a-zA-Z]+=\\s?[-0-9.]+") %>%
      unlist() %>%
      stringr::str_split("=", simplify = TRUE)
    for (i in seq_len(nrow(key_value_pairs))) {
      key <- stringr::str_trim(key_value_pairs[i, 1])
      value <- as.numeric(key_value_pairs[i, 2])
      assign(key, value, envir = .GlobalEnv)
      message("Saved parameter: ", key, " = ", value)
    }
  }

  # Step 5: Read genome list
  if (!file.exists(supporting_file)) {
    stop("Supporting file not found: ", supporting_file)
  }
  list <- readr::read_tsv(supporting_file, col_names = c("number", "strain"))
  x <- seq(1, nrow(list), by = 1)

  # Step 6: Read and preprocess Pan genome data
  if (!file.exists(pan_genome_file)) {
    stop("File pan_genome.txt not found in Supporting_files.")
  }
  pan_genome_data <- readr::read_tsv(pan_genome_file, show_col_types = FALSE)
  pan_genome_data <- pan_genome_data %>%
    dplyr::mutate(`Pan genome` = gsub("[^0-9.]", "", `Pan genome`)) %>%
    dplyr::mutate(`Pan genome` = as.numeric(`Pan genome`)) %>%
    dplyr::filter(!is.na(`Pan genome`)) %>%
    dplyr::mutate(`Number of genomes` = as.numeric(`Number of genomes`))

  # Step 7: Read and preprocess Core genome data
  if (!file.exists(core_genome_file)) {
    stop("File core_genome.txt not found in Supporting_files.")
  }
  core_genome_data <- readr::read_tsv(core_genome_file, show_col_types = FALSE)
  core_genome_data <- core_genome_data %>%
    dplyr::mutate(`Core genome` = gsub("[^0-9.]", "", `Core genome`)) %>%
    dplyr::mutate(`Core genome` = as.numeric(`Core genome`)) %>%
    dplyr::filter(!is.na(`Core genome`)) %>%
    dplyr::mutate(`Number of genomes` = as.numeric(`Number of genomes`))

  # Step 8: Compute Pan-Genome and Core-Genome
  a <- get("a", envir = .GlobalEnv)
  b <- get("b", envir = .GlobalEnv)
  c_param <- get("c", envir = .GlobalEnv)  # c는 이름 변경
  d_param <- get("d", envir = .GlobalEnv) 
  
  pan_genome <- a * x^b
  core_genome <- c_param * exp(d_param * x)

  # Generate data frame for visualization
  curve_data <- data.frame(
    `Number of genomes` = x,
    Pan_Genome = pan_genome,
    Core_Genome = core_genome,
    check.names = FALSE
  )

  # Visualization
  plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = curve_data, ggplot2::aes(x = `Number of genomes`, y = Pan_Genome, color = "Pan Genome"), size = 1) +
    ggplot2::geom_line(data = curve_data, ggplot2::aes(x = `Number of genomes`, y = Core_Genome, color = "Core Genome"), size = 1) +
    ggplot2::geom_boxplot(data = pan_genome_data,
                          ggplot2::aes(x = as.factor(`Number of genomes`), y = `Pan genome`, color = "Pan Genome"),
                          width = 0.3, alpha = 0.3) +
    ggplot2::geom_jitter(data = pan_genome_data,
                         ggplot2::aes(x = as.factor(`Number of genomes`), y = `Pan genome`, color = "Pan Genome"),
                         width = 0.1, size = 1.5, alpha = 0.5) +
    ggplot2::geom_boxplot(data = core_genome_data,
                          ggplot2::aes(x = as.factor(`Number of genomes`), y = `Core genome`, color = "Core Genome"),
                          width = 0.3, alpha = 0.3) +
    ggplot2::geom_jitter(data = core_genome_data,
                         ggplot2::aes(x = as.factor(`Number of genomes`), y = `Core genome`, color = "Core Genome"),
                         width = 0.1, size = 1.5, alpha = 0.5) +
    ggplot2::scale_x_discrete(
      name = "No of genomes",
      breaks = as.character(seq(min(curve_data$`Number of genomes`), max(curve_data$`Number of genomes`), by = 1)),
      labels = as.character(seq(min(curve_data$`Number of genomes`), max(curve_data$`Number of genomes`), by = 1))
    ) +
    ggplot2::scale_y_continuous(
      name = "Gene count",
      labels = scales::comma  # Add comma to y-axis numbers
    ) +
    ggplot2::scale_color_manual(values = c("Pan Genome" = "#3B81A1", "Core Genome" = "#D06461")) +
    ggplot2::labs(
      title = "Pan-core plot",
      color = "Type"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(angle = 90, size = 20),
      axis.text.x = ggplot2::element_text(color = "grey20", size = 15),
      axis.text.y = ggplot2::element_text(color = "grey20", size = 15),
      plot.title = ggplot2::element_text(size = 20)
    )

  # Save the plot if required
  if (save_plot) {
    output_file <- file.path(directory, plot_filename)
    ggplot2::ggsave(output_file, plot = plot, width = 10, height = 6, dpi = 300)
    message("Plot saved as: ", output_file)
  }

  # Print the plot
  print(plot)
}

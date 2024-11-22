#' Generate Phylogeny Plots and Color Palette Visualization
#'
#' This function generates phylogeny plots from tree files in the specified directory
#' and optionally saves the plots and color palette as PDF files.
#'
#' @param directory Path to the directory containing Supporting_files. Defaults to the current working directory.
#' @param save_pdf Logical. If `TRUE`, saves the combined phylogeny plots as a PDF. Defaults to `TRUE`.
#' @param save_palette Logical. If `TRUE`, saves the color palette visualization as a PDF. Defaults to `TRUE`.
#' @return A combined plot object of the phylogeny trees.
#' @export
#'
#' @examples
#' # Generate plots and save both PDFs
#' generate_phylogeny_plots(directory = "path/to/directory", save_pdf = TRUE, save_palette = TRUE)
#'
#' # Generate plots without saving any PDFs
#' generate_phylogeny_plots(directory = "path/to/directory", save_pdf = FALSE, save_palette = FALSE)

phylogeny_plot <- function(directory = NULL, save_plot = TRUE, save_palette = TRUE) {
  # Set directory to current working directory if not provided
  if (is.null(directory)) {
    directory <- getwd()
    message("Directory set to current working directory: ", directory)
  }

  # Define supporting file path
  supporting_file <- file.path(directory, "Supporting_files/")

  # List tree files
  tree_files <- list.files(supporting_file, pattern = "\\_MOD\\.ph$", full.names = TRUE)
  if (length(tree_files) == 0) {
    stop("No tree files found in the directory: ", supporting_file)
  }

  # Load and preprocess tree files
  for (i in seq_along(tree_files)) {
    tryCatch({
      assign(tools::file_path_sans_ext(basename(tree_files[i])), ape::read.tree(file = tree_files[i]))
    }, error = function(e) {
      stop("Error reading tree file: ", tree_files[i], "\n", e$message)
    })
  }

  # Iterate through the tree files and preprocess labels
  for (i in seq_along(tree_files)) {
    obj_name <- tools::file_path_sans_ext(basename(tree_files[i]))
    if (!exists(obj_name)) {
      stop("Tree object not found for: ", obj_name)
    }

    tree_obj <- get(obj_name)
    if (is.null(tree_obj$tip.label)) {
      stop("Tree object has no 'tip.label' attribute: ", obj_name)
    }

    tree_obj$tip.label <- stringi::stri_extract_first_regex(tree_obj$tip.label, "^[^-]*")
    tree_obj$tip.label <- as.numeric(tree_obj$tip.label)
    assign(obj_name, tree_obj)
  }

  # Generate Genome_summary_df
  if (!exists("Genbank_info_combiner")) {
    stop("Function 'Genbank_info_combiner()' not found. Please ensure it is defined.")
  }

  Genbank_info_combiner() # Ensure this function is defined elsewhere
  if (!exists("Genome_summary_df")) {
    stop("Object 'Genome_summary_df' not found after running 'Genbank_info_combiner()'.")
  }

  Genome_summary_df <- dplyr::mutate(
    Genome_summary_df,
    tip.label = 1:nrow(Genome_summary_df),
    genus = ifelse(is.na(SOURCE) | trimws(SOURCE) == "", NA_character_, qdap::beg2char(SOURCE, " ")),
    species = ifelse(is.na(SOURCE) | trimws(SOURCE) == "", NA_character_, stringr::word(SOURCE, 2)),
    strain = ifelse(is.na(SOURCE) | trimws(SOURCE) == "", NA_character_,
                    trimws(stringr::str_remove(SOURCE, paste0("^", genus, " ", species))))
  )

  # Ensure required columns exist
  required_columns <- c("genus", "species", "tip.label", "strain")
  missing_columns <- setdiff(required_columns, colnames(Genome_summary_df))
  if (length(missing_columns) > 0) {
    stop("The following required columns are missing in 'Genome_summary_df': ", paste(missing_columns, collapse = ", "))
  }

  # Define background colors
  background_color <- c(
    "#ff4947", "#f67e7d", "#f69c7d", "#f6b07d", "#AF967D", "#af7a46",
    "#dda46b", "#f2c9a1", "#f8c779", "#ffc661", "#ffb532", "#f59e03", "#f5b247", "#f8c924",
    "#f3e75b", "#fff098", "#c6dda6", "#bad780", "#919f7f",
    "#95b46a", "#A3BE8C", "#80d7c6", "#4fc6d0", "#4faad0", "#5c9ce4", "#337cd6",
    "#5565ca", "#5E81AC", "#AD8CAE", "#bc97ab", "#bb84a1",
    "#cb72a1", "#cc5293", "#a04876", "#484860", "#184860", "#c0c0c0"
  )

  # Assign colors
  assign_colors <- function(df, upper_col, sub_col, color_palette) {
    generate_color_ramp <- function(base_color, n) {
      grDevices::colorRampPalette(c(base_color, grDevices::adjustcolor(base_color, alpha.f = 0.8)))(n)
    }
    unique_upper <- unique(df[[upper_col]])
    upper_colors <- setNames(sample(color_palette, length(unique_upper), replace = TRUE), unique_upper)
    sub_colors <- list()
    for (level in unique_upper) {
      sub_levels <- df[[sub_col]][df[[upper_col]] == level]
      if (length(sub_levels) > 0) {
        base_color <- upper_colors[level]
        sub_colors[[level]] <- setNames(generate_color_ramp(base_color, length(sub_levels)), sub_levels)
      }
    }
    sub_colors_flat <- unlist(sub_colors)
    names(sub_colors_flat) <- stringr::str_remove(names(sub_colors_flat), ".*\\.")
    df$Color <- sub_colors_flat[as.character(df[[sub_col]])]
    df$Color[is.na(df$Color)] <- "#c0c0c0"
    return(df)
  }
  Genome_summary_df_with_colors <- assign_colors(Genome_summary_df, "species", "strain", background_color)

  # Save the color palette visualization
  if (save_palette) {
    grDevices::pdf(file = file.path(directory, "palette_color.pdf"), width = 10, height = 5)
    unikn::seecol(unikn::newpal(
      col = as.character(Genome_summary_df_with_colors$Color),
      names = as.character(paste(
        Genome_summary_df_with_colors$genus,
        Genome_summary_df_with_colors$species,
        Genome_summary_df_with_colors$strain
      ))
    ))
    grDevices::dev.off()
    cat("Palette color PDF saved as 'palette_color.pdf' in directory: ", directory, "\n")
  } else {
    cat("Palette color visualization skipped.\n")
  }

  # Generate plots with titles
  plot_list <- list()
  tree_files_short <- tools::file_path_sans_ext(basename(tree_files))

  for (i in seq_along(tree_files_short)) {
    obj_name <- tree_files_short[i]
    if (!exists(obj_name)) {
      stop("Tree object does not exist: ", obj_name)
    }

    OO <- ggtree::ggtree(get(obj_name), layout = "rectangular", size = 0.5, linetype = 1) + ggplot2::xlim(0, 1)
    OO$data <- dplyr::mutate(
      OO$data,
      label = as.numeric(label)
    ) %>%
      dplyr::left_join(Genome_summary_df_with_colors, by = c("label" = "tip.label"))
    plot <- OO +
      ggtree::geom_tiplab(
        ggplot2::aes(label = paste(OO$data$genus, OO$data$species, OO$data$strain), color = OO$data$Color),
        fontface = "italic",
        align = TRUE,
        linesize = 0.0,
        offset = 0.01,
        show.legend = FALSE
      ) +
      ggtree::geom_tiplab(
        ggplot2::aes(color = OO$data$Color, label = ""),
        fontface = "italic",
        align = TRUE,
        linesize = 0.5,
        linetype = "dashed",
        offset = 0.0,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_identity() +
      ggplot2::theme(
        legend.position = c(0.96, 0.5),
        legend.background = ggplot2::element_rect(fill = NA),
        legend.title = ggplot2::element_text(size = 7),
        legend.text = ggplot2::element_text(size = 6),
        legend.spacing.y = grid::unit(0.02, "cm")
      )

    format_title <- function(file_name) {
      if (grepl("PAN_PHYLOGENY", file_name, ignore.case = TRUE)) {
        return("Pan phylogeny")
      } else if (grepl("CORE_PHYLOGENY", file_name, ignore.case = TRUE)) {
        return("Core phylogeny")
      } else {
        return(file_name)  # Default: return the original file name
      }
    }

    # Add formatted title as a grob
    plot_with_title <- gridExtra::arrangeGrob(
      grobs = list(ggplotGrob(plot)),
      top = textGrob(
        label = format_title(tree_files[i]),
        gp = gpar(fontsize = 12, fontface = "bold"),
        just = "left",  # Left-align the title
        x = 0.05           # Align at the far left
      )
    )
    plot_list[[i]] <- plot_with_title
  }

  # Combine plots and optionally save as PDF
  if (save_plot) {
    grDevices::pdf(file = file.path(directory, "phylogeny_plots.pdf"), width = 16, height = 8)
    do.call(gridExtra::grid.arrange, c(plot_list, ncol = 2))
    grDevices::dev.off()
    cat("PDF saved as 'phylogeny_plots.pdf' in directory: ", directory, "\n")
  } else {
    cat("PDF saving skipped. Combined plots generated but not saved.\n")
  }

  # Return the combined plot
  return(do.call(gridExtra::grid.arrange, c(plot_list, ncol = 2)))
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_annotations.R
\name{process_annotations}
\alias{process_annotations}
\title{Process Gene and COG Annotations}
\usage{
process_annotations(directory = NULL, PanDB = NULL)
}
\arguments{
\item{directory}{A character string specifying the base directory for processing. If \code{NULL}, the current working directory
(\code{getwd()}) will be used, and the function will validate the presence of required subdirectories
(\code{"Results"} and \code{"Supporting_files"}).}

\item{PanDB}{A global environment list containing the processed GenBank files. If \code{NULL}, the function assumes the object
exists in the global environment.}
}
\value{
A list containing:
\itemize{
\item \code{Gene_conc}: A data frame of concatenated gene annotations.
\item \code{COG_conc}: A data frame of concatenated COG category annotations.
}
If no \code{.annotations.} files are found, the function returns \code{NULL}.
}
\description{
This function processes \code{.annotations.} files found in a specified directory, extracts gene and COG category annotations,
and concatenates the results into structured outputs. The outputs include gene and COG annotations grouped by clusters.
}
\details{
The function:
\itemize{
\item Checks for the presence of \code{.annotations.} files within the specified directory or its subdirectories.
\item Extracts gene and COG annotations by concatenating relevant information from the provided \code{PanDB}.
\item Outputs cleaned and concatenated gene (\code{Gene_conc.}) and COG category (\code{COG_conc.}) annotations.
\item Skips processing if no \code{.annotations.} files are detected.
}

The function validates the directory structure to ensure the required subdirectories are present. If they are missing, the
function stops with an error.
}
\examples{
# Example usage:
process_annotations(directory = "path/to/directory", PanDB = PanDB)

# Automatically detect and process using the current working directory:
process_annotations()

}
\seealso{
\code{\link{Genbank_organizer_batch}}, \code{\link{combine_gb_info}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EggNOG_organizer.R
\name{EggNOG_annotations}
\alias{EggNOG_annotations}
\title{EggNOG emapper table #}
\usage{
EggNOG_annotations(EggNOG_dir = NULL, save_output = FALSE)
}
\arguments{
\item{EggNOG_dir}{A string specifying the directory where the eggNOG-mapper output files are located. If NULL, the current working directory is used.}

\item{save_output}{A logical value indicating whether to save the output to an Excel file in XLSX format. Defaults to FALSE.}

\item{df}{A data frame containing the data from eggNOG-mapper (http://eggnog-mapper.embl.de).}
}
\value{
A list containing summary information of the processed data.
}
\description{
EggNOG emapper table #
}

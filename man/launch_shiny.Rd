% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/launch-shiny.R
\name{launch_shiny}
\alias{launch_shiny}
\title{EasyCircR - launch Shiny app}
\usage{
launch_shiny(shiny_host = "0.0.0.0", shiny_port = 3838)
}
\arguments{
\item{shiny_host}{the IPv4 address that the application should listen on. Defaults is "0.0.0.0".}

\item{shiny_port}{the TCP port that the application should listen on. Defaults is \code{3838}.}
}
\description{
Launch the Shiny app to visualize the connection between circRNA-miRNA-gene.
The app includes multiple types of filters, the ability to save filtered results in csv/excel and the
possibility to further investigate reconstructed circRNAs.
}
\examples{
launch_shiny(shiny_host ="0.0.0.0", shiny_port = 3838)

}
\author{
Luca Parmigiani, Antonino Aparo, Simone Avesani
}

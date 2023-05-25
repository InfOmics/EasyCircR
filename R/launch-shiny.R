#' @title EasyCircR - launch Shiny app
#'
#' @description Launch the Shiny app to visualize the connection between circRNA-miRNA-gene.
#' The app includes multiple types of filters, the ability to save filtered results in csv/excel and the 
#' possibility to further investigate reconstructed circRNAs.
#' 
#' @author Luca Parmigiani, Antonino Aparo, Simone Avesani
#'
#' @param shiny_host the IPv4 address that the application should listen on. Defaults is "0.0.0.0".
#'
#' @param shiny_port the TCP port that the application should listen on. Defaults is \code{3838}.
#'
#' @examples 
#' launch_shiny(shiny_host ="0.0.0.0", shiny_port = 3838)
#' 
#' @importFrom shiny runApp 
#' @export
launch_shiny <- function(shiny_host ="0.0.0.0",shiny_port = 3838) {
    .update_shiny_symlinks()
    options(shiny.port = shiny_port)
    options(shiny.host = shiny_host) 
    shiny::runApp(appDir = system.file("app", package = "EasyCircR"))
}

#' @importFrom R.utils createLink
.update_shiny_symlinks <- function () {
    dir.shiny.data <- system.file("app/data",package="EasyCircR")
    dir.create(dir.shiny.data, showWarnings = FALSE)

    filesToUpdate = c(file.path("EasyCirc/circRNA/miRNA", "circRNA_sequences.txt"),
                      file.path("EasyCirc/geneMirnaCirc", "geneMirnaCirc.rds"),
                      file.path("EasyCirc/circRNA/CIRI-Full/","predictedCircRNAs.rds"),
                      file.path("EasyCirc/circRNA/CIRI-Full/","countMatrix.rds"),
                      file.path("EasyCirc/circRNA/CIRI-Full/"))

    if (dir.exists("EasyCirc")) {
        for (file in filesToUpdate) {
            unlink(file.path(dir.shiny.data, basename(file)))
            R.utils::createLink(file.path(dir.shiny.data, basename(file)), 
                                file, overwrite=TRUE)
        }
    } else {
        stop("Directory EasyCirc is not present, alting.")
    }
}

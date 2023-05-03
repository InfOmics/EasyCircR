#' @importFrom shiny runApp
#' @export
launch_shiny <- function(shiny_host ="127.0.0.1",shiny_port = 7775) {
    .update_shiny_symlinks()
    options(shiny.port = shiny_port)
    options(shiny.host = shiny_host) #192.168.1.23
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

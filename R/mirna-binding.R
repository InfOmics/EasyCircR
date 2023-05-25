#' @title EasyCircR - miRNA response element prediction
#'
#' @description Predict miRNA response elements (MREs) from the full sequence of circRNA using TargetScan. 
#' Check if \code(circ_mirna.txt) is stored in "EasyCirc/circRNA/miRNA" otherwise execute the step (see \code{force} parameter).
#'
#' @author Luca Parmigiani, Antonino Aparo, Simone Avesani
#' 
#' @param circ_df the \code{data.frame} containing circRNAs features as results of \code{EasyCircR::read_ciri_output()}.
#'
#' @param force \code{logical(1)}. If \code{FALSE} the tool will search for the already stored RDS output files, if
#' there are none it will generate them. If \code{TRUE} it will force the redo of all the steps of the function.
#' 
#' @return  a \code{data.frame} storing associations between each circRNA and one or more miRNA families.
#' 
#' @examples 
#' # read CIRI output
#' circ <- read_ciri_output(samples_file)
#' names(circ)
#  [1] "circ_df"  "circ_mtx"
#' circ_df <- circ$circ_df
#' circ_mtx <- circ$circ_mtx
#' 
#' # identify DE circRNAs
#' design <- model.matrix(...) #define your model 
#' contr <- limma::makeContrasts(...) #define your contrasts 
#' circ_de <- de_circrna(...) #define your parameters 
#' 
#' #Select only DE circRNAs
#' circ_df_de <- circ_df[circ_df$bsj_id %in% rownames(circ_de), ]
#' 
#' 
#' circ_mirna <- get_mirna_binding(circ_df_de, force=FALSE)
#'
#' @importFrom R.utils createLink
#' @export
get_mirna_binding <- function(circ_df, force=FALSE) {
    outdir <- "EasyCirc/circRNA/miRNA"
    filename <- "circRNA_sequences.txt"
    fileOutTargetScan <- file.path(outdir, "circ_mirna.txt")

    if (file.exists(fileOutTargetScan) && !force) {
        cat("\n\tOutput file of targetscan", fileOutTargetScan, "already present,\n\tReturning it\n\n")
        return(read.delim(fileOutTargetScan, sep = "\t"))
    }
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # Creating input file for targetscan
    circrna_sequences <- data.frame(id = circ_df$bsj_id, tax = c(9606), seq = circ_df$seq)
    # 9606 taxonomy id for human
    write.table(circrna_sequences, file.path(outdir, filename), sep = "\t", 
        quote = FALSE, row.names = FALSE, col.names = FALSE)

    # Create link for shiny app to access it
    pathshiny <-file.path(system.file("app/data", package = "EasyCircR"), 
        filename)
    unlink(pathshiny)
    R.utils::createLink(pathshiny, file.path(outdir, filename), overwrite=TRUE)

    .run_targetscan(file.path(outdir, filename), fileOutTargetScan)

    return(read.delim(fileOutTargetScan, sep = "\t"))
}

.run_targetscan <- function(filename, outfile) {
    targetscan <- system.file("exec", "targetscan_70.pl", package = "EasyCircR")
    mirna_to_scan <- system.file("data", "miR.txt", package = "EasyCircR")
    cmd <- paste("perl", targetscan, mirna_to_scan, filename, outfile)
    system(cmd)
}

#.run_rnahybrid <- function(hybrid.dt, file) {
#  results.dir <- "results/rnahybrid"
#  # need to rm file first, otherwise appended!
#  if(file.exists(file.path(results.dir, file))) {
#    cat("The file already exists. Aborting.\n")
#  } else {
#    
#    # First create FASTA for each hybrid read and pass to RNAhybrid in turn
#    # (Cannot pass in bulk as otherwise compares every combination!)
#    hybrid.dt[, .(id, L_sequence, R_sequence)]
#    hybrid.dt[, fastaL := paste0(">", id, "_L\n", L_sequence)]
#    hybrid.dt[, fastaR := paste0(">", id, "_R\n", R_sequence)]
#    
#    for (n in 1:nrow(hybrid.dt)) {
#      
#      writeLines(hybrid.dt$fastaL[n], file.path(results.dir, "RNAhybridL.fa"))
#      writeLines(hybrid.dt$fastaR[n], file.path(results.dir, "RNAhybridR.fa"))
#      
#      cmd <- paste0("RNAhybrid -c -s 3utr_human -m 1000 -n 1000 -q ", results.dir, "/RNAhybridL.fa -t ", results.dir, "/RNAhybridR.fa >> ", results.dir, "/", file)
#      system(cmd)
#      
#    }
#  }
#}

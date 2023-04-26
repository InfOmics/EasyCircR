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
    pathshiny <-file.path(system.file("app/data", package = "EasyCirc"), 
        filename)
    unlink(pathshiny)
    R.utils::createLink(pathshiny, file.path(outdir, filename), overwrite=TRUE)

    .run_targetscan(file.path(outdir, filename), fileOutTargetScan)

    return(read.delim(fileOutTargetScan, sep = "\t"))
}

.run_targetscan <- function(filename, outfile) {
    targetscan <- system.file("exec", "targetscan_70.pl", package = "EasyCirc")
    mirna_to_scan <- system.file("data", "miR.txt", package = "EasyCirc")
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

#--------------------------------------------------------------------------------
# TESTING
#--------------------------------------------------------------------------------
.testReadCIRIoutput <- function() {
    #setwd('/home/luca/Work/IOR/work/EasyCirc/R')
    #sampleFile <- system.file("extdata","samples_VL51.txt", package="EasyCirc")
    sampleFile <- system.file("extdata","samples_TMD8_PQR.txt", package="EasyCirc")
    genomeFile <- "/home/luca/Data/Bio/ensmbl/hg38.fa"
    genomeAnnotation <- "/home/luca/Data/Bio/ensmbl/hg38.gtf"
    trimReadsLength <- 130
    trimReads(sampleFile, trimReadsLength)
    run_ciri_full(sampleFile, genomeFile, genomeAnnotation)
    #--------------------------------------------------------------------------------

    circ <- read_ciri_output(sampleFile)
    names(circ)
    circ_df <- circ$circ_df 
    circ_mtx <- circ$circ_mtx
    head(circ_df)
    head(circ_mtx)
    nrow(circ_mtx)
    nrow(circ_df)
    head(circ)
    mirna <- get_mirna_binding(circ, force=T)
}

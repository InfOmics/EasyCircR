#' @title EasyCircR - plot circRNA structure
#' 
#' @description Show circRNA results of CIRI-vis, which is part of the CIRI-Full algorithm. 
#' The function returns a plot without launching the shiny app and checks if there are some possible 
#' isoforms for the same circRNA (Back-Splice junction ID, \code{bsj_id}.)
#' If multiple samples have the same \code{bsj_id}, the function plots the first sample in the \code{sample_file}.
#'  
#' @author Luca Parmigiani, Antonino Aparo, Simone Avesani
#' 
#' @param circ_df the \code{data.frame} containing circRNAs features as results of \code{EasyCircR::read_ciri_output()}.
#'
#' @param bsj_id the Back-Splice junction ID of the circRNA user wants to display.
#' @param samplename the name of the sample (as reported in the \code(sample_file)) for which user wants to show the circRNA.
#' 
#' @examples 
#' circ <- read_ciri_output(sampleFile)
#' circ_df <- circ$circ_df 
#' circ_mtx <- circ$circ_mtx
#' 
#' bsj_id = "1:117402186|117405645:+:2"
#' plot_circ(circ_df, bsj_id)
#' plot_circ(circ_df, bsj_id, "TMD8_DMSO_8h")
#'
#' @importFrom Biobase openPDF
#' @export
plot_circ = function (circ_df, bsj_id, samplename=NULL) {
    #samplename=NULL
    query <- circ_df[circ_df$bsj_id == bsj_id,]
    head(circ_df[grepl(",",circ_df$samples),])
    if (nrow(query) > 0) {

        if (is.null(samplename)) {
            if (any(grepl(",", query$samples))) {
                cat("Multiple samples have the same bsj_id (",query$samples,").  \nPlotting the first in the list.\nChoose a specific samplename otherwise.\n")
            }
            # Extract sample name and stout_id of circ
            preproc_stout_id = strsplit(query$stout_id, split=",")[[1]][1]
            stout_id = strsplit(preproc_stout_id, split=":")[[1]][2]
            samplename = strsplit(query$samples, split=",")[[1]][1]
            Biobase::openPDF(normalizePath(file.path("EasyCirc/circRNA/CIRI-Full", samplename, 
                                       "CIRI-vis_out", paste0(stout_id,".pdf"))))
        } else {
            if (any(grepl(samplename, query$samples))) {
                preproc_stout_id = strsplit(query$stout_id, split=",")[[1]]
                preproc_stout_id = preproc_stout_id[grepl(samplename,preproc_stout_id)]
                stout_id = strsplit(preproc_stout_id, split=":")[[1]][2]
                Biobase::openPDF(normalizePath(file.path("EasyCirc/circRNA/CIRI-Full", samplename, 
                                           "CIRI-vis_out", paste0(stout_id,".pdf"))))
            } else {
                cat("samplename:", samplename, "not found.\nPossible samplenames for",bsj_id,"\n")
                print(strsplit(query$samples,split=",")[[1]])
            }
        }

        query <- circ_df[circ_df$bsj_id == bsj_id,]
    } else {
        isoform <- circ_df[circ_df$bsj == bsj_id,]
        if (nrow(isoform) > 0) {
            cat("There are some possible isoforms for the bsj:", bsj_id,"\n")
            print(isoform$bsj_id)
        }
        stop("bsj_id: ", bsj_id, " not found\n")
    }
}

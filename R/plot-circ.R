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
    #--------------------------------------------------------------------------------
    # Test plot
    #Multiple isoform
    bsj_id = "3:158122103|158123991:+:1"
    plot_circ(circ_df, bsj_id)
    #Multiple sample
    bsj_id = "1:117402186|117405645:+:2"
    plot_circ(circ_df, bsj_id)
    plot_circ(circ_df, bsj_id, "TMD8_DMSO_8h")
    plot_circ(circ_df, bsj_id, "TMD8_DMSO_12h")
    plot_circ(circ_df, bsj_id,"TMD8_PQR_12h")
}


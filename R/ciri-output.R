#' @import dplyr
#' @importFrom magrittr %>%
#' @export
read_ciri_output <- function(sampleFile=NULL, chrMT=T) {
    ciridir <- "EasyCirc/circRNA/CIRI-Full"
    if (is.null(sampleFile)) {
        dirs <- list.dirs(full.names=FALSE, ciridir)
        dirs <- dirs[grepl("CIRI-vis_out", dirs)]
        samplenames <- dirname(dirs)
    } else {
        samples <- read_samplefile(sampleFile, showWarnings=FALSE)
        samplenames <- samples$SampleName
    }

    #--------------------------------------------------------------------------------
    # Extracting circ_df
    #--------------------------------------------------------------------------------
    cat("Extracting circRNAs data.frame ... ")
    circ.list <- vector("list")
    for (samplename in samplenames) {
        circ.list[[samplename]] <- .read_cirivis_output(samplename)
    }
    circ_df <- do.call("rbind", circ.list)
    aggregate_what <- grep("stout_id|sample", colnames(circ_df))
    circ_df <- aggregate(circ_df[,aggregate_what], circ_df[,-aggregate_what], FUN = paste, collapse = ",")
    colnames(circ_df)[ncol(circ_df)] <- "samples"

    # Number of isoform per bsj
    #iso_group_by = circ_df %>% dplyr::group_by(bsj) %>% dplyr::count(sort=TRUE)
    iso_group_by <- circ_df %>% dplyr::group_by(bsj) %>% dplyr::count() 

    # Sort by bsj and length
    circ_df <- circ_df %>% dplyr::arrange(bsj, length)

    # Create unique bsj_id, considering mulitple isoform for a single bsj
    circ_df$bsj_id <- circ_df$bsj
    bsj_isoform <- iso_group_by$bsj[iso_group_by$n > 1]

    # Append to bsj_id a ":" if there is an isoform
    circ_df[circ_df$bsj_id %in% bsj_isoform,"bsj_id"] <- sapply(circ_df[circ_df$bsj_id %in% bsj_isoform,"bsj_id"], paste0, ":")

    # Add isoform number based on length
    j = 1
    for (i in 1:nrow(circ_df)) {
        #Ends with :, so has an isoform
        if (grepl(":$", circ_df$bsj_id[i])) {
            circ_df$bsj_id[i] <- paste0(circ_df$bsj_id[i],j)
            j = j + 1
        } else {
            j = 1
        }
    }
    #sum(duplicated(circ_df$bsj_id))

    cat("OK\n")
    #--------------------------------------------------------------------------------
    # Extracting circ_mtx
    #--------------------------------------------------------------------------------
    cat("Extracting circRNAs count matrix ... ")
    circ_mtx <- .read_circ_count_mtx(samplenames, circ_df)
    cat("OK\n")
    
    if(chrMT == F){
        cat("Removing circRNAs from chromosome MT ...")
        circ_df = circ_df[circ_df$chr != "MT", ]
        circ_mtx = circ_mtx[rownames(circ_mtx) %in% circ_df$bsj_id,]
    }
    circ_val = .validate_circ(circ_df)
    circ_df = merge(circ_df, circ_val[,c(7,11)], by.x = "bsj_id", by.y = "V7", all.x = T)
    colnames(circ_df)[13] = "validated_name"
    unlink(file.path(ciridir, "predictedCircRNAs.rds"))
    #saveRDS_shiny(circ_df[, c("bsj", "seq", "length")], file = file.path(ciridir, "predictedCircRNAs.rds"))
    .saveRDS_shiny(circ_df, file = file.path(ciridir, "predictedCircRNAs.rds"))
    .saveRDS_shiny(circ_mtx, file = file.path(ciridir, "countMatrix.rds"))

    return(list(circ_df = circ_df, circ_mtx = circ_mtx))
}


.read_cirivis_output <- function (samplename) {
    circ.list <- vector("list")

    outfile <- file.path("EasyCirc/circRNA/CIRI-Full", 
                         samplename, "CIRI-vis_out/stout.list_circle.fa")
    path.stout.list <- file.path("EasyCirc/circRNA/CIRI-Full", 
                                 samplename, "CIRI-vis_out/stout.list") 

    columns <- c("Image_ID","Circle_ID","Chr","start","end",
                 "total_exp","isoform_number","isoform_exp",
                 "isoform_length","isoform_state","strain",
                 "gene_id","isoform_cirexon")

    if (!file.exists(outfile)) {
        cat("\nWarning: CIRI-Full was not able to reconstruct any circRNA from", 
            samplename, "\n")
        circ_df <- data.frame(matrix(ncol = length(columns), nrow = 0))
        colnames(circ_df) <- columns
        return(circ_df)
    } else {
        con = file(outfile, "r")
        i = 1
        while (TRUE) {
            line = readLines(con, n = 1)
            if (length(line) == 0) {
              break
            }
            if (substr(line,1,1) == ">") {
                #line <- ">stout_173#6:72748093|72748213 length=121 1/1 -"
                stout_id <- strsplit(line, split="#")[[1]][1]
                stout_id <- substr(stout_id, 2,nchar(stout_id))
                info_circ <- strsplit(strsplit(line, split="#")[[1]][2], split=" ")[[1]]
                strand <- info_circ[4]
                bsj <- paste0(info_circ[1],":",strand)
                chr <- strsplit(bsj,split=":")[[1]][1]
                start_end <- strsplit(bsj,split=":")[[1]][2]
                start <- strsplit(start_end, split="|",fixed=T )[[1]][1]
                end <- strsplit(start_end, split="|",fixed=T )[[1]][2]
                length <- as.numeric(gsub("length=","",info_circ[2]))

                circ.list[[i]] <- data.frame(bsj=bsj, chr=chr, start=start, end=end, strand=strand,
                                             seq=NA, length=length, stout_id=paste0(samplename,":",stout_id), 
                                             sample=samplename)
            } else {
                circ.list[[i]]["seq"] <- gsub("T","U",line)
                i <- i + 1
            }
        }
        close(con)

        circ.fa.df <- do.call("rbind", circ.list)
        stout.list.df <- read.csv(path.stout.list, sep="\t", header=FALSE, skip=1)
        colnames(stout.list.df) = columns
        stout.list.df$bsj  <- paste0(stout.list.df$Circle_ID,":",stout.list.df$strain)
        stout.list.df <- stout.list.df %>% dplyr::filter(isoform_state == "Full")

        circ_df <- dplyr::left_join(stout.list.df, circ.fa.df, by=c("bsj"="bsj"))
        colnames(circ_df)[grep("start.y",colnames(circ_df))] <- "start"
        colnames(circ_df)[grep("end.y",colnames(circ_df))] <- "end"
        return(circ_df[,c(colnames(circ.fa.df),"isoform_cirexon", "gene_id")])
    }
}

.read_circ_count_mtx <- function(samplenames, circ_df) {
    #--------------------------------------------------------------------------------
    # Join CIRI-Full results
    #--------------------------------------------------------------------------------
    j <- 1
    samplename <- samplenames[j]
    circ_mtx <- .read_sample_circ_count_mtx(samplename, circ_df)
    while (nrow(circ_mtx) == 0) {
        j <- j + 1
        samplename <- samplenames[j]
        circ_mtx <- .read_sample_circ_count_mtx(samplename, circ_df)
    }

    for (i in (j+1):length(samplenames)) {
        samplename <- samplenames[i]
        circ_mtx.tmp <- .read_sample_circ_count_mtx(samplename, circ_df)
        if (nrow(circ_mtx.tmp)>0) {
            circ_mtx <- dplyr::full_join(circ_mtx, circ_mtx.tmp, by=c("bsj_id"="bsj_id"))
        }
    }
    #circ_mtx <- aggregate(. ~ bsj_id, circ_mtx, sum)
    circ_mtx[is.na(circ_mtx)] <- 0
    rownames(circ_mtx) <- circ_mtx$bsj_id
    circ_mtx$bsj_id <- NULL
    head(circ_mtx)

    return(circ_mtx)
}

.read_sample_circ_count_mtx = function (samplename, circ_df) {

    path <- "EasyCirc/circRNA/CIRI-Full"

    stout.list.output <- file.path(path, samplename,"CIRI-vis_out/stout.list")

    columns <- c("Image_ID","Circle_ID","Chr","start","end",
                 "total_exp","isoform_number","isoform_exp",
                 "isoform_length","isoform_state","strain",
                 "gene_id","isoform_cirexon")

    info = file.info(stout.list.output)

    if (info$size != 0) {
        ciri.output <- read.table(stout.list.output, header = FALSE, 
                                  sep = "\t", row.names=NULL, skip=1)
        colnames(ciri.output) <- columns

        # Filter by Full reconstruct circRNAs
        ciri.output <- ciri.output[ciri.output$isoform_state=="Full",]
        #Create BSJ_id adding the strain
        ciri.output$bsj <- paste0(ciri.output$Circle_ID,":",ciri.output$strain)
        ciri.output <- dplyr::left_join(ciri.output, 
                                        circ_df[,c("bsj_id", "bsj", "isoform_cirexon")], 
                                        by=c("bsj"="bsj", 
                                             "isoform_cirexon" = "isoform_cirexon"))
        circ_mtx <- data.frame(bsj_id = ciri.output$bsj_id, x = ciri.output$total_exp)
        colnames(circ_mtx) <- c("bsj_id", samplename)
        return(circ_mtx)
    } else {
        cat("\nWarning: Omitting", samplename, "from count matrix\n")
        circ_mtx <- data.frame(matrix(ncol = 2, nrow = 0))
        colnames(circ_mtx) <- c("bsj_id", samplename)
        return(circ_mtx)
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
    #trimReads(sampleFile, trimReadsLength)
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
}

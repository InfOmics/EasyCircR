#' @importFrom parallel detectCores
#' @export
run_ciri_full <- function(samples_file, genome_file, genome_annotation_file, 
                          percent_of_reads_after_trim=NULL, 
                          trim_reads_length=NULL, force=FALSE, 
                          n_core=parallel::detectCores(), remove_temporary_files=FALSE) {
    samples <- read_samplefile(samples_file)

    trimdir <- "EasyCirc/circRNA/trimmed_fastq"
    indir <- "EasyCirc/circRNA/CIRI-Full/"
    outdir <- "EasyCirc/circRNA"

    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    #--------------------------------------------------------------------------------
    # Run CIRI-Full for each samples
    #--------------------------------------------------------------------------------
    for (i in 1:nrow(samples)) {
        samplename <- samples$SampleName[i]

        fq1 <- samples$FileName1[i]
        fq2 <- samples$FileName2[i]

        if (force && dir.exists(file.path(outdir,"CIRI-Full",samplename))) {
            unlink(file.path(outdir,"CIRI-Full",samplename),recursive=TRUE)
        }

        # Trim Reads
        if (!is.null(trim_reads_length)) {
            cat("--------------------------------------------------------------------------------")
            cat("\nTrimming", samplename, "\n")
            .trim_reads(fq1=fq1,
                      fq2=fq2,
                      trim_reads_length=trim_reads_length, 
                      outdir=file.path(outdir, "trimmed_fastq"),
                      samplename=samplename,
                      force=force,
                      n_core=n_core)

            fq1 <- file.path(trimdir, paste0(samplename,"_R1.trimmed.fq.gz"))
            fq2 <- file.path(trimdir, paste0(samplename,"_R2.trimmed.fq.gz"))
        } 

        if (!dir.exists(file.path(outdir,"CIRI-Full",samplename))) {
            # Run CIRI-Fulk
            cat("--------------------------------------------------------------------------------")
            cat("\nRunning CIRI-Full:", samplename, "\n")
            .run_ciri_full(fq1=fq1,
                         fq2=fq2,
                         outdir=file.path(outdir, "CIRI-Full"),
                         samplename=samplename,
                         genome_file=genome_file,
                         genome_annotation_file=genome_annotation_file,
                         force=force)

            # Clean CIRI-Full folders
            if (remove_temporary_files) {
                clean_easycirc_folder(samples_file)
            }
        } else {
            cat(samplename, "has already a CIRI-Full output directory ... skipping\n")
        }

        if (!dir.exists(file.path(outdir,"CIRI-Full",samplename,"CIRI-vis_out"))) {
            # Run CIRI-Vis
            cat("--------------------------------------------------------------------------------")
            cat("\nRunning CIRI-Vis:", samplename, "\n")
            .run_ciri_vis(samplename, genome_file, indir, force)

        } else {
            cat(samplename, "has already a CIRI-vis_out directory ... skipping\n")
        }

    }
    
}

.run_ciri_full <- function(fq1, fq2, outdir, samplename, genome_file, genome_annotation_file, force = FALSE) {
    cat("--------------------------------------------------------------------------------")
    cat("\nCIRI-Full: ", samplename, "\n")
    cat("--------------------------------------------------------------------------------")

    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    cirifull_jar <- system.file("java","CIRI-full.jar", package="EasyCirc")
    command <- paste("java -jar", cirifull_jar, "Pipeline",
                     "-1", fq1, "-2", fq2, 
                     "-a", genome_annotation_file, "-r", genome_file, 
                     "-d", file.path(outdir,samplename), "-o", samplename, "-0")

    cat(command,"\n")
    system(command)
}

.run_ciri_vis <- function (samplename, genome_file, indir, force=FALSE) {
    cat("--------------------------------------------------------------------------------")
    cat("\nCIRI-Vis: ", samplename, "\n")
    cat("--------------------------------------------------------------------------------")

    cirivis_jar <- system.file("java","CIRI-vis.jar", package="EasyCirc")
    outdir <- file.path("EasyCirc/circRNA/CIRI-Full/", samplename, "CIRI-vis_out")

    if (force && dir.exists(outdir) ) {
        unlink(outdir,recursive=TRUE)
    }

    command <- paste("java -jar", cirivis_jar,
                     "-i", file.path(indir, samplename, "CIRI-full_output", 
                                     paste0(samplename,"_merge_circRNA_detail.anno")),
                     "-l", file.path(indir, samplename, "CIRI-AS_output", 
                                     paste0(samplename,"_library_length.list")),
                     "-r", genome_file, 
                     "-d", outdir,
                     "-min 1")
    Sys.unsetenv("DISPLAY") #debug for Rstudio Server java exception "cannot connect to X11 server using ":0" as DISPLAY variable  
    cat(command,"\n")
    system(command)
}

#' @importFrom tools file_path_sans_ext
.trim_reads <- function (fq1, fq2, trim_reads_length, outdir, 
                       samplename=basename(tools::file_path_sans_ext(fq1)),
                       force=FALSE,
                       n_core=parallel::detectCores()) {

    cat("--------------------------------------------------------------------------------")
    cat("\nTrimming:", samplename, "\n")
    cat("--------------------------------------------------------------------------------")

    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    fq1_paired   <- file.path(outdir, paste0(samplename,"_R1.trimmed.fq.gz"))
    fq1_unpaired <- file.path(outdir, paste0(samplename,"_R1un.trimmed.fq.gz"))
    fq2_paired   <- file.path(outdir, paste0(samplename,"_R2.trimmed.fq.gz"))
    fq2_unpaired <- file.path(outdir, paste0(samplename,"_R2un.trimmed.fq.gz"))
    if (!force && file.exists(fq1_paired) 
               && file.exists(fq1_unpaired)
               && file.exists(fq2_paired)
               && file.exists(fq2_unpaired) ) {
        cat("Already trimmed (use force=TRUE if you want to trim with different length) ... OK\n")
    } else {
        trimmomatic_jar <- system.file("java","trimmomatic-0.39.jar", package="EasyCirc")

        cutfilter = paste0("CROP:",trim_reads_length," MINLEN:",trim_reads_length)

        cat("Running trimmomatic ... \n")
        command <- paste("java -jar",trimmomatic_jar,
                         "PE -phred33 -threads",n_core,
                         fq1,fq2,
                         fq1_paired,fq1_unpaired,
                         fq2_paired,fq2_unpaired,
                         cutfilter)
        cat(command,"\n")
        system(command)
        cat("\nOK\n")
        #cat("Compressing ... ")
        #R.utils::gzip(fq1_paired)
        #R.utils::gzip(fq2_paired)
        #R.utils::gzip(fq1_unpaired)
        #R.utils::gzip(fq2_unpaired)
        #cat("OK\n")
    }
}

clean_easycirc_folder = function (samples_file) {
    samples <- read_samplefile(samples_file)
    #--------------------------------------------------------------------------------
    # Clean all the CIRI-Full paths
    #--------------------------------------------------------------------------------
    for (i in 1:nrow(samples)) {
        samplename <- samples$SampleName[i]
        .clean_easycirc_folder(samplename)
    }
}

.clean_easycirc_folder = function (samplename) {
    folder <- "EasyCirc/circRNA"

    ciri_full_intermediate_folder <- c("CIRI-AS_output", "CIRI-full_output",
                                        "CIRI_output", "sam")

    path_sample <- file.path(folder,"CIRI-Full",samplename)
    if (dir.exists(path_sample)) {
        cat("Removing intermediate files from ", path_sample, " folder\n")
        for (ciri_folder in ciri_full_intermediate_folder) {
            unlink(file.path(path_sample, ciri_folder), recursive=TRUE)
        }
    }
}

#--------------------------------------------------------------------------------
# Testing CIRI-FULL
#--------------------------------------------------------------------------------
.testCIRIFull <- function() {
    #setwd('/home/luca/Work/IOR/work/EasyCirc/R')
    #samples_file <- system.file("extdata","samples_VL51.txt", package="EasyCirc")
    samples_file <- system.file("extdata","samples_TMD8_PQR.txt", package="EasyCirc")
    genome_file <- "/home/luca/Data/Bio/ensmbl/hg38.fa"
    genome_annotation_file <- "/home/luca/Data/Bio/ensmbl/hg38.gtf"
    trim_reads_length <- 130
    #trim_reads(samples_file, trim_reads_length)
    run_ciri_full(samples_file, genome_file, genome_annotation_file)
}

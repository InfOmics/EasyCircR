.get_reads_length <- function (fastq, reads_length=rep(0,1000)) {
    con = file(fastq, "r")

    i = 3
    N = length(reads_length)
    while ( TRUE ) {
        line = readLines(con, n = 1)
        if ( length(line) == 0 ) {
          break
        }
        if (i == 0) {
            #print(line)
            len = nchar(line)
            if (len > N) {
                reads_length = c(reads_length, rep(0, length(reads_length)-len))
            }
            reads_length[len] = reads_length[len] + 1
        }
        i = (i + 1) %% 4
    }

    close(con)

    return(reads_length)
}
#' @importFrom R.utils decompressFile
#' @import ggplot2
#' @export
plot_readlen_hist = function (samples_file, threshold_read_length=0.8, fast=TRUE) {

    samples <- read_samplefile(samples_file)

    cat("Reading fastq ... ")
    N = nrow(samples)
    # If fast read only first sample
    if (fast)
        N = 1
    for (i in 1:N) {
        fq1gz <- samples$FileName1[i]
        fq1=gsub(sprintf("[.]%s$", "gz"), "", fq1gz, ignore.case = TRUE)
        if (!file.exists(fq1)) 
            R.utils::decompressFile(fq1gz, destname=fq1, remove=F, ext="gz", FUN=gzfile)
        len = .get_reads_length(fq1)
        file.remove(fq1)

        fq2gz <- samples$FileName2[i]
        fq2=gsub(sprintf("[.]%s$", "gz"), "", fq2gz, ignore.case = TRUE)
        if (!file.exists(fq1)) 
            R.utils::decompressFile(fq2gz, destname=fq2, remove=F, ext="gz", FUN=gzfile)
        len = .get_reads_length(fq2, len)
        file.remove(fq2)
    }
    cat("OK\n")
    max_length = 0
    for (i in length(len):2) {
        if (len[i] != 0) 
            break
    }
    max_length=i
    len = len[1:max_length]
    cumsum(len)
    #print(len)

    if (sum(len==0) == max_length - 1) {
        cat("All the reads have the same length:", which(len!=0), "bp\n")
    } 
    else {
        tot_reads <- length(len)
        tot <- 0
        for (i in length(len):1) {
            tot <- tot + len[i]
            if (tot/tot_reads >= threshold_read_length)
                break
        }
        vertical_line = i

        #df = data.frame(length_reads=max_length:1, n =rev(len) )
        df = data.frame(length_reads=1:max_length, n =len )
        p <- ggplot2::ggplot(df,  ggplot2::aes(x=length_reads,y=n)) + 
             ggplot2::geom_bar(stat='identity', color="black", fill="white") +
             ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5)) +
             ggplot2::scale_x_continuous(name=df$x, breaks = 1:nrow(df)) + 
             ggplot2::geom_vline(ggplot2::aes(xintercept=max_length), 
                                 color="blue", linetype="dashed", size=1) 
             #ggplot2::scale_x_continuous(1:max_length)
             #ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
             #ggplot2::theme(axis.text.x= ggplot2::element_text(df$length_reads))
        p
    }
}

#get_readlen_to_trim <- function () {
#}

#--------------------------------------------------------------------------------
# Testing
#--------------------------------------------------------------------------------
.test_plot_readlen_hist <- function () {
    #samples_file <- system.file("inst/extdata","samples_VL51.txt", package="EasyCirc")
    samples_file <- "/home/luca/@focus/_playground/datatest/samples_test.txt"
    genome_file <- "/home/luca/@focus/_playground/datatest/genome/hg38_chr3_chr7.fa"
    genome_annotation_file <- "/home/luca/@focus/_playground/datatest/genome/hg38_chr3_chr7.gtf"
    plot_readlen_hist(samples_file)
}


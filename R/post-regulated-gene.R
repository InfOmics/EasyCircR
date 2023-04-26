#' @param samples_file path to tab separated file with three named columns indicating:
#' two pair-end fastq files and sample name.  Column names are FileName1,
#' FileName2 and SampleName.  
#' @param genome_file path to genome file.
#' @param genome_annotation_file path to annotation file.
#' @param condition \code{numeric}, \code{character} or \code{factor} with two levels
#'   that groups the samples (columns of \code{cntEx} and \code{cntIn}) into two
#'   conditions. The contrast will be defined as secondLevel - firstLevel.
#' @param outdir 
#' @param aligner selects the aligner program to be used for aligning the reads for the eisa step. Currently, only “Rbowtie” and “Rhisat2” (default) are supported.
#' @param force \code{logical(1)}. If \code{FALSE} it will search for the already stored RDS output files, if
#' there are none it will generate them. If \code{TRUE} it will force the redo of all 
#' the steps of the pipeline.
#' @param stranded_data logical(1). If TRUE, the RNA-seq data is assumed to be strand-specific, and therefore only overlapping genes that are on the same strand will be filtered out during the eisa step. If FALSE, also genes overlapping on opposite strands will be filtered out.
#' @param n_core number of core to be used. Default is 1.
#' @param cacheDir 
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom eisaR getRegionsFromTxDb runEISA
#' @importFrom parallel makeCluster
#' @importFrom QuasR qAlign qCount
#' @import dplyr
#' @export
get_postregulated_gene = function (samples_file, genome_file, genome_annotation_file, condition,
                                   outdir='.', aligner="Rhisat2", force = FALSE, 
                                   stranded_data=TRUE, n_core=1,
                                   cacheDir=NULL) {

    output <- file.path(outdir,"EasyCirc/postGene","postGene.rds")
    if (force | !file.exists(output)) {
        txdb <- GenomicFeatures::makeTxDbFromGFF(genome_annotation_file)
        regions <- eisaR::getRegionsFromTxDb(txdb = txdb, strandedData = stranded_data)
        
        # Align
        cl <- parallel::makeCluster(n_core)
        if (!is.null(cacheDir))
            dir.create(cacheDir, showWarnings=F)
        proj <- QuasR::qAlign(sampleFile  = samples_file,
                              genome = genome_file,
                              aligner = aligner, splicedAlignment = TRUE,
                              cacheDir = cacheDir, clObj=cl) 

        #Count alignments in exons and gene bodies
        cntEx <- QuasR::qCount(proj, regions$exons, orientation = "any")
        cntGb <- QuasR::qCount(proj, regions$genebodies, orientation = "any")
        cntIn <- cntGb - cntEx

        Rex <- cntEx[, colnames(cntEx) != "width"]
        Rin <- cntIn[, colnames(cntIn) != "width"]

        # run EISA
        res <- eisaR::runEISA(Rex, Rin, condition)
        ExIn <- res$tab.ExIn
        ExIn <- ExIn[order(ExIn$PValue),]
        head(ExIn)

        dir.create(file.path(outdir,"EasyCirc/postGene"), showWarnings = FALSE , 
                   recursive=TRUE)
        saveRDS(ExIn, file = output) 
        unlink(cacheDir, recursive = TRUE)
    } else {
        cat("Post regulated genes already found, returning them\n")
    }
    ExIn <- readRDS(output)
    return(ExIn)
}


#' @importFrom msigdbr msigdbr
.addGeneInfo <- function (gene_mirna_circ, category, subcategory=NULL) {
    extra_info <- msigdbr::msigdbr(species = "Homo sapiens", 
                                  category=category, 
                                  subcategory=subcategory)

    extra_info <- extra_info %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
    colnames(extra_info) <- c(paste0("category_",category), "target_symbol")
    return(dplyr::left_join(gene_mirna_circ, extra_info, by="target_symbol"))
}

#--------------------------------------------------------------------------------
# TESTING
#--------------------------------------------------------------------------------

.testGENE <- function () {
    samples_file <- system.file("extdata","samples_rna_paired.txt", package="EasyCirc")
    genome_file <- "/home/luca/Data/Bio/ensmbl/hg38.fa"
    genome_annotation_file <- "/home/luca/Data/Bio/ensmbl/hg38.gtf"
    # create condition factor (contrast will be TN - ES)
    condition <- factor(c("ES", "ES", "TN", "TN"))
    levels(condition)
    aligner <- "Rhisat2" #or Rbowtie

    post_gene <- get_postregulated_gene(samples_file, genome_file, condition, genome_annotation_file, aligner=aligner, force=TRUE)
    head(post_gene)
    nrow(post_gene)
    outdir='.'
    stranded_data=TRUE
    force=TRUE
}

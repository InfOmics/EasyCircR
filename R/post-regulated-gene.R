#' @title EasyCircR - Quantify post-transcriptional regulation of gene expression
#'
#' @description  Measure changes across different experimental conditions to quantify post-transcriptional regulation of gene expression (EISA).
#'  
#' @author Luca Parmigiani, Antonino Aparo, Simone Avesani
#' 
#' @param samples_file path to tab separated file with three named columns indicating:
#' two pair-end fastq files and sample name.  Column names are FileName1,
#' FileName2 and SampleName. A number of samples greater than 2 is required for each condition.  
#' @param genome_file path to genome file.
#' @param genome_annotation_file path to annotation file.
#' @param condition \code{numeric}, \code{character} or \code{factor} with two levels
#'   that groups the samples into two conditions. The contrast will be defined as secondLevel - firstLevel.
#' @param outdir path to the output directory 'EasyCirc' that store the final results. The three major folders ("circRNA", 
#' "geneMirnaCirc", "postGene") contain each one or two `.rds` file storing the results of the pipeline, 
#' e.g. `postGene.rds` the results of `get_postregulated_gene` or `countMatrix.rds` contains the counts of circRNA found.
#' @param aligner selects the aligner program to be used for aligning the reads for the EISA step. Currently, only "Rbowtie" and 
#' "Rhisat2" (default) are supported.
#' @param force \code{logical(1)}. If \code{FALSE} the tool will search for the already stored RDS output files, if
#' there are none it will generate them. If \code{TRUE} it will force the redo of all the steps of the function.
#' @param stranded_data \code{logical(1)}. If \code{TRUE}, the RNA-seq data is assumed to be strand-specific, and therefore only overlapping genes that are on the same strand will be filtered out during the EISA analysis If \code{FALSE}, also genes overlapping on opposite strands will be filtered out.
#' @param n_core number of cores to be used. Default is \code{1}.
#' @param cacheDir specifies the location to store (potentially huge) aligner temporary files . If set to \code{NULL} (default), the temporary directory of the current R session as returned by \code{tempdir()} will be used.
#' 
#' @return a \code{data.frame} with elements that stores statisical results for differential changes 
#' between exonic and intronic contrast, an indication for post-transcriptional regulation.
#' 
#' @examples
#' samples_file <- "samples.txt"
#' genome_file <- "genome/hg38.fa"
#' genome_annotation_file <- "genome/hg38.gtf"
#' condition <- factor(rep(c("DMSO","PQR"),2))
#' 
#' post_gene <- get_postregulated_gene(samples_file, genome_file, genome_annotation_file, 
#'                                    condition, aligner="Rbowtie",
#'                                    force=FALSE, n_core=10)
#'                                    
#' post_gene <- post_gene[post_gene$FDR <= 0.1,] 
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

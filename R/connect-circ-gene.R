#' @importFrom multiMiR multimir_dbInfoVersions multimir_switchDBVersion get_multimir
#' @export
connect_circ_gene <- function (circMirnas=NULL, postGenes=NULL, tabletype="validated", 
                             only_significant_genes=TRUE, predicted.cutoff=35, force=FALSE) {
    outdir <- "EasyCirc/geneMirnaCirc"
    fileGeneMirnaCirc <- file.path(outdir, "geneMirnaCirc.rds")
    mirna_family_path <- system.file("data", "miR_human.txt", package = "EasyCircR")
    mirna_family <- read.csv(mirna_family_path, sep="\t")
    mirna_family <- mirna_family[,c("MiRBase.ID","miR.family")]
    category_extra <- "H"
    subcategory_extra <- ""

    # Check if output file is already present
    if (file.exists(fileGeneMirnaCirc) && !force) {
        cat("\n\tOutput file", fileGeneMirnaCirc, "already present,\n\tReturning it\n\n")
        return(readRDS(fileGeneMirnaCirc))
    }
    if (is.null(circMirnas) || is.null(postGenes)) {
        cat("\n\tOutput file", fileGeneMirnaCirc, "not present,\n\tmissing arguments: circMirnas, postGenes\n\n")
        return(NULL)
    }
    #predicted.cutoff=35
    #tabletype="validated"
    if (only_significant_genes) {
        postGenes <- postGenes[postGenes$FDR <= 0.1,]
    }

    # Setup multiMiR database
    vers_table <- multiMiR::multimir_dbInfoVersions()
    curr_vers  <- vers_table[1, "VERSION"]  # current version
    multiMiR::multimir_switchDBVersion(db_version = curr_vers)

    #example1 <- get_multimir(mirna = 'hsa-miR-18a-3p', summary = TRUE)
    circMirnas$mirna_family <- circMirnas$miRNA_family_ID
    circMirnas$miRNA_family_ID <- NULL
    
    colnames(circMirnas)[colnames(circMirnas) == 'a_Gene_ID'] <- 'circRNA_id'
    colnames(mirna_family) <- c("mature_mirna_id", "mirna_family")
    circMirnas <- merge(circMirnas, mirna_family, by=c("mirna_family"))
    length(circMirnas$mature_mirna_id %>% unique())
    #length(circMirnas$mirna_family %>% unique())
    #library(dplyr)
        
    query <- multiMiR::get_multimir(org = "hsa",
                                    #mirna   = 'hsa-miR-18a-3p',
                                    mirna = (circMirnas$mature_mirna_id %>% unique()),
                                    table = tabletype,
                                    predicted.cutoff = predicted.cutoff,
                                    predicted.cutoff.type = "p",
                                    predicted.site = "all")
    mirnaFrame <- query@data
    mirnaFrame <- mirnaFrame[!duplicated(mirnaFrame),]

    #TODO: Accept ensmbl, entrezid and target_symbol
    postGenes$target_ensembl <- rownames(postGenes)
    print(dim(postGenes))
    geneMirna <- merge(postGenes, mirnaFrame, by="target_ensembl")
    print(geneMirna)
    head(geneMirna)
    #geneMirna[is.na(geneMirna$mature_mirna_id),] # TRUE

    head(circMirnas)
    gene_mirna_circ <- merge(geneMirna, circMirnas, by=c("mature_mirna_id"))
    head(gene_mirna_circ)
    nrow(gene_mirna_circ)

    gene_mirna_circ$chr <- sapply(gene_mirna_circ$circRNA_id, function(x){
                                  unlist(strsplit(x, split = ":"))[1]
                                 })
    
    gene_mirna_circ$circ_pos_start <- sapply(gene_mirna_circ$circRNA_id, function(x){
                                               fpart <- unlist(strsplit(x, split = "[|]"))[1]
                                               unlist(strsplit(fpart, split = ":"))[2]
                                            })
    
    gene_mirna_circ$circ_pos_end <- sapply(gene_mirna_circ$circRNA_id, function(x){
                                            epart <- unlist(strsplit(x, split = "[|]"))[2]
                                            unlist(strsplit(epart, split = ":"))[1]
                                           })
  
    # Change some columns name
    colnames(gene_mirna_circ)[colnames(gene_mirna_circ)=='UTR_start'] <- "Binding_start"
    colnames(gene_mirna_circ)[colnames(gene_mirna_circ)=='UTR_end'] <- "Binding_end"

    # Add Hallmarks
    gene_mirna_circ <- .addGeneInfo(gene_mirna_circ, category=category_extra)

    returncolumns <- c("target_symbol","mature_mirna_id", "mirna_family", "type", "circRNA_id",
                      "logFC","PValue","database", "support_type","Binding_start",
                      "Binding_end","Site_type", paste0("category_",category_extra), "target_ensembl", "chr", "circ_pos_start", "circ_pos_end")
    #adding pubmed_id and experiment remove the duplicates
    gene_mirna_circ <- gene_mirna_circ[!duplicated(gene_mirna_circ[,returncolumns]),]
    nrow(gene_mirna_circ)
    
    
    # Save output
    cat("Saving output to ", fileGeneMirnaCirc, "...")
    dir.create(outdir, showWarnings = FALSE , recursive=TRUE)
    .saveRDS_shiny(gene_mirna_circ[,returncolumns], file = fileGeneMirnaCirc) 
    cat("done\n")

    return(gene_mirna_circ[,returncolumns])
}

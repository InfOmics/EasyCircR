#Save X as an RDS file and create a soft link in app/data
#in order for the shiny app to access it
#' @importFrom R.utils createLink
.saveRDS_shiny <- function (X, file) {
    
    dir.shiny.data <- system.file("app/data",package="EasyCirc")
    dir.create(dir.shiny.data, showWarnings = FALSE)
    linkpath <- file.path(dir.shiny.data, basename(file))

    if (!length(grep(".rds$", file))) file <- paste0(file,".rds")
    unlink(linkpath, force=T)

    saveRDS(X, file=file)
    # Create link for shiny app to access it
    R.utils::createLink(linkpath, file, overwrite=TRUE)
}

#validate circRNA
#' @importFrom bedr check.binary bedr
.validate_circ <- function(circ){
  reference <- system.file("data/CircReferenceDb_hg38.bed",package="EasyCirc")
  outdir <- "EasyCirc/circRNA/CIRI-Full"
  testcirc <- file.path(outdir, "ciri_results.bed")
  if(!file.exists(reference)){
    cat("No circRNAs reference file found in package")
    return(NULL)
  }
  if(is.null(circ)){
    cat("No circRNA found in input")
    return(NULL)
  }
  if(!all(c("chr","start","end","stout_id","strand") %in% names(circ))){
    cat("Provide a valid dataframe in input\n(colnames \"chr\", \"start\", \"end\", \"stout_id\" and \"strand\" are mandatory)")
    return(NULL)
  }
  circ$score = 0
  write.table(circ[,c("chr","start","end","stout_id","score","strand", "bsj_id")], testcirc, sep = "\t", col.names = F, row.names = F, quote = F)
  
  if (bedr::check.binary("bedtools")){
    ciri_intersect <- bedr::bedr(
      input = list(a = testcirc, b = reference), 
      method = "intersect", 
      params = "-wo -f 0.99 -r -s"
    );
    ciri_results = do.call(rbind, by(ciri_intersect, ciri_intersect$V4, function(x) x[which.min(x$V14), ] ))
    return(ciri_results)
  }
  else{
    cat("You need to install bedtools before proceeding")
    return(NULL)
  }
}

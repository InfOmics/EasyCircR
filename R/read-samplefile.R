#' @importFrom tools file_path_as_absolute file_ext
#' @export
read_samplefile = function (sampleFile, showWarnings=TRUE) {

    if(!file.exists(sampleFile)){stop("File does not exists",sampleFile,call.=FALSE)}

    samples <- read.delim(sampleFile,header=TRUE,colClasses="character")

    if(nrow(samples)<2){stop(sampleFile," is either empty or there is a header missing or there are less than 2 samples",call.=FALSE)}

    # get absolute path of the sample.txt file
    sampleFileAbsolutePath <- dirname(tools::file_path_as_absolute(sampleFile))

    # ONLY PAIR END ARE ALLOWED
    # If only two columns are present sampleFile they must direct to bam files.
    if(ncol(samples)==2){
    if(all(colnames(samples) == c("FileName","SampleName"))){
      # this is a valid single end samples file. check if listed files exist
      for(i in 1:nrow(samples)){
        pathRet <- .pathAsAbsoluteRedirected(samples$FileName[i],sampleFileAbsolutePath)
            if(!is.null(pathRet)){
                if(samples$SampleName[i]==""){stop(samples$FileName[i]," listed in ",sampleFile," has no sample name",call.=FALSE)}
            samples$FileName[i] <- pathRet
        }else{
            #stop(samples$FileName[i]," listed in ",sampleFile," does not exist",call.=FALSE)
            if (showWarnings) {cat("Warning:", samples$FileName[i]," listed in ",sampleFile," does not exist\n")}        
         }
      }
      # deterime format of the files (fa,fasta,fna)
      samplesFormat <- .determineSamplesFormat(samples$FileName)

      if(samplesFormat == "bam"){
        # samples are provided as bam files
        #if(is.null(paired)){stop("Paired status of the provided bam files cannot be determined automatically, please set the paired parameter",call.=FALSE)}
      }else{
        # samples are provided as reads
        if(length(unique(samples$FileName)) != nrow(samples)){
          stop("There are duplicate files in sampleFile. This is not allowed as it would result in non-unique alignment file names",call.=FALSE)
        }

      }
    } else {stop(sampleFile," should contain column names FileName,SampleName (single end data)",call.=FALSE)} # incorrect format
    }
    # Only pair-reads are allowed for CIRIquant
    else if(ncol(samples)==3){
    if(all(colnames(samples) == c("FileName1","FileName2","SampleName"))){
        # this is a valid paired end samples file. check if listed files exist
        for(i in 1:nrow(samples)){
            pathRet <- .pathAsAbsoluteRedirected(samples$FileName1[i],sampleFileAbsolutePath)
            if(!is.null(pathRet)){
                if(samples$SampleName[i]==""){
                    stop(samples$FileName1[i]," listed in ",sampleFile," has no sample name",call.=FALSE)
                }
            samples$FileName1[i] <- pathRet
        }else{
            #stop(samples$FileName1[i]," listed in ",sampleFile," does not exist",call.=FALSE)
            if (showWarnings) {cat("Warning", samples$FileName1[i]," listed in ",sampleFile," does not exist\n")}
         }
        }
        for(i in 1:nrow(samples)){
            pathRet <- .pathAsAbsoluteRedirected(samples$FileName2[i],sampleFileAbsolutePath)
            if(!is.null(pathRet)){
                if(samples$SampleName[i]==""){
                    stop(samples$FileName2[i]," listed in ",sampleFile," has no sample name",call.=FALSE)
                }
            samples$FileName2[i] <- pathRet
        }else{
            #stop(samples$FileName2[i]," listed in ",sampleFile," does not exist",call.=FALSE)
            if (showWarnings) {cat("Warning", samples$FileName2[i]," listed in ",sampleFile," does not exist\n")}
         }
      }
      # deterime format of the files (fa,fasta,fna)
      samplesFormat <- .determineSamplesFormat(c(samples$FileName1,samples$FileName2))

      if(samplesFormat == "bam"){stop("Bam files need to be listed in a two column file: ",sampleFile,call.=FALSE)}
      if(length(unique(samples$FileName1)) != nrow(samples)){
        stop("There are duplicate files in sampleFile. This is not allowed as it would result in non-unique alignment file names",call.=FALSE)
      }

    }else{stop(sampleFile," should contain the column names FileName1,FileName2,SampleName (paired end data)",call.=FALSE)} # incorrect format
    # The format of the file is not supported
    }else{stop(sampleFile," should be a tab delimited file with either 2 or 3 columns",call.=FALSE)}

    return(samples)
}

# helper function that converts filepaths to absolute filepaths in the case where the working
# directory needs to be temporarily redirected. This is needed if e.g. the samples.txt file
# is not located in the working directory.
# returns the absolute path if the file exists
# returns NULL if the file does not exist
.pathAsAbsoluteRedirected <- function(fileName,redirectPath){
  curwd <- getwd() # store the original working directory required to jump back
  on.exit(setwd(curwd)) # make sure that it will jump back in a case of error or not

  setwd(redirectPath) # jump to the redirected position
  if(file.exists(fileName)){
    return(tools::file_path_as_absolute(fileName))
  }else{return(NULL);}
}

# given a vector of filenames, the function determines the final format 
# taking into account all the different types of extensions for 
# sequence files. Throughs an error if one ore more samples do not contain 
# a valid file extension
.determineSamplesFormat <- function(filename){
  fileExtension <- .consolidateFileExtensions(filename,compressed=TRUE) 

  validExtsSel <- fileExtension %in% c("fasta","fastq","bam","csfasta","csfastq")

  if(all(validExtsSel)){
    if(length(unique(fileExtension))==1){
      return(unique(fileExtension))
    }else{stop("Mixed sample file formats are not allowed: ",paste(unique(fileExtension),collapse=","),call.=FALSE)}
  }else{
    stop(filename[!validExtsSel][1]," does not have a valid file extension (fa,fna,fasta,fq,fastq,bam,csfasta,csfastq)",call.=FALSE)
  } 
}

.consolidateFileExtensions <- function(filename, compressed=FALSE) {
    # convert various sequence file extension to one representative
    if(compressed)
        fileExtension <- tolower(tools::file_ext(sub("[.](gz|bz2|xz)$","",filename)))
    else
        fileExtension <- tolower(tools::file_ext(filename))

    fileExtension[ fileExtension %in% c("fa", "fna", "fasta") ] <- "fasta"
    fileExtension[ fileExtension %in% c("fq", "fastq") ] <- "fastq"
    return(fileExtension)
}

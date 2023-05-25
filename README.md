Table of Contents
=================
<!--ts-->
* [Table of Contents](#table-of-contents)
* [EasyCircR](#easycircr)
   * [Overview](#overview)
   * [Installation](#installation)
      * [Installation from R](#installation-from-r)
      * [Installation with Docker](#installation-with-docker)
   * [Testing dataset](#testing-dataset)
   * [Pipeline](#pipeline)
      * [1. EISA](#1-eisa)
      * [2. CIRI-Full](#2-ciri-full)
      * [3. circRNA DE](#3-circrna-de)
      * [4. miRNA binding](#4-mirna-binding)
      * [5. Interaction between circRNAs and genes](#5-interaction-between-circrnas-and-genes)
      * [6. Further analysis and Shiny](#6-further-analysis-and-shiny)
   *  [Example of analysis with docker](#example-of-analysis-with-docker)
<!-- Added by: runner, at: Fri Jun 11 08:39:18 UTC 2021 -->

<!--te-->

# EasyCircR

## Overview 

EasyCircR combines the detection and reconstruction of circular RNAs (circRNAs) with the exon-intron split analysis (EISA) by the means of miRNA response element (MREs) prediction.

Starting from RNA-seq data, it uses CIRI-full for the reconstruction of
full-length circRNAs and EISA for the detection of genes which are
post-transcriptional regulated. In order to combine the two results, miRNA
response elements (MREs), with respective binding miRNA, are predicted using
TargetScan for each circRNA.   
The relationship miRNA-gene is retrieved using MultimiR package and results are
then filtered by genes which where previously found post-transcriptional
regulated.

The pipeline is divided in 6 major steps.

1. Measure changes across different experimental conditions to 
quantify post-transcriptional regulation of gene expression (EISA).
2. Reconstruct and quantify full-length circRNAs with CIRI-Full.
3. Differentially expressed analysis of detected circRNAs (with limma).
4. Predict MREs from the full sequence of circRNA using TargetScan.
5. Finally, connecting circRNA to genes.
6. Further analysis and visualization with Shiny app.

## Installation

We prodide two ways to download and install EasyCircR package:
 
- Installation from R
- Installation with Docker  


##  Installation from R

the package relies on several dependencies of libraries and tools that must be installed before proceeding with package installation including:

- R ( >=4.0 ) 
- java (suggested: openjdk-11-jre-headless) 
- bwa (https://bio-bwa.sourceforge.net/)
- bedtools (https://bedtools.readthedocs.io/en/latest/)

Ubuntu users need to install also: 
- curl
- libssl-dev
- libfontconfig1-dev
- libcurl4-openssl-dev
- libxml2-dev
- libpoppler-cpp-dev
- libharfbuzz-dev
- libfribidi-dev
- libfreetype6-dev
- libpng-dev
- libtiff-dev
- libjpeg-dev
- libgit2-dev

Once that all the dependencies are installed EasyCircR can be installed from github using  "devtools" package.

```R
# install devtools package
install.packages('devtools')

# install EasyCircR
devtools::install_github('InfOmics/EasyCircR')
```

## Installation with Docker

The package is also available inside a dedicated docker container. Follow the following links to install Docker on [Linux](https://docs.docker.com/desktop/install/linux-install/), [Windows](https://docs.docker.com/docker-for-windows/install/) or [macOS](https://docs.docker.com/desktop/install/mac-install/), and follow the on-screen instructions.


To pull the EasyCircR docker image run the following command:

```bash
docker pull savesani/easycircr
```

An example on how to run EasyCircR analysis with docker is provided at the [Example of analysis with docker](#example-of-analysis-with-docker) section. 

## Testing dataset
Tests were performed using the data described in: https://doi.org/10.3390/ncrna7020026.
Analysis was performed on 6 RNA-seq samples (pair-end) of TMD8 cell line collected form  diffuse large B cell lymphoma (DLBCL) cells treated with DMSO (`control`) and PQR-309 (`treatment`). 
The script used to download and index reference genome and annotation files is available on the supplementary branch. 

NB: Before running CIRI-Full make sure to index the genome with `bwa`.


## Pipeline

In order to use as input the fastq files, their path must be stored in a tab separated file with three named columns like the one present in the supplementary branch.  
For each sample you need to indicate the two pair-end fastq files and a sample name.

EasyCircR requires a significant amount of disk memory to store temporary files that will be removed at the end of each pipeline step. It is recommended to start the analysis having at least 200 GB of free disk space.

```
FileName1	FileName2	SampleName
fastq/S1_DMSO_4h_R1.fq.gz	fastq/S1_DMSO_4h_R2.fq.gz	S1_DMSO_4h
fastq/S2_PQR_4h_R1.fq.gz	fastq/S2_PQR_4h_R2.fq.gz	S2_PQR_4h
fastq/S3_DMSO_8h_R1.fq.gz	fastq/S3_DMSO_8h_R2.fq.gz	S3_DMSO_8h
fastq/S4_PQR_8h_R1.fq.gz	fastq/S4_PQR_8h_R2.fq.gz	S4_PQR_8h
fastq/S5_DMSO_12h_R1.fq.gz	fastq/S5_DMSO_12h_R2.fq.gz	S5_DMSO_12h
fastq/S6_PQR_12h_R1.fq.gz	fastq/S6_PQR_12h_R2.fq.gz	S6_PQR_12h
```

Since many of the steps in the pipeline are very long (even 5/6h per sample) we
store the final results in a folder `EasyCircR`, in the same directory of the code.  
The final structure of `EasyCircR` is something similar to this:

```bash
EasyCirc
├── circRNA
│   ├── CIRI-Full
│   │   ├── countMatrix.rds
│   │   ├── predictedCircRNAs.rds
│   │   ├── S5_DMSO_12h
│   │   │   └── CIRI-vis_out
│   │   ├── S1_DMSO_4h
│   │   │   └── CIRI-vis_out
│   │   ├── S3_DMSO_8h
│   │   │   └── CIRI-vis_out
│   │   ├── S6_PQR_12h
│   │   │   └── CIRI-vis_out
│   │   ├── S2_PQR_4h
│   │   │   └── CIRI-vis_out
│   │   └── S4_PQR_8h
│   │       └── CIRI-vis_out
│   ├── miRNA
│   │   ├── circ_mirna.txt
│   │   └── circRNA_sequences.txt
│   └── trimmed_fastq
├── geneMirnaCirc
│   └── geneMirnaCirc.rds
└── postGene
    └── postGene.rds
```

The three major folders (`circRNA`, `geneMirnaCirc`, `postGene`) contain each
one or two `.rds` file storing the results of the pipeline, e.g.
`countMatrix.rds` contains the counts of circRNA found or `postGene.rds` the
results of `get_postregulated_gene`.
The `CIRI-Full` folder contains a folder for each sample in which are stored 
all the circRNAs in a pdf form.  
Most of the functions in the package also work as a loading function, checking if
there is the respective `.rds` file and returning it, instead of recomputing it.  
With these functions we can be bypass this behaviour with the parameter `force=TRUE`.

### 1. EISA

Eisa requires a `samples_file`, the genome with its annotations and the vector `condition`.

```R
samples_file <- "samples.txt"
genome_file <- "genome/hg38.fa"
genome_annotation_file <- "genome/hg38.gtf"
condition <- factor(rep(c("DMSO", "PQR"),3))
```

Important to note that in `condition`, the i-th position represents the i-th
sample in the tab separated file.

The constrat matrix for the DE in eisa is build automatically by 
`eisaR` based on the `condition` vector defined as secondLevel - firstLevel.  
In this case we do not need to change the level order of the `PQR` and `DMSO`
since we want to `PQR` - `DMSO`.

```R
levels(condition)
#[1] "DMSO" "PQR"
```

The following runs eisa on the fastq files. This will also return non statistically significant genes. So we should filter them based on `FDR`.

```R
post_gene <- get_postregulated_gene(samples_file, genome_file, 
                                    genome_annotation_file,
                                    condition, aligner="Rhisat2", 
                                    force=FALSE, n_core=1)

post_gene <- post_gene[post_gene$FDR <= 0.1,]
```

```R
head(post_gene)

                     logFC    logCPM        F       PValue          FDR
ENSG00000083845  1.123991 10.977869 119.44677 5.834036e-09 6.069731e-05
ENSG00000188846  1.335570  9.392296  85.82574 6.176369e-08 2.554220e-04
ENSG00000127824 -2.588764  6.016761  83.68345 7.365108e-08 2.554220e-04
ENSG00000156508  1.434746 12.825495  73.82974 1.744079e-07 3.506286e-04
ENSG00000137154  1.149152 10.647291  72.52710 1.968497e-07 3.506286e-04
ENSG00000117395 -1.183689  6.821934  72.24068 2.022080e-07 3.506286e-04

```



### 2. CIRI-Full

Since CIRI-Full relies on CIRI-AS, which was developed before CIRI2, 
it cannot process sequencing reads with different lengths.   
The function `run_ciri_full` accept the parameter `trim_reads_length` which
can take a number that define at which length reads must be cutted. Reads smaller
than that size are discarded and longer are trimmed.  
In order to keep most of the reads and as long as possible, the function
`plot_readlen_hist` can help.  
In the dataset all the reads have already been trimmed and the function returns
the length of the reads.

```R
plot_readlen_hist(samples_file)
#All the reads have the same length: 122 bp
```

Otherwise it would have returned an histogram like the following:

![Image](img/2021-05-26-03:52:10.png)

Where the blue dotted line represents at which length cutting the reads will
result in keeping 80% of all the original reads.

We can now run CIRI-Full with the following functions.

```R
run_ciri_full(samples_file, genome_file, genome_annotation_file, 
              trim_reads_length=130, force=FALSE)
```

After CIRI-Full has finished we can read the output using:

```R
circ <- read_ciri_output(samples_file, chrMT = T)
names(circ)
#[1] "circ_df"  "circ_mtx"
circ_df <- circ$circ_df 
circ_mtx <- circ$circ_mtx
```

`chrMT` parameter can be set to FALSE if the user wants to exclude circRNAs detected on MT chromosome.
It returns two objects: `circ_df` and `circ_mtx`.  
`circ_df` contains information like the chromosome, the start and end of the
back-spliced junction (BSJ), the length of the reconstructed circRNA and the
reconstructed sequence.  
`circ_mtx` is the count matrix of circRNAs.

### 3. circRNA DE

EasyCirc provide the function `de_circrna` to perfom the differential expression (DE) analysis  on circRNAs with `limma`.  
The user can of course use any other program he might prefer.  

First of all we create the design and constrast matrix.

```R
design <- model.matrix(~0+condition)
contr <- limma::makeContrasts("PQRvsDMSO"  = conditionPQR - conditionDMSO, levels = colnames(design))
```

[A guide to creating design matrices for gene expression experiments](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)

Than we run the DE analysis with the following:

```R
circ_de <- de_circrna(circ_mtx, condition, design, contr, lfc=0, p.value=0.05,
                      voomWithQualityWeights=FALSE, min_num_samples=2, plotVoomResult=FALSE, ...)
```

This function also perform filtering based on the parameter `min_num_samples`,
removing all the circRNAs that are not present in at least a specified number of samples.

We filter `circ_df` only for DE circRNAs.

```R
circ_df_de <- circ_df[circ_df$bsj_id %in% rownames(circ_de), ]
```

### 4. miRNA binding

`get_mirna_binding` accept a `data.frame` of circRNAs and for each of them
 predicts the possible MREs by the means of TargetScan.

```R

circ_mirna <- get_mirna_binding(circ_df_de, force=FALSE)

head(circ_mirna[,c(1,2)])

             a_Gene_ID          miRNA_family_ID
 12:116230533|116237705:-:13      miR-101-5p
 12:116230533|116237705:-:13      miR-105-5p
 12:116230533|116237705:-:13      miR-10b-3p
 12:116230533|116237705:-:13      miR-1180-5p
 12:116230533|116237705:-:13      miR-1184
 12:116230533|116237705:-:13      miR-1205
```

The result is a `data.frame` that associates each circRNA to one or more
miRNA families.

### 5. Interaction between circRNAs and genes

The final step requires to join the predicted relations circRNA-miRNA with
known/predicted miRNA-gene interactions. The miRNA-gene relation is downloaded using `mulitMiR` from predicted and validated miRNA–target interactions databases and/or from not validated ones.

```R
gene_mirna_circ <- connect_circ_gene(circ_mirna, post_gene, only_significant_genes=TRUE, force=FALSE, tabletype="validated")

```
`tabletype` parameter can be set to `validated` to retrieve only validated miRNA-gene interactions or `all` to retrieve all the potential interactions.
The parameter `only_significant_genes` does not make any difference if in the
first step of eisa we filtered `post_gene` by `FDR`.  

### 6. Further analysis and Shiny  

EasyCirc provides a Shiny app to visualize the connection between
circRNA-miRNA-gene. The app includes multiple types of filter, the ability to
save filtered results in csv/excel and the possibility to further investigate
reconstructed circRNAs.

```R
launch_shiny(shiny_host ="0.0.0.0", shiny_port = 3838)
```

We can also show the circRNA structure directly from R without launching the
Shiny app using `plot_circ`

```R
plot_circ(circ_df,"12:116230533|116237705:-:13")
```

That shows the results of CIRI-vis, which is part of the CIRI-Full algorithm.

Finally we can plot the regulatory nerwork.

```R
plot_regulatory_net("12:116230533|116237705:-:13")
```


## Example of analysis with docker 

We provide in the supplementary branch an R script called EasyCircR_test.R that could be downloaded, modidified and executed inside the docker container. 
Remember that to run the pipeline, is necessary to have the reference genome indexed.

```bash 
wget https://github.com/InfOmics/EasyCircR/tree/supplementary/Supplementary_data/EasyCircR_test.R
```

Docker container can be used to run directly the script:

```bash 
# run container
docker run --rm -p 9000:3838 -v /volume_dir:/volume_dir/ savesani/easycircr Rscript EasyCirc_test.R 

```
Or we can run the docker container in the interactive mode and execute each step of the analysis individually 


```bash 
# run container in the interactive mode
docker run -ti --rm -p 9000:3838 -v /volume_dir:/volume_dir/ savesani/easycircr  
```

After that the analysis is completed, Shiny app can be run keeping the container active and open the following link in a browser page: http://0.0.0.0:9000





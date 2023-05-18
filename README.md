Table of Contents
=================
<!--ts-->
* [Table of Contents](#table-of-contents)
* [EasyCirc](#easycirc)
   * [Overview](#overview)
   * [Installation](#installation)
      * [1. Installation via Docker](#1-installation-via-docker)
      * [2. Installation via R](#2-installation-via-R)
   * [Testing dataset](#testing-dataset)
   * [Pipeline](#pipeline)
      * [1. EISA](#1-eisa)
      * [2. CIRI-Full](#2-ciri-full)
      * [3. circRNA DE](#3-circrna-de)
      * [4. miRNA binding](#4-mirna-binding)
      * [5. Interaction between circRNAs and genes](#5-interaction-between-circrnas-and-genes)
      * [6. Further analysis and Shiny](#6-further-analysis-and-shiny)

<!-- Added by: runner, at: Fri Jun 11 08:39:18 UTC 2021 -->

<!--te-->

# EasyCirc

## Overview 

EasyCirc combines the detection and reconstruction of circular RNAs (circRNAs) with the
exon-intron split analysis (
) by the means of miRNA response element (MREs)
prediction.

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
4. Predict MREs from the full sequence of circRNA (TargetScan, miRanda,
   RNA-hybrid).
5. Finally, connecting circRNA to genes.
6. Further analysis and visualization with shiny app.

## Installation

```R
install.packages('devtools')
devtools::install_github('InfOmics/EasyCircR')
```

It also require bwa already installed to run CIRI-Full.

### 1. Installation via Docker
asd

### 2. Installation via R
asd

## Testing dataset

Fastq files to test the package can be download [here](blank).  
The zip files contains 6 RNA-seq samples (pair-end) of ABC DLBCL cell line TMD8
treated with DMSO (`control`) or PQR-309 (`treatment`).  
Other than fastq files, a filtered ensmbl HG38 genome (only `chr5`, `chr7`) and
filtered genome annotation file are included in the zip.

Before running CIRI-Full make sure to index the genome with `bwa`.

```bash
bwa index hg38_chr3_chr7.fa
```

## Pipeline

In order to pass the fastq files, their path must be stored in a tab separated 
file with three named columns like the one present in the testing dataset.  
For each sample you need to indicate the two pair-end fastq files and a sample name.

```
FileName1	FileName2	SampleName
fastq/TMD8_DMSO_4h_R1.fq.gz	fastq/TMD8_DMSO_4h_R2.fq.gz	TMD8_DMSO_4h
fastq/TMD8_PQR_4h_R1.fq.gz	fastq/TMD8_PQR_4h_R2.fq.gz	TMD8_PQR_4h
fastq/TMD8_DMSO_8h_R1.fq.gz	fastq/TMD8_DMSO_8h_R2.fq.gz	TMD8_DMSO_8h
fastq/TMD8_PQR_8h_R1.fq.gz	fastq/TMD8_PQR_8h_R2.fq.gz	TMD8_PQR_8h
fastq/TMD8_DMSO_12h_R1.fq.gz	fastq/TMD8_DMSO_12h_R2.fq.gz	TMD8_DMSO_12h
fastq/TMD8_PQR_12h_R1.fq.gz	fastq/TMD8_PQR_12h_R2.fq.gz	TMD8_PQR_12h
```

Since many of the steps in the pipeline are very long (even 5/6h per sample) we
store the final results in a folder `EasyCirc`, in the same directory of the code.  
The final structure of `EasyCirc` is something similar to this:

```bash
EasyCirc
├── circRNA
│   ├── CIRI-Full
│   │   ├── countMatrix.rds
│   │   ├── predictedCircRNAs.rds
│   │   ├── TMD8_DMSO_12h
│   │   │   └── CIRI-vis_out
│   │   ├── TMD8_DMSO_4h
│   │   │   └── CIRI-vis_out
│   │   ├── TMD8_DMSO_8h
│   │   │   └── CIRI-vis_out
│   │   ├── TMD8_PQR_12h
│   │   │   └── CIRI-vis_out
│   │   ├── TMD8_PQR_4h
│   │   │   └── CIRI-vis_out
│   │   └── TMD8_PQR_8h
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
samples_file <- "test/samples_test.txt"
genome_file <- "genome/hg38_chr3_chr7.fa"
genome_annotation_file <- "genome/hg38_chr3_chr7.gtf"
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

The following runs eisa on the fastq files.
```R
post_gene <- get_postregulated_gene(samples_file, genome_file, genome_annotation_file, 
                                    condition, aligner="Rhisat2", 
                                    force=FALSE, n_core=1)
```

```R
head(post_gene)
                     logFC    logCPM        F       PValue          FDR
ENSG00000188846  1.2034417 11.734689 73.89589 3.098571e-08 3.383639e-05
ENSG00000085662  1.1478014  9.166539 36.27856 6.207244e-06 3.389155e-03
ENSG00000163931  0.7139964 11.284703 31.49309 1.566549e-05 5.702238e-03
ENSG00000128536  0.6426320 11.436869 29.49626 2.367210e-05 6.462483e-03
ENSG00000121879 -0.6312345 10.713384 26.97729 4.087103e-05 8.926233e-03
ENSG00000163584  0.7476950  9.533056 20.71156 1.837729e-04 3.344667e-02
nrow(post_gene)
#[1] 1092
```

This will also return non statistically significant genes.  
So we should filter based on `PValue` and `logFC`.

```R
post_gene <- post_gene[post_gene$PValue <= 0.05 & abs(post_gene$logFC) >= 1,]
nrow(post_gene)
#[1] 15
```

### 2. CIRI-Full

Since CIRI-Full relies on CIRI-AS, which was developed before CIRI2, 
it cannot process sequencing reads with different lengths.   
The function `run_ciri_full` accept the parameter `trimReadsLength` which
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
              trimReadsLength=NULL, force=FALSE)
```

Keep in mind that this is the slowest step in the pipeline and it takes on
average ~5/6h for a sample composed by two pair-end fastqs of ~2.5GB.
For the test dataset, since they are only 130MB each, it will only take ~10m per
sample.  
(Based on my CPU Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz and 16GB of RAM).  
The slowest step of the CIRI-Full algorithm is when is 
`Waiting for the result of CIRI-AS` but since it is not meant to print anything in the
process it can give the false idea that is crushed, while it is not.


After CIRI-Full has finished we can read the output using
```R
circ <- read_ciri_output(samples_file, chrMT = T)
names(circ)
#[1] "circ_df"  "circ_mtx"
circ_df <- circ$circ_df 
circ_mtx <- circ$circ_mtx
```

It returns two `data.frame`: `circ_df` and `circ_mtx`.  
`circ_df` contains information like the chromosome, the start and end of the
back-spliced junction (BSJ), the length of the reconstructed circRNA and the
reconstructed sequence.  
`circ_mtx` is the count matrix of circRNAs.

### 3. circRNA DE

EasyCirc provide the function `de_circrna` to perom the DE on circRNAs with `limma`.  
The user can of course use any other program he might prefer.  

First of all we create the design and constrast matrix.
```R
design <- model.matrix(~0+condition)
contr <- limma::makeContrasts("PQRvsDMSO"  = conditionPQR - conditionDMSO, levels = colnames(design))
```

[A guide to creating design matrices for gene expression experiments](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)

Than we run the DE analysis with the following:
```R
circ_de <- de_circrna(circ_mtx, condition, design, contr, lfc=1, p.value=0.05,
                      voomWithQualityWeights=FALSE, min_num_samples=2)
```

This function also perform filtering based on the parameter `min_num_samples`,
removing all the circRNAs that are not present in at least a specified number of
samples.

We filter `circ_df` only for DE circRNAs.

```R
circ_df_de <- circ_df[circ_df$bsj_id %in% rownames(circ_de), ]
```

### 4. miRNA binding

`get_mirna_binding` accept a `data.frame` of circRNAs and for each of them
it predicts the possible MREs by the means of TargetScan.

```R
circ_mirna <- get_mirna_binding(circ_df_de, force=FALSE)
head(circ_mirna[,c(1,2)])
#.             a_Gene_ID miRNA_family_ID
#7:115974170|115984513:-     miR-1185-5p
#7:115974170|115984513:-     miR-1224-5p
#7:115974170|115984513:-     miR-1252-5p
#7:115974170|115984513:-        miR-1276
#7:115974170|115984513:-     miR-1304-5p
#7:115974170|115984513:-      miR-15b-3p
```

The result is a `data.frame` that associates each circRNA to one or more
miRNA families.

### 5. Interaction between circRNAs and genes

The final step requires to join the predicted relations circRNA-miRNA with
known/predicted miRNA-gene. The miRNA-gene is downloaded using `mulitMiR` from
predicted miRNA–target interactions databases like miRanda, miRDBmiR and
TargetScan, and/or from validated miRNA–target interactions databases like
TarBase or miRTarBase.

```R
gene_mirna_circ <- connect_circ_gene(circ_mirna, post_gene, only_significant_genes=TRUE, force=FALSE, tabletype="validated")
```

The parameter `only_significant_genes` does not make any difference if in the
first step of eisa we filtered `post_gene` by `PValue` and `logFC`.  

### 6. Further analysis and Shiny

EasyCirc provides a shiny app to visualize the connection between
circRNA-miRNA-gene. The app includes multiple types of filter, the ability to
save filtered results in csv/excel and the possibility to further investigate
reconstructed circRNAs.

```R
launch_shiny(shiny_host ="127.0.0.1", shiny_port = 7775)
```

We can also show the circRNA structure directly from R without launching the
shiny app using `plot_circ`

```R
plot_circ(circ_df,"7:116092137|116112038:-:9")
```

That shows the results of CIRI-vis, which is part of the CIRI-Full algorithm.

Finally we can plot the regulatory nerwork.

```R
plot_regulatory_net("7:116092137|116112038:-:9")
```

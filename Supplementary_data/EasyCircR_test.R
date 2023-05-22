
library(EasyCircR)

# set your working directory 
setwd("/my_dir/")

samples_file <- "samples.txt"
genome_file <- "genome/hg38.fa"
genome_annotation_file <- "genome/hg38.gtf"
condition <- factor(rep(c("DMSO", "PQR"),3))



#---------------------EISA
post_gene <- get_postregulated_gene(samples_file, genome_file, genome_annotation_file, 
                                     condition, aligner="Rhisat2",
                                     force=FALSE, n_core=15)

post_gene <- post_gene[post_gene$FDR <= 0.1,]


#---------------------CIRI-full
run_ciri_full(samples_file, genome_file, genome_annotation_file, force=FALSE, trim_reads_length=130, n_core= 15, remove_temporary_files=TRUE)

circ <- read_ciri_output(samples_file)
circ_df <- circ$circ_df
circ_mtx <- circ$circ_mtx


#---------------------circDE
design <- model.matrix(~0+condition)
contr <- limma::makeContrasts("PQRvsDMSO"  = conditionPQR - conditionDMSO, levels = colnames(design))

circ_de <- de_circrna(circ_mtx, condition, design, contr, lfc=0, p.value=0.05, 
                      voomWithQualityWeights=FALSE, min_num_samples=2)

circ_df_de <- circ_df[circ_df$bsj_id %in% rownames(circ_de), ]


#---------------------miRNA binding
circ_mirna <- get_mirna_binding(circ_df_de, force=F)
head(circ_mirna[,c(1,2)])


#---------------------Interaction between circRNAs and genes
gene_mirna_circ <- connect_circ_gene(circ_mirna, post_gene, only_significant_genes=TRUE, force=F, tabletype="all")


#---------------------Run Shiny app

launch_shiny(shiny_host ="0.0.0.0", shiny_port = 3838)

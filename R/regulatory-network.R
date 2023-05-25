#' @title EasyCircR - plot the regulatory network
#'
#' @description Plot the circRNA-miRNA-gene regulatory network.
#' 
#' @author Luca Parmigiani, Antonino Aparo, Simone Avesani
#' 
#' @param bsj_id the Back-Splice junction ID of the circRNA user wants to display.
#'
#' @param gene_mirna_circ the \code{data.frame} as results of \code{EasyCircR::connect_circ_gene(...)}.
#' By default (\code{NULL}) the function check if the \code(geneMirnaCirc.rds) is stored in "EasyCirc/geneMirnaCirc",
#' otherwise execute \code{EasyCircR::connect_circ_gene()}.
#' 
#' @return 
#' 
#' @examples 
#' plot_regulatory_net("12:116230533|116237705:-:13")
#'
#' @importFrom igraph graph_from_data_frame as_adjacency_matrix
#' @importFrom network network network.vertex.names
#' @import ggnet
#' @export
plot_regulatory_net = function (bsj_id, gene_mirna_circ=NULL) {
  library(network)
  library(ggnet)
  library(devtools)
  library(ggplot2)
  if (is.null(gene_mirna_circ)) {
    gene_mirna_circ <- connect_circ_gene()
  }
  gene_mirna_circ <-  gene_mirna_circ[gene_mirna_circ$circRNA_id == bsj_id ,]
  genes <- gene_mirna_circ$target_symbol %>% unique()
 # genes <- genes[genes != ""]
  mirnas <- gene_mirna_circ$mature_mirna_id %>% unique()
  circrnas <- gene_mirna_circ$circRNA_id %>% unique()
    
  mirnaGene <- gene_mirna_circ[,c("mature_mirna_id","target_symbol")] %>% unique()
  circMirna <- gene_mirna_circ[,c("circRNA_id", "mature_mirna_id")]
  colnames(mirnaGene) = c("u","v")
  colnames(circMirna) = c("u","v")

  vertices <- data.frame(id = c(genes, mirnas, circrnas), 
                         type=c(rep("gene",length(genes)), 
                                rep("miRNA",length(mirnas)),
                                rep("circRNA",length(circrnas))))
  
  edges <- rbind(mirnaGene, circMirna)
  G <- igraph::graph_from_data_frame(d=edges,vertices=vertices,directed=F)
  m <- igraph::as_adjacency_matrix(G)
  net = network::network(as.matrix(m), directed = FALSE)

  net %v% "type" = ifelse( network::network.vertex.names(net) %in% mirnas, "miRNA", "gene")
  set.vertex.attribute(net,"type", value = "circRNA", v = which(network.vertex.names(net) %in% circrnas))
  network.vertex.names(net)[(length(network::network.vertex.names(net)))] = ""
  col = c("circ" = "black", "mirna" = "gold", "gene" = "lightblue")  

  ggnet::ggnet2(net, node.color = "type",  palette = "Set2", alpha = 0.75,
         edge.alpha = 1/3,
         max_size = 10, node.label = T,
         legend.size = 12) +
         coord_equal() + 
         ggtitle(bsj_id) +
         theme(plot.title = element_text(color="black", size=14, face="bold.italic"))
}
   


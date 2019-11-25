


get_net_ggplot <- function (plotcord, color, name){
  pl <- ggplot(plotcord)  +
    geom_segment(
      data = edges,
      aes_(
        x =  ~ X1,
        y =  ~ Y1,
        xend =  ~ X2,
        yend =  ~ Y2
      ),
      size = 0.5,
      alpha = 0.9,
      colour = "gray25"
    ) +
    geom_point(aes_(
      x =  ~ X1,
      y =  ~ X2,
      size =  ~ WeightedDegree,
      alpha =  ~ WeightedDegree
    ),
    color = color) +
    geom_label_repel(
      aes_(
        x =  ~ X1,
        y =  ~ X2,
        label =  ~ vertex.names,
        color =  ~ Type
      ),
      box.padding = unit(1, "lines"),
      data = function(x) {
        x[x$shouldLabel,]
      }
    ) +
    scale_colour_manual(values = c("Gene" = "#005E87",
                                   "Topic" = "#540814")) +
    labs(title = name) +
    ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white",
                                               colour = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  return(pl)
}
get_sum_of_weights <- function(plot_counts){
  sum_of_weights1 <-
    plot_counts %>% group_by(Genes) %>%   summarise(Weight = sum(count))
  sum_of_weights2 <-
    plot_counts %>% group_by(Topic) %>%   summarise(Weight = sum(count))
  force_bind = function(df1, df2) {
    colnames(df2) = colnames(df1)
    rbind(df1, df2)
  }
  
  sum_of_weights <- force_bind(sum_of_weights1, sum_of_weights2)
  WeightedDegree <- sum_of_weights$Weight
  names(WeightedDegree) <- sum_of_weights$Genes
  return(WeightedDegree)
}
get_fruchterman_reingold_coordinates <- function(network_object_adjacency_matrix){
  return(data.frame(sna::gplot.layout.fruchtermanreingold(network_object_adjacency_matrix, NULL)))
}
get_edges_from_network_edgelist <- function(network_object_edgelist){
  data.frame(plotcord[network_object_edgelist[, 1], ], plotcord[network_object_edgelist[, 2], ])
}
extract_attribute_to_coordinates <- function(plotcord, network_object, attribute){
  plotcord[attribute] <-
    as.factor(network::get.vertex.attribute(network_object, attribute))
  return(plotcord)
}
set_coordinate_types <- function(plot_counts, plotcord){
  plotcord[, "Type"] <- ""
  plotcord[which(plotcord$vertex.names %in% plot_counts$Genes), "Type"] <-
    "Gene"
  plotcord[which(plotcord$vertex.names %in% plot_counts$Topic), "Type"] <-
    "Topic"
  
  return(plotcord)
}
set_genes_to_label <- function(plotcord, max_n){
  plotcord[, "shouldLabel"] <- FALSE
  int_hubs <- names(sort(WeightedDegree, decreasing = TRUE))[seq_len(max_n)]
  int_bool <- plotcord[, "vertex.names"] %in% int_hubs
  
  mod_genes <- names(WeightedDegree[names(WeightedDegree) %in% plot_counts$Genes])
  sel_genes <-
    names(sort(WeightedDegree[names(WeightedDegree) %in% plot_counts$Genes], decreasing = TRUE))[seq_len(max_n)]
  sel_vertex <-
    c(sel_genes, names(WeightedDegree[names(WeightedDegree) %in% plot_counts$Topic]))
  plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <-
    TRUE
  
  return(plotcord)
}

#' #'plot_literature_graph
#'
#'Plot a graph inspired in CEMiTool's graphs
#'
#'
#' @param plot_counts The dataframe returned from the get_literature_score function
#' @param name The name of the plot.
#' @param color The color of the plot. Defaults to a shade of red ("#B30000FF").
#' @param max_number_of_labels The max number of gene labels to show.  Defaults to 10.
#' @import igraph
#' @import ggrepel
#' @import ggplot2
#' @return A ploty/ggplot2 object is either returned or directly plotted
#' @export
#' @examples
#'   gene <- c('CD4','CD14', "AIF1", "ACVR1", "CDY2A")
#'   terms_of_interest <- c("CD4 T cell", "CD14+ Monocyte", "B cell", "CD8 T cell",
#'                             "FCGR3A+ Monocyte", "NK cell", "Dendritic cell", "Megakaryocyte", 'immunity')
#'   literature_counts <- get_literature_score(gene, terms_of_interest)
#'   pl <- plot_literature_graph(literature_counts, name = 'test')
#'   pl
 

plot_literature_graph <-
  function(plot_counts,
           name,
           color = "#B30000FF",
           max_number_of_labels = 10) {
    
    
    plot_counts <- plot_counts[plot_counts$count > 0,]
    
    WeightedDegree <- get_sum_of_weights(plot_counts)
    max_n <- min(max_number_of_labels, sum(names(WeightedDegree) %in% plot_counts$Genes))
    
    degrees <- igraph::degree(igraph_object, normalized = FALSE)
    
    igraph_object <- igraph::graph.data.frame(plot_counts)
    igraph_object <-
    igraph::set_vertex_attr(igraph_object, "WeightedDegree", value = WeightedDegree)
    
    network_object <- intergraph::asNetwork(igraph_object)
    network_object_adjacency_matrix <-  network::as.matrix.network.adjacency(network_object) 
    network_object_edgelist <- network::as.matrix.network.edgelist(network_object)
    
    plotcord <- get_fruchterman_reingold_coordinates(network_object_adjacency_matrix)
      
    colnames(plotcord) <- c("X1", "X2")

    edges <- get_edges_from_network_edgelist(network_object_edgelist)
    
    plotcord <- extract_attribute_to_coordinates(plotcord, network_object, attribute = "vertex.names")
    plotcord <- extract_attribute_to_coordinates(plotcord, network_object, attribute = "WeightedDegree")
    plotcord$WeightedDegree <- as.numeric(plotcord$WeightedDegree)
    plotcord <- set_coordinate_types(plot_counts, plotcord)
    plotcord <- set_genes_to_label(plotcord, max_n)
    
    plotcord$Degree_cut <-
      cut(plotcord$WeightedDegree,
          breaks = 3,
          labels = FALSE)
    colnames(edges) <-  c("X1", "Y1", "X2", "Y2")
    plotcord$in_mod <- TRUE
    not_in <- setdiff(plotcord[, "vertex.names"], mod_genes)
    plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <-
      FALSE
    
    pl <- get_net_ggplot(plotcord, color, name)
    return(pl)
  }

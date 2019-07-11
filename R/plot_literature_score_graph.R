#' #'plot_literature_graph
#'
#'Plot a graph inspired in CEMiTool's graphs
#'
#'
#' @param plot_counts The dataframe returned from the get_literature_score function
#' @param name The name of the plot.
#' @param color The color of the plot. Defaults to a shade of red ("#B30000FF").
#' @param n The max number of gene labels to show.  Defaults to 10.
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
           n = 10) {
    plot_counts <- plot_counts[plot_counts$count > 0,]
    ig_obj <- igraph::graph.data.frame(plot_counts)
    sum_of_weights1 <-
      plot_counts %>% group_by(Genes) %>%   summarise(Weight = sum(count))
    sum_of_weights2 <-
      plot_counts %>% group_by(Topic) %>%   summarise(Weight = sum(count))
    
    
    force_bind = function(df1, df2) {
      colnames(df2) = colnames(df1)
      rbind(df1, df2)
    }
    
    sum_of_weights <- force_bind(sum_of_weights1, sum_of_weights2)
    
    force_bind = function(df1, df2) {
      colnames(df2) = colnames(df1)
      rbind(df1, df2)
    }
    
    sum_of_weights <- force_bind(sum_of_weights1, sum_of_weights2)
    
    
    
    
    
    
    
    
    degrees <- igraph::degree(ig_obj, normalized = FALSE)
    sumwts <- sum_of_weights$Weight
    names(sumwts) <- sum_of_weights$Genes
    ig_obj <-
      igraph::set_vertex_attr(ig_obj, "sumwts", value = sumwts)
    net_obj <- intergraph::asNetwork(ig_obj)
    m <-
      network::as.matrix.network.adjacency(net_obj) # get sociomatrix
    # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
    plotcord <-
      data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
    # or get it them from Kamada-Kawai's algorithm:
    # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
    colnames(plotcord) <- c("X1", "X2")
    edglist <- network::as.matrix.network.edgelist(net_obj)
    edges <-
      data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
    plotcord$vertex.names <-
      as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
    plotcord$WeightedDegree <-
      network::get.vertex.attribute(net_obj, "sumwts")
    plotcord[, "shouldLabel"] <- FALSE
    plotcord[, "Type"] <- ""
    
    max_n <- min(n, sum(names(sumwts) %in% plot_counts$Genes))
    int_hubs <- names(sort(sumwts, decreasing = TRUE))[seq_len(max_n)]
    int_bool <- plotcord[, "vertex.names"] %in% int_hubs
    plotcord[which(plotcord$vertex.names %in% plot_counts$Genes), "Type"] <-
      "Gene"
    plotcord[which(plotcord$vertex.names %in% plot_counts$Topic), "Type"] <-
      "Topic"
    mod_genes <- names(sumwts[names(sumwts) %in% plot_counts$Genes])
    sel_genes <-
      names(sort(sumwts[names(sumwts) %in% plot_counts$Genes], decreasing = TRUE))[seq_len(max_n)]
    sel_vertex <-
      c(sel_genes, names(sumwts[names(sumwts) %in% plot_counts$Topic]))
    colnames(edges) <-  c("X1", "Y1", "X2", "Y2")
    #edges$midX  <- (edges$X1 + edges$X2) / 2
    #edges$midY  <- (edges$Y1 + edges$Y2) / 2
    plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <-
      TRUE
    plotcord$Degree_cut <-
      cut(plotcord$WeightedDegree,
          breaks = 3,
          labels = FALSE)
    plotcord$in_mod <- TRUE
    #mod_genes <- cem@module[cem@module$modules==name,]$genes
    not_in <- setdiff(plotcord[, "vertex.names"], mod_genes)
    plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <-
      FALSE
    
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

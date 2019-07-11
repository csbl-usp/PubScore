# lubianat 28/09/2018

#' plot_literature_score
#'
#' Plots a non-clusterized heatmap of the article counts for the combination of
#' gene list and list of terms
#' NOTE: the object has to be exactly the one returned by get_literature_score.R .
#' Otherwise ggplot2 will not be able to identify the correct columns.
#'
#' @param plot_counts The dataframe returned from the get_literature_score function
#' @param return_ggplot If TRUE, returns a ggplot2 object instead of plotting. Defaults to FALSE.
#' @param is_plotly If TRUE, a interactive plot is plotted in the place o static ggplot. Defaults to FALSE.
#' @import ggplot2
#' @import graphics
#' @return A ggplot2 object is either returned or directly plotted
#' @export
#' @examples
#'   gene <- c('CD4','CD14', "AIF1", "ACVR1", "CDY2A")
#'   terms_of_interest <- c("CD4 T cell", "CD14+ Monocyte", "B cell", "CD8 T cell",
#'                           "FCGR3A+ Monocyte", "NK cell", "Dendritic cell", "Megakaryocyte", 'immunity')
#'   literature_counts <- get_literature_score(gene, terms_of_interest)
#'   P <-plot_literature_score(literature_counts, return_ggplot = TRUE)
#'   plot(P)

plot_literature_score <-
  function(plot_counts,
           return_ggplot = FALSE,
           is_plotly = FALSE) {
    plot_counts$breaks <-
      cut(
        plot_counts[, 3],
        breaks = c(-0.01, 0.01, 10, 50, 100, 500, Inf),
        right = FALSE
      )
    breaks = c("[-0.01,0.01)",
               "[0.01,10)",
               "[10,50)",
               "[50,100)",
               "[100,500)",
               "[500,Inf)")
    labels = c("0", "1-10", "11-50", "51-100", "100-500", "500+")
    plot_counts$labels <- labels[match(plot_counts$breaks, breaks)]
    plot_counts$number_of_articles <- plot_counts[, 3]
    
    lbs <- c("0", "1-10", "11-50", "51-100", "100-500", "500+")[which(breaks %in% levels(factor(plot_counts$breaks)))]
    vls <-c(
      "0" = "black",
      "1-10" = "rosybrown2",
      "11-50" = "lightsalmon1",
      "51-100" = "salmon2",
      "100-500" = "indianred3",
      "500+" = 'red3'
    )[which(breaks %in% levels(factor(plot_counts$breaks)))]
    p <-
      ggplot(plot_counts, aes(Genes, Topic, fill = labels, label = number_of_articles)) +
      geom_tile(aes(fill = labels)) +
      theme(
        panel.background = element_rect(
          fill = "gray",
          colour = "gray",
          size = 0.5,
          linetype = "solid"
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_fill_manual(
        values = vls,
        name = "Article counts"
      ) 
    if (return_ggplot) {
      return(p)
    } else{
      if (is_plotly) {
        plotly::ggplotly(p)
      }
      if (!is_plotly) {
        plot(p)
      }
    }
  }

# lubianat 28/09/2018

#' literature_score_plot
#'
#' Plots a non-clusterized heatmap of the article counts for the combination of
#' gene list and list of terms
#'
#' @param plot_counts The dataframe returned as the second object ($counts) in the
#' list output of get_literature_score function
#' @param return_ggplot If TRUE, returns a ggplot2 object instead of plotting. Defaults to FALSE.
#' @param is_plotly If TRUE, a interactive plot is plotted in the place o static ggplot. Defaults to FALSE.
#' @export
#' @examples
#'gene <- 'CD4'
#'terms_of_interest <- c("CD4 T cell", "CD14+ Monocyte", "B cell", "CD8 T cell",
#'                       "FCGR3A+ Monocyte", "NK cell", "Dendritic cell", "Megakaryocyte", 'immunity')
#' literature_list <- get_literature_score(gene, terms_of_interest, max.score = 500)
#' literature_score_plot(literature_score_list$counts)




literature_score_plot <- function(plot_counts, return_ggplot=F, is_plotly = F){
  require(ggplot2)
  require(plotly)
  p <-  ggplot(plot_counts, aes(Var1, Var2)) +
    geom_tile(aes(fill = get(colnames(plot_counts)[3])), alpha = 0.6) +
    #   geom_text(aes(label = round(count, 1)))+
    theme(
      panel.background = element_rect(fill = "gray",
                                      colour = "gray",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_gradient(low = "white", high = "red")
  if (return_ggplot){
    return(p)
  } else{
  if (is_plotly){
    ggplotly(p)
  }
  if (!is_plotly){
    plot(p)
  }
  }
}

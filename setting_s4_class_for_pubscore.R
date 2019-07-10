#' @import rentrez
#' @import progress
library(dplyr)

setOldClass('gg')
setOldClass('ggplot')
setOldClass('gtable')

#' An S4 class to represent PubScore results
#' @slot terms_of_interest A list of terms of interest related to the topic you want to find the relevance.
#' @slot genes The genes to which you want to calculate and visualize the literature score.
#' @slot max_score The max score to use when calculating literature scores
#' @slot date The date when the object was initialized. PubScore counts will likely increase with time. 
#' @slot counts A data frame with the counts retrieved on PubMed
#' @slot network A visualization of the results found in a network
#' @slot heatmap A visualization of the results found in a heatmap
setClass(Class = 'PubScore', 
         slots = list(genes="character",
                      terms_of_interest = "character",
                      max_score = "numeric",
                      date = "Date",
                      counts = "data.frame",
                      network = 'gg',
                      heatmap = 'gg'))


setMethod('initialize', signature('PubScore'),
          function(.Object, genes, terms_of_interest){
            cat("~~~ Initializing PubScore ~~~ \n")
            .Object@genes <- genes
            .Object@terms_of_interest <- terms_of_interest  
             cts <- get_literature_score(genes, terms_of_interest)
            .Object@counts <- cts
            .Object@date <- Sys.Date()
            .Object@max_score <- Inf
            .Object@heatmap <- plot_literature_score(cts,
                       return_ggplot = TRUE,
                       is_plotly = FALSE)
            .Object@network <- plot_literature_graph(cts,
                         name = 'PubScore Network',
                         color = "#B30000FF",
                         n = 10)
            return(.Object)
            
          })
pub <- new(Class = "PubScore",genes = c('cd4','cd8'),terms_of_interest = c('blabla','immunity'))

#' Full PubScore analysis for visualization
#' @param terms_of_interest A list of terms of interest related to the topic you want to find the relevance for
#' @param genes A vector with multiple genes.
#' 
#' @return Object of class \code{PubScore}
#' @export

pubscore <- function(terms_of_interest, genes){
   results <- new(Class = "PubScore", genes, terms_of_interest)
   return(results)
}

#' Retrieve the heatmap attribute
#' @param pub Object of class \code{PubScore}
#' @return A "gg" object, from ggplot2, containing a heatmap from the counts table.
#' @examples
#' Create a new pubscore object
#' pub <- pubscore(genes = c('cd4','cd8'),terms_of_interest = c('blabla','immunity'))
#' plot(heatmapViz(pub))

#' @rdname heatmapViz
#' @export
setGeneric("heatmapViz", function(pub) {
  standardGeneric("heatmapViz")
})

#' @rdname heatmapViz
setMethod("heatmapViz", signature("PubScore"),
          function(pub){
            return(pub@heatmap)
          })


#' Retrieve the network attribute
#' @param pub Object of class \code{PubScore}
#' @return A "gg" object, from ggplot2, containing a network from the counts table.
#' @examples
#' Create a new pubscore object
#' pub <- pubscore(genes = c('cd4','cd8'),terms_of_interest = c('blabla','immunity'))
#' plot(networkViz(pub))

#' @rdname networkViz
#' @export
setGeneric("networkViz", function(pub) {
  standardGeneric("networkViz")
})

#' @rdname networkViz
setMethod("networkViz", signature("PubScore"),
          function(pub){
            return(pub@network)
          })




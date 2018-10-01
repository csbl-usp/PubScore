# Pairwise search on pubmed for detection of relevance to a certain topic
# lubianat 28/09/2018


#' iscore.gene
#'
#' Calculates the importance score for a given gene.
#' The importance score is the total counts of articles in the pubmed
#' database that are a result for that gene AND each term of a list
#'
#' @param terms_of_interest A list of terms of interest related to the topic you want to find the relevance for
#' @param gene The gene which you want to calculate the iscore for
#' @export
#' @examples
gene <- 'CD4'
terms_of_interest <- c("CD4 T cell", "CD14+ Monocyte", "B cell", "CD8 T cell",
                       "FCGR3A+ Monocyte", "NK cell", "Dendritic cell", "Megakaryocyte", 'immunity')


require(data.table)
require(rentrez)
all_combinations <- expand.grid(gene, terms_of_interest)
all_combinations$count <- -1


get_iscore<- function(love=TRUE){
  if(love==TRUE){
    print("I love cats!")
  }
  else {
    print("I am not a cool person.")
  }
}

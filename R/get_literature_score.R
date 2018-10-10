# Pairwise search on pubmed for detection of relevance to a certain topic
# lubianat 28/09/2018

#' get_literature_score
#'
#' Calculates the importance score for a given gene.
#' The importance score is the total counts of articles in the pubmed
#' database that are a result for that gene AND each term of a list
#'
#' @param terms_of_interest A list of terms of interest related to the topic you want to find the relevance for
#' @param gene The gene which you want to calculate the literature_score for, or a vector with multiple genes
#' @param max_score A cutoff for maximum numbers of search. Useful for avoiding outlier relations
#'  weighting too much. Defaults to 500.
#' @param wait_time How long should be the interval between two requests to the ENTREZ database when it fails.
#' Defaults to 0. In seconds.
#' @param is.list If you are searching a single gene or a list of genes. Defaults to false (single gene)
#' @export
#' @examples
#'gene <- 'CD4'
#'terms_of_interest <- c("CD4 T cell", "CD14+ Monocyte", "B cell", "CD8 T cell",
#'                       "FCGR3A+ Monocyte", "NK cell", "Dendritic cell", "Megakaryocyte", 'immunity')
#'get_literature_score(gene, terms_of_interest, max_score = 500)
#'get_literature_score(gene, terms_of_interest, max_score = Inf)




get_literature_score<- function(gene, terms_of_interest, is.list = F, max_score = 500, verbose = T,wait_time =0){
  require(data.table)
  require(rentrez)
  all_combinations <- expand.grid(gene, terms_of_interest)
  all_combinations$count <- -1
  for (index in 1:nrow(all_combinations)){

    print(index)
    search_topic <- paste0(all_combinations$Var1[index], ' AND ', all_combinations$Var2[index])
    tryCatch(
      s <- entrez_search(db = "pubmed",
                         term = search_topic,
                         retmax = 3,
                         use_history = T
                         ),

      error=function(e){
        tryCatch(
          print("Query failed, but I'm trying again"),
          Sys.sleep(wait_time),
          s <- entrez_search(db = "pubmed",
                             term = search_topic,
                             retmax = 3,
                             use_history = T
                             ),

          error=function(e){
            print("Query failed! Not your lucky day, but I'm trying again")
            Sys.sleep(wait_time)
            s <- entrez_search(db = "pubmed",
                               term = search_topic,
                               retmax = 3)
            print("Actually, you are safe. ")
          })
      })
    all_combinations$count[index] <- s$count
  }
  all_combinations$count[all_combinations$count>max_score]<-max_score
  colnames(all_combinations)[3] <- paste0('counts_', Sys.Date())
  if (!is.list) {
    gene_literature_score <- sum(all_combinations$count)
    literature_score <- list(gene_literature_score = gene_literature_score, counts = all_combinations)
    return(literature_score)
  }
  if (is.list) {
    list_literature_score <- sum(all_combinations$count)/length(gene)
    literature_score <- list(list_literature_score = list_literature_score, counts = all_combinations)
    return(literature_score)

    return(all_combinations)
  }

}

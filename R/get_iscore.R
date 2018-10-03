# Pairwise search on pubmed for detection of relevance to a certain topic
# lubianat 28/09/2018

#' get_iscore
#'
#' Calculates the importance score for a given gene.
#' The importance score is the total counts of articles in the pubmed
#' database that are a result for that gene AND each term of a list
#'
#' @param terms_of_interest A list of terms of interest related to the topic you want to find the relevance for
#' @param gene The gene which you want to calculate the iscore for, or a vector with multiple genes
#' @param max.score A cutoff for maximum numbers of search. Useful for avoiding outlier relations
#'  weighting too much. Defaults to 500.
#' @param is.list If you are searching a single gene or a list of genes. Defaults to false (single gene)
#' @export
#' @examples
#'gene <- 'CD4'
#'terms_of_interest <- c("CD4 T cell", "CD14+ Monocyte", "B cell", "CD8 T cell",
#'                       "FCGR3A+ Monocyte", "NK cell", "Dendritic cell", "Megakaryocyte", 'immunity')
#' get_iscore(gene, terms_of_interest, max.score = 500)
#'get_iscore(gene, terms_of_interest, max.score = Inf)




get_iscore<- function(gene, terms_of_interest, is.list = F, max.score = 500, verbose = T){
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
          Sys.sleep(1.2),
          s <- entrez_search(db = "pubmed",
                             term = search_topic,
                             retmax = 3,
                             use_history = T
                             ),

          error=function(e){
            print("Query failed! Not your lucky day, but I'm trying again")
            Sys.sleep(1)
            s <- entrez_search(db = "pubmed",
                               term = search_topic,
                               retmax = 3)
            print("Actually, you are safe. ")
          })
      })
    all_combinations$count[index] <- s$count
  }
  all_combinations$count[all_combinations$count>max.score]<-max.score
  colnames(all_combinations)[3] <- paste0('counts_', Sys.Date())
  if (!is.list) {
    gene_iscore <- sum(all_combinations$count)
    iscore <- list(gene_iscore = gene_iscore, counts = all_combinations)
    return(iscore)
  }
  if (is.list) {
    list_iscore <- sum(all_combinations$count)/length(gene)
    iscore <- list(list_iscore = list_iscore, counts = all_combinations)
    return(iscore)

    return(all_combinations)
  }

}

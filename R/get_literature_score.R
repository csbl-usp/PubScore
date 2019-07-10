

#' .query_pubmed
#'
#' Auxiliary function for getting the list score
#' @param search_topic Item to search on PubMed via rentrez
#' @param wait_time Time between searches
.query_pubmed <- function (search_topic, wait_time = 0) {
  out <- tryCatch({
    s <- rentrez::entrez_search(
      db = "pubmed",
      term = search_topic,
      retmax = 1,
      use_history = TRUE
    )
    return(s)
  },
  
  error = function(e) {
    tryCatch({
      print("Query failed, but I'm trying again")
      Sys.sleep(wait_time)
      s <- rentrez::entrez_search(
        db = "pubmed",
        term = search_topic,
        retmax = 1,
        use_history = TRUE
      )
      return(s)
      
    },
    
    error = function(e) {
      print("Query failed! Not your lucky day, but I'm trying again")
      Sys.sleep(wait_time)
      s <- rentrez::entrez_search(db = "pubmed",
                                  term = search_topic,
                                  retmax = 1)
      print("Actually, you are safe.")
      return(s)
    })
  })
  return(out)
}

# Pairwise search on pubmed for detection of relevance to a certain topic
# lubianat 28/09/2018

#' get_literature_score
#'
#' Calculates the importance score for a given gene.
#' The importance score is the total counts of articles in the pubmed
#' database that are a result for that gene AND each term of a list
#'
#' @param terms_of_interest A list of terms of interest related to the topic you want to find the relevance for
#' @param genes A vector with multiple genes.
#' @param wait_time How long should be the interval between two requests to the ENTREZ database when it fails.
#' Defaults to 0. In seconds.
#' @param show_progress If TRUE, a progress bar is displayed. Defaults to TRUE.
#' @param verbose If TRUE, will display the index of the search occuring. Defaults to FALSE.
#' @import rentrez
#' @import progress
#' @return A dataframe with the literature scores.
#' @export
#' @examples
#' genes <- c('CD8', 'CD4')
#' terms_of_interest <-
#'   c(
#'     "CD4 T cell",
#'     "CD14+ Monocyte",
#'     "B cell",
#'     "CD8 T cell",
#'     "FCGR3A+ Monocyte",
#'     "NK cell",
#'     "Dendritic cell",
#'     "Megakaryocyte",
#'     'immunity'
#'   )
#' get_literature_score(genes, terms_of_interest)

get_literature_score <-
  function(genes,
           terms_of_interest,
           wait_time = 0,
           show_progress = TRUE,
           verbose = FALSE) {
    all_combinations <- expand.grid(genes, terms_of_interest)
    all_combinations$count <- -1
    colnames(all_combinations) <- c('Genes', 'Topic', 'count')
    number_of_combinations <- nrow(all_combinations)
    # This shows a progress bar to the user.
    pb <-
      progress::progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", total = number_of_combinations)
    if (show_progress) {
      message('Running PubScore:')
    }
    for (index in seq_len(nrow(all_combinations))) {
      if (show_progress) {
        pb$tick()
      }
      if (verbose) {
        print(index)
      }
      search_topic <-
        paste0(all_combinations$Genes[index], ' AND ', all_combinations$Topic[index])
      
      s <- .query_pubmed(search_topic)
      all_combinations$count[index] <- s$count
    }
    return(all_combinations)
  }

.getSimulation_test <- function(literature_object, ambiguous = c()) {
  simulation_of_literature_null <-
    literature_object$all_gene_combinations[!literature_object$all_gene_combinations$Gene %in% ambiguous, ]
  
  n_genes <- length(levels(literature_object$counts$Gene[!literature_object$counts$Gene %in% ambiguous ]))
  total_genes <- levels(droplevels(simulation_of_literature_null$Gene))
  message('Running 100000 simulations')
  
  distribution_of_scores <- c()
  for (i in seq_len(100000)) {
    genes_to_sample_now <- sample(total_genes, n_genes)
    simu_now <-
      simulation_of_literature_null[simulation_of_literature_null$Gene %in% genes_to_sample_now, ]$count
    list_score <- sum(simu_now) / n_genes
    distribution_of_scores <- c(distribution_of_scores, list_score)
    
  }
  
  distribution_of_scores <- data.frame(d = distribution_of_scores)
  return(distribution_of_scores)
  
}

#  Estimates p-value for a literature score provided an object 
#from the get_literature_score function
# lubianat 28/09/2018

#' test_literature_score
#'
#' gets a pvalue, by simulation, for a literature score.
#'
#' @param literature_object The object returned by the get_literature_score function
#' @param total_genes A list of all the possible genes in your study. 
#' Usually all the names in the rows of an "exprs" object.
#' @param show_progress If TRUE, a progress bar is displayed. Defaults to True.
#' @param verbose If TRUE, will display the index of the search occuring. Defaults to false.
#' @param remove_ambiguous If TRUE, ambiguously named genes (such as "MARCH") will be removed. Defaults to TRUE.  
#' @import rentrez
#' @import progress
#' @return A dataframe with the literature scores.
#' @export
#' @examples
#'      genes <- c('CD14')
#'      terms_of_interest <-
#'        c(
#'          "CD4 T cell",
#'          "CD14+ Monocyte",
#'          "B cell",
#'          "CD8 T cell",
#'          "FCGR3A+ Monocyte",
#'          "NK cell",
#'          "Dendritic cell",
#'          "Megakaryocyte",
#'          'immunity'
#'        )
#'    literature_object <- PubScore::get_literature_score(genes, terms_of_interest)
#'    total_genes <- c('CD14', 'notagenemock1','notagenemock2','notagenemock3', "PROC", "IMPACT")
#'    literature_object <- test_literature_score(literature_object, total_genes)
#'  
test_literature_score <-
  function(literature_object,
           total_genes,
           show_progress = TRUE,
           remove_ambiguous = TRUE,
           verbose = FALSE) {
    terms_of_interest <- levels(literature_object$counts$Topic)
    genes_to_sample <- length(levels(literature_object$counts$Gene))
    
    simulation_of_literature_null <- data.frame(2,2,2)
    simulation_of_literature_null <- simulation_of_literature_null[-1,]

    message('Running PubScore for all genes. Might take a while!')
    pb_test_score <-
      progress::progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", total = length(total_genes))
    
    for (i in total_genes){
      pb_test_score$tick()

      tryCatch({
        ps_object <- PubScore::get_literature_score(i, terms_of_interest, 
                                                    max_score =  literature_object$max_score,
                                                    show_progress = FALSE)
        new_line <- data.frame(i, ps_object$counts)
        simulation_of_literature_null <- rbind(simulation_of_literature_null, new_line)
        Sys.sleep(0.2)}, error = function(e){print(e)
        })
    }
      simulation_of_literature_null <- simulation_of_literature_null[simulation_of_literature_null$i != "",]  
      
      obj <- literature_object
      obj[['all_gene_combinations']] = simulation_of_literature_null
      

    if (remove_ambiguous == FALSE){
      
      distribution_of_scores <- .getSimulation_test(obj)
      score <- mean(literature_object$counts$count)
      pvalue <-sum(distribution_of_scores[,1] >= score)/length(distribution_of_scores[,1])
      
      print('The p-value by simulation is:')
      print(pvalue) 
      
      obj[['p_value']] = pvalue
      return(obj)
    }


    if (remove_ambiguous == TRUE){
    ambiguous_terms <- c("PC", "JUN", "IMPACT", "ACHE", "SRI", "SET", "CS", "PROC", 
                         "MET", "SHE", "CAD", "DDT", "PIGS", "SARS", "REST", "GC", "CP", 
                         "STAR", "SI", "GAN", "MARS", "SDS", "AGA", "NHS", "CPE", "POR", 
                         "MAX", "CAT", "LUM", "ANG", "POLE", "CLOCK", "TANK", "ITCH", 
                         "SDS", "AES", "CIC", "FST", "CAPS", "COPE", "F2", "AFM", "SPR", 
                         "PALM", "C2", "BAD", "GPI", "CA2", "SMS", "INVS", "WARS", "HP", 
                         "GAL", "SON", "AFM", "BORA", "MBP", "MAK", "MALL", "COIL", "CAST ")
    
    distribution_of_scores <- .getSimulation_test(obj, ambiguous =  ambiguous_terms)
    
    score <- mean(literature_object$counts$count[!literature_object$counts$Gene %in% ambiguous_terms])
    pvalue <-sum(distribution_of_scores[,1] >= score)/length(distribution_of_scores[,1])
    
    print('The p-value by simulation removing ambiguous gene names is:')
    print(pvalue) 
    
    obj[['p_value']] = pvalue
    return(obj)
    }
  }
    

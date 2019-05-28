

.getSimulation_retest <- function(literature_object, ambiguous = c(),  n_simulations ) {
  simulation_of_literature_null <-
    literature_object$all_gene_combinations[!literature_object$all_gene_combinations$Gene %in% ambiguous, ]
  n_genes <- length(levels(literature_object$counts$Gene))
  total_genes <- simulation_of_literature_null$Gene
  message(paste0('Running ', n_simulations,' simulations'))
  distribution_of_scores <- c()
  for (i in seq_len(n_simulations)) {
    genes_to_sample_now <- sample(total_genes, n_genes)
    simu_now <-
      simulation_of_literature_null[simulation_of_literature_null$Gene %in% genes_to_sample_now, ]$count
    list_score <- sum(simu_now) / n_genes
    distribution_of_scores <- c(distribution_of_scores, list_score)
    
  }
  
  distribution_of_scores <- data.frame(d = distribution_of_scores)
  return(distribution_of_scores)
  
}

#  Re-Estimates p-value for a literature score provided an object 
# from the get_literature_score function
# lubianat 24/05/2019

#' retest_literature_score
#'
#' Regets a pvalue, by simulation, for a literature score.
#'
#' @param literature_object The object returned by the test_literature_score function
#' @param new_ambigous_term_list A vector containing a list of genes deemed ambiguous by the user. 
#' If left blank, a default list is used.
#' @param remove_ambiguous If TRUE, ambiguously named genes (such as "MARCH") will be removed. 
#' Defaults to TRUE.  
#' @param nsim The number of simulations to run. Defaults to 100000.
#' @import rentrez
#' @import progress
#' @return A dataframe with the literature scores.
#' @export
#' @examples
#' 
#' genes <- c('CD14')
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
#' literature_object <- PubScore::get_literature_score(genes, terms_of_interest)
#' total_genes <- c('CD14','notagenemock1','notagenemock2','notagenemock3', "PROC", "IMPACT")
#' literature_object <- test_literature_score(literature_object, total_genes, remove_ambiguous = FALSE)
#' literature_object <- retest_literature_score(literature_object, remove_ambiguous = TRUE)
#' 
retest_literature_score <-
  function(literature_object,
           new_ambigous_term_list = FALSE,
           remove_ambiguous = FALSE,
           nsim = 100000) {
    genes_to_sample <- length(levels(literature_object$counts$Gene))
    
    
    if (any(new_ambigous_term_list != FALSE)) {
      distribution_of_scores <-
        .getSimulation_retest(literature_object, new_ambigous_term_list, n_simulations = nsim)
      score <-
        sum(literature_object$counts$count[!literature_object$counts$Gene %in% new_ambigous_term_list])/genes_to_sample
      pvalue <-
        sum(distribution_of_scores[, 1] >= score) / length(distribution_of_scores[, 1])
      literature_object[['p_value']] <- pvalue
      print('The p-value by simulation removing the new list of ambiguous genes is:')
      print(pvalue) 
      
      return(literature_object)
    }
    
    if (remove_ambiguous == FALSE) {
      distribution_of_scores <- .getSimulation_retest(literature_object,  n_simulations = nsim)
      
      score <- sum(literature_object$counts$count)/genes_to_sample
      pvalue <-
        sum(distribution_of_scores[, 1] >= score) / length(distribution_of_scores[, 1])
      
      print('The p-value by simulation  is:')
      print(pvalue) 
      
      literature_object[['p_value']] <- pvalue
      
      return(literature_object)
      
      
      
    }
    
    
    
    if (remove_ambiguous == TRUE) {
      ambiguous_terms <-
        c(
          "PC",
          "JUN",
          "MARCH",
          "IMPACT",
          "ACHE",
          "SRI",
          "SET",
          "CS",
          "PROC",
          "MET",
          "SHE",
          "CAD",
          "DDT",
          "PIGS",
          "SARS",
          "REST",
          "GC",
          "CP",
          "STAR",
          "SI",
          "GAN",
          "MARS",
          "SDS",
          "AGA",
          "NHS",
          "CPE",
          "POR",
          "MAX",
          "CAT",
          "LUM",
          "ANG",
          "POLE",
          "CLOCK",
          "TANK",
          "ITCH",
          "SDS",
          "AES",
          "CIC",
          "FST",
          "CAPS",
          "COPE",
          "F2",
          "AFM",
          "SPR",
          "PALM",
          "C2",
          "BAD",
          "GPI",
          "CA2",
          "SMS",
          "INVS",
          "WARS",
          "HP",
          "GAL",
          "SON",
          "AFM",
          "BORA",
          "MBP",
          "MAK",
          "MALL",
          "COIL",
          "CAST "
        )
      distribution_of_scores2 <-
        .getSimulation_retest(literature_object, ambiguous_terms,  n_simulations = nsim)
      
      score <-
        sum(literature_object$counts$count[!literature_object$counts$Gene %in% ambiguous_terms])/genes_to_sample
      pvalue <-
        sum(distribution_of_scores2[, 1] >= score) / length(distribution_of_scores2[, 1])
      literature_object[['p_value']] <- pvalue
      
      print('The p-value by simulation after removing the default list of ambiguous genes is:')
      print(pvalue) 
      
      return(literature_object)
    }
  }


#  Estimates p-value for a literature score provided an object 
#from the get_literature_score function
# lubianat 28/09/2018

#' get_literature_score
#'
#' Calculates the importance score for a given gene.
#' The importance score is the total counts of articles in the pubmed
#' database that are a result for that gene AND each term of a list
#'
#' @param literature_object The object returned by the get_literature_score function
#' @param total_genes A list of all the possible genes in your study. 
#' Usually all the names in the rows of an "exprs" object.
#' @param show_progress If TRUE, a progress bar is displayed. Defaults to True.
#' @param verbose If TRUE, will display the index of the search occuring. Defaults to false.
#' @import rentrez
#' @import progress
#' @return A dataframe with the literature scores.
#' @export
#' @examples
#'    genes <- c('CD4','CD8')
#'    terms_of_interest <-
#'      c(
#'        "CD4 T cell",
#'        "CD14+ Monocyte",
#'        "B cell",
#'        "CD8 T cell",
#'        "FCGR3A+ Monocyte",
#'        "NK cell",
#'        "Dendritic cell",
#'        "Megakaryocyte",
#'        'immunity'
#'      )
#'  literature_object <- PubScore::get_literature_score(gene, terms_of_interest, 500)
#'  total_genes <- c('CD4', 'CD8','CD7','ASDASD','ADSFASD','ADSFASDS', 'asdas')
#'  literature_object <- test_literature_score(literature_object, total_genes)

test_literature_score <-
  function(literature_object,
           total_genes,
           show_progress = TRUE,
           verbose = FALSE) {
    terms_of_interest <- levels(literature_object$counts$Var2)
    genes_to_sample <- length(levels(literature_object$counts$Var1))
    
    simulation_of_literature_null <- data.frame(2,2,2)
    simulation_of_literature_null <- simulation_of_literature_null[-1,]

    message('Running PubScore for all genes. Might take a while!')
    for (i in total_genes){
      pb_test_score$tick()

      tryCatch({
        ps_object <- PubScore::get_literature_score(i, terms_of_interest, 
                                                    max_score =  literature_object$max_score,
                                                    show_progress = F)
        new_line <- data.frame(i, ps_object$counts)
        simulation_of_literature_null <- rbind(simulation_of_literature_null, new_line)
        Sys.sleep(0.2)}, error = function(e){print(e)
        })
    }
      simulation_of_literature_null <- simulation_of_literature_null[simulation_of_literature_null$i != "",]  
      
      message('Running 100000 simulations')
      distribution_of_scores <- c()
      for (i in 1:100000){    
        genes_to_sample_now <- sample(total_genes, genes_to_sample)
        simu_now <- simulation_of_literature_null[simulation_of_literature_null$Var1 %in% genes_to_sample_now,]$count
        list_score <- sum(simu_now) / genes_to_sample
        distribution_of_scores <- c(distribution_of_scores,list_score )
        
      }
      
      distribution_of_scores<- data.frame(d = distribution_of_scores)
      
      
      library(ggplot2)
      ggplot(distribution_of_scores, aes(x = d))+
        geom_density() +
        ggtitle('Density of the distribution of null scores(10^6 resamplings)') +
        xlab('Average number of counts on PubMed') +
        ylab('Proportion of gene-lists for each category')
      
      
      score <- literature_object$list_literature_score
      

      pvalue <-sum(distribution_of_scores[,1] >= score)/length(distribution_of_scores[,1])
    
    
    print('The p-value by simulation is:')
    print(pvalue) 
    obj <- literature_object
    obj[['all_gene_combinations']] = simulation_of_literature_null
    obj[['p_value']] = pvalue
    tested_score <- obj
    return(tested_score)

    
  }

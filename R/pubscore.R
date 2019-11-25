#' @import rentrez
#' @import progress
#' @import dplyr
#' @importFrom methods new
#' @importFrom utils data

setOldClass('gg')
setOldClass('ggplot')
setOldClass('gtable')

####### SECTION 1 - THE PUBSCORE CLASS  #######


#' An S4 class to represent PubScore results
#'
#' The S4 class to PubScore and its basic initializate and show methods.
#' @slot terms_of_interest A list of terms of interest related to the topic you want to find the relevance.
#' @slot genes The genes to which you want to calculate and visualize the literature score.
#' @slot date The date when the object was initialized. PubScore counts will likely increase with time.
#' @slot gene2pubmed Boolean noting if gene to pubmed was used or not.
#' @slot counts A data frame with the counts retrieved on PubMed
#' @slot network A visualization of the results found in a network
#' @slot heatmap A visualization of the results found in a heatmap
setClass(
  Class = 'PubScore',
  slots = list(
    genes = "character",
    terms_of_interest = "character",
    literature_score = "numeric",
    date = "Date",
    counts = "data.frame",
    network = 'gg',
    heatmap = 'gg',
    all_counts = "data.frame",
    total_genes = "character",
    p_value = "numeric",
    gene2pubmed = "logical"
  )
)

#'initialize
#'
#'@param gene2pubmed Boolean defining if gene2pubmed db is going to be used.
#'@return A object of the PubScore class
setMethod('initialize', signature('PubScore'),
          function(.Object, genes, terms_of_interest, gene2pubmed = FALSE) {
            cat("~~~ Initializing PubScore Object ~~~ \n")
            .Object@genes <- genes
            .Object@terms_of_interest <- terms_of_interest
            
            if (gene2pubmed){
              cts <- get_literature_score(genes, terms_of_interest, gene2pubmed = gene2pubmed)
              all_counts <- get_literature_score(genes, terms_of_interest, gene2pubmed = gene2pubmed, return_all = TRUE)
              
              .Object@all_counts <- all_counts
              
            } else{
              cts <- get_literature_score(genes, terms_of_interest)
              .Object@counts <- cts
              .Object@all_counts <- data.frame()
            }

            .Object@date <- Sys.Date()
            .Object@heatmap <- plot_literature_score(cts,
                                                     return_ggplot = TRUE,
                                                     is_plotly = FALSE)
            .Object@network <- plot_literature_graph(cts,
                                                     name = 'PubScore Network',
                                                     color = "#B30000FF",
                                                     max_number_of_labels = 10)
            .Object@literature_score <-
              sum(cts$count) / (length(genes) * length(terms_of_interest))
            .Object@total_genes <- 'empty'
            .Object@gene2pubmed <- gene2pubmed
            .Object@p_value <- Inf
            return(.Object)
            
          })

setMethod(
  'show',
  signature = 'PubScore',
  definition = function(object) {
    cat("An object of class ", class(object), "\n", sep = "")
    cat(" - Number of genes: ", length(object@genes), "\n", sep = "")
    if (length(object@genes) < 6) {
      cat(object@genes, sep = ", ")
      cat('\n')
    } else {
      cat(object@genes[1:5], sep = ", ")
      cat(" and other ", (length(object@genes) - 5), " more.", "\n")
    }
    cat(" - Number of terms of interest: ",
        length(object@terms_of_interest),
        "\n",
        sep = "")
    if (length(object@terms_of_interest) < 6) {
      cat(object@terms_of_interest, sep = ", ")
      cat('\n')
    } else {
      cat(object@terms_of_interest[1:5], sep = ", ")
      cat(" and other ", (length(object@terms_of_interest) - 5), " more.", "\n")
    }
    if (nrow(pub@all_counts) == 0) {
      cat(
        "- P-value for association: ",
        " Not calculated. Run 'test_score' to get one.",
        "\n",
        sep = ""
      )
    } else {
      cat(
        "- P-value for association: ",
        object@p_value ,
        ", estimated by simulation.",
        "\n",
        sep = ""
      )
    }
  }
)

#' PubScore fundamental analysis
#'
#' Runs the initalization and the basic functions for querying pubmed and getting the literature scores.
#'
#' @param terms_of_interest A list of terms of interest related to the topic you want to find the relevance for
#' @param genes A vector with multiple genes.
#' @param gene2pubmed Boolean defining if gene2pubmed db is going to be used. Defaults to FALSE.
#' @return Object of class \code{PubScore}
#' @export

pubscore <- function(terms_of_interest, genes, gene2pubmed = FALSE) {
  results <- new(Class = "PubScore", genes, terms_of_interest, gene2pubmed)
  return(results)
}


####### SECTION 2 - METHODS TO RETRIEVE ATTRIBUTES  #######


#' Retrieve the heatmap attribute
#' 
#' @param pub Object of class \code{PubScore}
#' @return A "gg" object, from ggplot2, containing a heatmap from the counts table.
#' @examples
#' #Create a new pubscore object
#' pub <- pubscore(genes = c('cd4','cd8'),terms_of_interest = c('blabla','immunity'))
#' plot(heatmapViz(pub))
#' @rdname heatmapViz
#' @export
#' 
setGeneric("heatmapViz", function(pub) {
  standardGeneric("heatmapViz")
})

#' @rdname heatmapViz
setMethod("heatmapViz", signature("PubScore"),
          function(pub) {
            return(pub@heatmap)
          })


#' Retrieve the network attribute
#' 
#' @param pub Object of class \code{PubScore}
#' @return A "gg" object, from ggplot2, containing a network from the counts table.
#' @examples
#' # Create a new pubscore object
#' pub <- pubscore(genes = c('cd4','cd8'),terms_of_interest = c('blabla','immunity'))
#' plot(networkViz(pub))
#' @rdname networkViz
#' @export
#' 
setGeneric("networkViz", function(pub) {
  standardGeneric("networkViz")
})
#' @rdname networkViz
setMethod("networkViz", signature("PubScore"),
          function(pub) {
            return(pub@network)
          })


#' Retrieve the literature_score attribute
#' 
#' @param pub Object of class \code{PubScore}
#' @return A "numeric" with the literature score for this gene x term combination
#' @examples
#' # Create a new pubscore object
#' pub <- pubscore(genes = c('cd4','cd8'),terms_of_interest = c('blabla','immunity'))
#' plot(networkViz(pub))
#' @rdname getScore
#' @export
setGeneric("getScore", function(pub) {
  standardGeneric("getScore")
})

#' @rdname getScore
setMethod("getScore", signature("PubScore"),
          function(pub) {
            return(pub@literature_score)
          })

#' Retrieve the all_counts attribute
#' 
#' @param pub Object of class \code{PubScore}
#' @return A dataframe containing the counts table for all genes.
#' @examples
#' # Create a new pubscore object
#' pub <- pubscore(genes = c('cd4','cd8'),terms_of_interest = c('blabla','immunity'))
#' plot(networkViz(pub))
#' @rdname AllCounts
#' @export
setGeneric("AllCounts", function(pub, ...) {
  standardGeneric("AllCounts")
})

#' @rdname AllCounts
setMethod("AllCounts", signature("PubScore"),
          function(pub) {
            return(pub@all_counts)
          })

#' Auxiliary function for the test method
#' @param pub An object of class \code{PubScore}
#' @param ambiguous A character vector with possible ambiguous gene names
#' @param n_simulations The number of simulations to run.
#' @return A data-frame with a simulation of literature scores for random samplings
.getSimulation_test <-
  function(pub, ambiguous = c(), n_simulations) {
    simulation_of_literature_null <-
      pub@all_counts[!pub@all_counts$Genes %in% ambiguous,]
    
    n_genes <- length(pub@genes[!pub@genes %in% ambiguous])
    total_genes <-
      levels(droplevels(simulation_of_literature_null$Genes))
    message(paste0('Running ', n_simulations, ' simulations'))
    distribution_of_scores <- c()
    
    for (i in seq_len(n_simulations)) {
      genes_to_sample_now <- sample(total_genes, n_genes)
      simu_now <-
        simulation_of_literature_null[simulation_of_literature_null$Genes %in% genes_to_sample_now,]$count
      list_score <-
        sum(simu_now) / (length(pub@genes) * length(pub@terms_of_interest))
      distribution_of_scores <- c(distribution_of_scores, list_score)
      
    }
    
    distribution_of_scores <- data.frame(d = distribution_of_scores)
    return(distribution_of_scores)
    
  }

#######  SECTION 3 - ADDITIONAL METHODS ####### 

#' Test the literature enrichment score
#' 
#' @param pub Object of class \code{PubScore}
#' @param total_genes A list of all the possible genes in your study.
#' Usually all the names in the rows of an "exprs" object.
#' @param show_progress If TRUE, a progress bar is displayed. Defaults to True.
#' @param verbose If TRUE, will display the index of the search occuring. Defaults to false.
#' @param remove_ambiguous If TRUE, ambiguously named genes (such as "MARCH") will be removed. Defaults to TRUE.
#' @param nsim The number of simulations to run. Defaults to 100000.
#' @param ambiguous_terms A character vector of the ambiguous terms to use instead of the default. 
#' The default includes 60 genes pre-selected as ambiguous (as IMPACT, MARCH and ACHE).
#' @return A "gg" object, from ggplot2, containing a network from the counts table.
#' @examples
#' # Create a new pubscore object
#' pub <- pubscore(genes = c('cd4','cd8'),terms_of_interest = c('blabla','immunity'))
#' pub <- test_score(pub, total_genes = c('notagene1', 'notagene2', 'cd4', 'cd8'),remove_ambiguous = TRUE)
#' @rdname test_score
#' @export
setGeneric("test_score", function(pub,
                                  total_genes,
                                  show_progress = TRUE,
                                  remove_ambiguous = TRUE,
                                  verbose = FALSE,
                                  nsim = 100000,
                                  ambiguous_terms = c(
                                    "PC",
                                    "JUN",
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
                                  )) {
  standardGeneric("test_score")
})

#' @rdname test_score
setMethod("test_score", signature("PubScore"),
          function(pub,
                   total_genes,
                   show_progress = TRUE,
                   remove_ambiguous = TRUE,
                   verbose = FALSE,
                   nsim = 100000,
                   ambiguous_terms = c(
                     "PC",
                     "JUN",
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
                   )) {
            if (pub@gene2pubmed){
              stop("Test score not available when gene2pubmed is used... yet!")
            }
            terms_of_interest <- pub@terms_of_interest
            genes_to_sample <- pub@genes
            simulation_of_literature_null <- data.frame(2, 2, 2)
            simulation_of_literature_null <-
              simulation_of_literature_null[-1, ]
            
            
            if (length(pub@all_counts) == 0) {
              pub@total_genes <- total_genes
              message('Running the PubScore test for all genes. Might take a while to finish queries!')
              pb_test_score <-
                progress::progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", total = length(total_genes))
              
              for (i in total_genes) {
                pb_test_score$tick()
                
                tryCatch({
                  pub_counts <-
                    get_literature_score(i, terms_of_interest, show_progress = FALSE)
                  new_line <- data.frame(i, pub_counts)
                  simulation_of_literature_null <-
                    rbind(simulation_of_literature_null, new_line)
                  Sys.sleep(0.2)
                }, error = function(e) {
                  print(e)
                })
              }
              simulation_of_literature_null <-
                simulation_of_literature_null[simulation_of_literature_null$i != "", ]
              pub@all_counts <- simulation_of_literature_null[, -1]
            }
            if (remove_ambiguous == FALSE) {
              distribution_of_scores <-
                .getSimulation_test(pub, n_simulations = nsim)
              score <- pub@literature_score
              pvalue <-
                sum(distribution_of_scores[, 1] >= score) / length(distribution_of_scores[, 1])
              
              print('The p-value by simulation is:')
              print(pvalue)
              pub@p_value <- pvalue
              return(pub)
            }
            
            
            if (remove_ambiguous == TRUE) {
              distribution_of_scores <-
                .getSimulation_test(pub, ambiguous =  ambiguous_terms, n_simulations = nsim)
              score <- pub@literature_score
              pvalue <-
                sum(distribution_of_scores[, 1] >= score) / length(distribution_of_scores[, 1])
              
              print('The p-value by simulation is:')
              print(pvalue)
              pub@p_value <- pvalue
              return(pub)
            }
          })


######### Methods to insert attributes  ######### 


#' Set the all_counts attribute
#' 
#' @param pub Object of class \code{PubScore}
#' @return A dataframe containing the counts table for all genes.
#' @examples
#' terms_of_interest <- c('Dengue')
#' pub <- pubscore(terms_of_interest = terms_of_interest, genes = c("CD4", "CD8", "CD14") )
#' print(getScore(pub))
#' data("all_counts")
#' set_all_counts(pub) <- all_counts
#' @rdname set_all_counts
#' @export
setGeneric("set_all_counts<-", function(pub, value) {
  standardGeneric("set_all_counts<-")
})

#' @rdname set_all_counts
setMethod("set_all_counts<-", signature("PubScore"),
          function(pub, value) {
            pub@all_counts <- value
            return(pub)
          })

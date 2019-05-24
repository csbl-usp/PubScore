
#' literature_object after pvalue testing
#' Result of the test literature score function
#'

#' @name literature_object
#' @docType data
#' @usage data(literature_object)
#' @format A list
#' @keywords datasets, pubmed, test, literature
#' @examples
#' data(literature_object)
#' print(literature_object$p_value)
#' retest_literature_score(literature_object, remove_ambiguous = TRUE, nsim = 10000 )
"literature_object"

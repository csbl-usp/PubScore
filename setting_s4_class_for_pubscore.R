setClass(Class = 'PubScore', 
         representation = representation(genes="character", terms_of_interest = "character", max_score = "numeric"))

pubscore <- new('PubScore', genes = c('CD4', 'CD8'),  terms_of_interest = c('immunology', 'disease'))

slotNames("PubScore")

getClass('PubScore')

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

setMethod(
  f="search",
  signature="PubScore",
  definition = function(x){cat(x@genes)}
)

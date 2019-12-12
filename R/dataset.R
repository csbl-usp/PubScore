#' hgcn_entrez_reference
#' 
#' 
#' Contains the result of a query to the biomaRt service done in October, 2019. 
#' 
#' 2 columns: entrezgene_id (containing the Entrez ids) and 
#' hgnc_symbol (containing gene symbols from the HUGO gene nomenclature consortium)
#' 
#' 20491 rows, for the mapping between the two nomenclatures for human genes.
#'
#'
#' @name hgcn_entrez_reference
#' @docType data
#' @references Mapping identifiers for the integration of genomic datasets with 
#' the R/Bioconductor package biomaRt. Steffen Durinck,Paul T. Spellman, Ewan Birney and Wolfgang Huber, 
#' Nature Protocols 4, 1184-1191 (2009).
#' @usage data(hgcn_entrez_reference)
#' @format An object of class \code{data.frame}
#' @keywords datasets, pubmed, test, literature
"hgcn_entrez_reference"





#' human genes on gene2pubmed_db
#'
#' A subset of the gene2pubmed database downloaded via FTP from 
#' ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz.
#' #' The subset contains only the rows  corresponding to humans (#tax_id = 906)
#' The table was downloaded in October 2019. 
#'
#' Contains:
#' 3 columns: 
#' #tax_id: The reference ID for the taxon. All are 9606 (humans).
#' GeneID: The Entrez ID code for a given gene.
#' PubMedID: A PubMed ID for a paper that mentions the gene in the "Gene ID" column.
#' 
#' 1335548 rows: gene-paper associations in the gene2pubmed database.
#'
#'
#' @name gene2pubmed_db
#' @docType data
#' @references Maglott, Donna, et al. 'Entrez Gene: gene-centered information at NCBI.'
#'  Nucleic acids research 33.suppl_1 (2005): D54-D58.
#' @usage data(gene2pubmed_db)
#' @format An object of class \code{data.frame}
#' @keywords datasets, pubmed, test, literature
"gene2pubmed_db"



#' all_counts
#'
#' A dataframe with all pubmed counts  for the genes in the Dengue dataset in relation 
#' to the term "Dengue". 
#' 
#' Outcome of the test_score method of the pubscore class.
#' As this function may take a long time, this dataset speeds up the vignette.
#' 
#' Contains:
#' 3 columns: 
#' #tax_id: The reference ID for the taxon. All are 9606 (humans).
#' GeneID: The Entrez ID code for a given gene.
#' PubMedID: A PubMed ID for a paper that mentions the gene in the "Gene ID" column.
#' 
#' 1335548 rows: gene-paper associations in the gene2pubmed database.
#' @name all_counts
#' @docType data
#' @usage data(all_counts)
#' @format An object of class \code{data.frame}
#' @keywords datasets, pubmed, test, literature
"all_counts"

library(PubScore)
library(testthat)

test_that("get literature works", {
  selected_genes <- c('CD79A', 'CD14', 'NKG7', 'CST3', 'AIF1')
  
  # These genes were selected from a panel containing the following genes
  total_genes <- c('CD79A', 'CD14', 'NKG7', 'CST3', 'AIF1', 'FOXA1', 'PPT2', 'ZFP36L1','AFF4', 'ANTXR2', "HDAC8", "VKORC1" )
  
  terms_of_interest <- c("B cells", "macrophages", "NK cells")
  pub <- pubscore(terms_of_interest = terms_of_interest, genes = selected_genes )
  score <- getScore(pub)
  expect_that(score, is_more_than(500))
  
})

test_that("test literature  score works", {
  selected_genes <- c('CD79A', 'CD14', 'NKG7', 'CST3', 'AIF1')
  
  # These genes were selected from a panel containing the following genes
  total_genes <- c('CD79A', 'CD14', 'NKG7', 'CST3', 'AIF1', 'FOXA1', 'PPT2', 'ZFP36L1','AFF4', 'ANTXR2', "HDAC8", "VKORC1" )
  
  terms_of_interest <- c("B cells", "macrophages", "NK cells")
  pub <- pubscore(terms_of_interest = terms_of_interest, genes = selected_genes )
  pub <- test_score(pub, total_genes =  total_genes)
  expect_that(pub@p_value, is_less_than(0.01))
  
})

test_that("test literature  score works when removing ambiguous terms", {
  selected_genes <- c('CD79A', 'CD14', 'NKG7', 'CST3', 'AIF1')
  
  # These genes were selected from a panel containing the following genes
  total_genes <- c('CD79A', 'CD14', 'NKG7', 'CST3', 'AIF1', 'PC', "ACHE", "IMPACT", "SHE", 'FOXA1', 'PPT2', 'ZFP36L1','AFF4', 'ANTXR2', "HDAC8", "VKORC1" )
  
  terms_of_interest <- c("B cells", "macrophages", "NK cells")
  pub <- pubscore(terms_of_interest = terms_of_interest, genes = selected_genes )
  pub <- test_score(pub, total_genes =  total_genes, remove_ambiguous = T, nsim = 1e4)
  p1 <- pub@p_value
  expect_that(pub@p_value, is_less_than(0.01))
  pub <- test_score(pub, total_genes =  total_genes, remove_ambiguous = T, ambiguous_terms = c('PC', "ACHE", "IMPACT"), nsim = 1e4)
  expect_that(pub@p_value, is_less_than(0.2))
  expect_that(pub@p_value, is_more_than(0.05))
})



test_that("test get literature  score when using gene2pubmed", {
  selected_genes <- c('CD79A', 'CD14', 'NKG7', 'CST3', 'AIF1')
    # These genes were selected from a panel containing the following genes
  total_genes <- c('CD79A', 'CD14', 'NKG7', 'CST3', 'AIF1', 'PC', "ACHE", "IMPACT", "SHE", 'FOXA1', 'PPT2', 'ZFP36L1','AFF4', 'ANTXR2', "HDAC8", "VKORC1" )
  
  terms_of_interest <- c("B cells", "macrophages", "NK cells")
  pub <- pubscore(terms_of_interest = terms_of_interest, genes = selected_genes,gene2pubmed = TRUE )
  pub <- test_score(pub, total_genes =  total_genes, remove_ambiguous = T, nsim = 1e4)
  p1 <- pub@p_value
  expect_that(pub@p_value, is_less_than(0.01))
  pub <- test_score(pub, total_genes =  total_genes, remove_ambiguous = T, ambiguous_terms = c("SHE"), nsim = 1e4)
  expect_that(pub@p_value, is_less_than(0.01))
})




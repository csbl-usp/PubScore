test_that("get literature works", {
  gene <- 'CD4'
  terms_of_interest <-
    c(
      "CD4 T cell",
      "CD14+ Monocyte",
      "B cell",
      "CD8 T cell",
      "FCGR3A+ Monocyte",
      "NK cell",
      "Dendritic cell",
      "Megakaryocyte",
      'immunity'
    )
  obj <- get_literature_score(gene, terms_of_interest, max_score = 500)
  expect_that(obj$max_score, equals(500))
  
})

# simulation of null distribution for the l-score
# tiabgo lubiana 03/10/2018

library(iscore)

load("data/pbmc_tutorial.rds")
exprs <- as.data.frame(pbmc@scale.data)
terms_of_interest <- c("CD4 T cell", "CD14+ Monocyte", "B cell", "CD8 T cell", 
                       "FCGR3A+ Monocyte", "NK cell", "Dendritic cell", "Megakaryocyte", 'immunity')


rm(list = ls()[!(ls() %in% c('exprs'))])
set.seed(2018)
all_genes <- sample(rownames(exprs))

simulation_of_literature_null_genes <- data.frame(2,2)
simulation_of_literature_null_genes <- simulation_of_literature_null_genes[-1,]

#getting scores for all genes it can before tomorrow morning

for (i in all_genes){
  tryCatch({
   bla <- get_iscore(i, terms_of_interest)
   new_line <- data.frame(i, bla$gene_iscore)
   simulation_of_literature_null_genes <- rbind(simulation_of_literature_null_genes, new_line)
   Sys.sleep(0.2)}, error = function(e){print(e)
  })
   
}

distribution_of_scores <- c()
for (i in 1:100000){
  list_score <- mean(sample(simulation_of_literature_null_genes$bla.gene_iscore, 80))
  distribution_of_scores <- c(distribution_of_scores,list_score )
  
}

distribution_of_scores<- as.data.frame(distribution_of_scores)

save(simulation_of_literature_null_genes, file = './int/null_score_genes.RData')
MAX_SCORE = max(distribution_of_scores)
MIN_SCORE = min(distribution_of_scores)

library(ggplot2)
ggplot(distribution_of_scores, aes(x = distribution_of_scores))+
  geom_density() +
  ggtitle('Density of the distribution of null scores(10^6 resamplings)') +
  xlab('Average number of counts on PubMed') +
  ylab('Proportion of gene-lists for each category')

ggsave(filename = './results/density_null_scores.pdf')


ggplot(distribution_of_scores, aes(x = distribution_of_scores))+
  geom_boxplot(aes(y=distribution_of_scores)) +
  ggtitle('Boxplot of the distribution of null scores(10^6 resamplings)') +
  xlab('80-gene list') + 
  ylab('Average number of counts on PubMed')+
  annotate('text', x = 150, y = 150, label = paste('Max score:', MAX_SCORE))+
  annotate('text', x = 150, y = 120, label = paste('Min score:', MIN_SCORE))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename = './results/boxplot_null_scores.pdf')


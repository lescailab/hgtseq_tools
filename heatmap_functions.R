library(tidyverse)
library(furrr)

dataset = readRDS("/home/lescailab/COLLABS/malacrida/glossina_hgtseq_run02/classified_reads_collated_single.rds")
small_dataset = head(dataset, 20)

table1 = tibble(
  "read name" = small_dataset$read_name,
  "taxID" = small_dataset$taxID,
  "k-mers info" = small_dataset$`k-mers info`
)

###########
## ideal ##
###########

nodes_rds_object = readRDS("/home/lescailab/local_REFS/taxonomy_sql/nodes_parsed.rds")

parseTaxAssignment <- function(kmersinfo, current_taxid, nodes=nodes_rds_object){
  tax_tree = calcTree(current_taxid, nodes)
  tax_score = evalKrakenQual(kmersinfo, tax_tree)
  tax_heatmap = heatmapList(kmersinfo, tax_tree)
  results = list(tax_tree, tax_score, tax_heatmap)
  return(results)
}

heatmap_data = table1 %>% 
  mutate(
    heatmap_var = map2(`k-mers info`, `taxID`, parseTaxAssignment)
  ) 

heatmap_table = heatmap_data %>% 
  select(heatmap_var) %>% 
  unnest_wider(heatmap_var, names_sep = "_") %>% 
  as.matrix()

hmap = as.matrix(heatmap_table[,3])
row.names(hmap) <- heatmap_data$`read name`

x1 = lapply(hmap, "length<-", max(lengths(hmap)))
names(x1) = heatmap_data$`read name`

final_matrix = t(bind_rows(x1, .id = NULL ))
class(final_matrix) <- "numeric"
heatmap(final_matrix, Colv = NA, Rowv = NA, scale = "column")

#####################
### tree function ###
#####################

calcTree <- function(current_taxid, nodes){
  tree = c(as.character(current_taxid))
  parent = 0
  while (length(parent)>0 & parent!=1){
    parent = nodes %>% 
      filter(id == current_taxid) %>% 
      pull(parent)
    if(length(parent)>0 & parent != 1){
      tree = c(tree, as.character(parent))
      current_taxid = parent
    }
  }
  return(tree)
}

######################
### score function ###
######################

evalKrakenQual <- function(kmersinfo, tree){
  kmersTab <- as_tibble(
    matrix(ncol = 2, byrow = TRUE,
           data = unlist(strsplit(
             unlist(strsplit(kmersinfo, " ")), 
             ":"))),
    .name_repair = "minimal"
  )
  names(kmersTab) <- c("taxid", "kmers_found")
  kmersTab$kmers_found = as.numeric(kmersTab$kmers_found)
  sum_all_kmers = sum(kmersTab$kmers_found)
  kmersTab = kmersTab %>%
    filter(taxid != "0")
  
  sum_tax_classified_kmers= kmersTab %>% 
    filter(
      taxid %in% tree
    ) %>% 
    mutate(
      kmers_found = as.numeric(kmers_found),
    ) %>% 
    pull(kmers_found) %>% 
    sum()
  
  taxid_score = sum_tax_classified_kmers/sum_all_kmers
  
  return(taxid_score)
}

######################
## heatmap function ##
######################

heatmapList <- function(kmersinfo, tree){
  kmersTab <- as_tibble(
    matrix(ncol = 2, byrow = TRUE,
           data = unlist(strsplit(
             unlist(strsplit(kmersinfo, " ")), 
             ":"))),
    .name_repair = "minimal"
  )
  names(kmersTab) <- c("taxid", "kmers_found")
  kmersTab$kmers_found = as.numeric(kmersTab$kmers_found)
  lista = vector()
  for (l in 1:nrow(kmersTab)){
    row_example = kmersTab[l,]
    tax = row_example$taxid
    kmer = row_example$kmers_found
    for (i in 1:kmer){
      tax_tree = ifelse(tax %in% tree, 1, 0)
      lista = c(lista, tax_tree)
    }
  }
  return(lista)
}

#############################################################################
# multicore with furrr - tested as script (rstudio does not like multicore) #
#############################################################################
save.image("tests_mapping.RData")

performance_tests = tibble(
  lines = c(rep(2000,5), rep(5000,5), rep(2000,5), rep(5000,5)),
  cores = rep(c(2,4,8,16,32),4),
  method = c(rep("tidyverse", 10), rep("furrr", 10)),
  runtime =c(
    c(2.768, 2.749, 2.77, 2.797, 2.738), ## tidyverse with 2000 lines
    c(6.22, 6.242, 6.28, 6.24, 6.201), ## tidyverse with 5000 lines
    c(2.209, 1.657, 1.332, 1.267, 1.435), ## furrr with 2000 lines
    c(5.223, 4.246, 3.391, 3.047, 3.06) ## furrr with 5000 lines
  )
)

ggplot(performance_tests, 
       aes(x=cores, y=runtime, colour = method))+
  geom_point()+
  geom_smooth(method = "loess")+
  geom_vline(xintercept = 8, colour = "blue")+
  facet_wrap(lines~., scales = "free")



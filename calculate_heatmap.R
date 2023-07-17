library(tidyverse)
library(furrr)

### how to run this code ###
# Rscript calculate_heaptmap.R /path/to/classified_reads.rds ifparallel numbercores #

nodes_rds_object = readRDS("/home/lescailab/local_REFS/taxonomy_sql/nodes_parsed.rds")

args            = commandArgs(trailingOnly = TRUE)
dataset         = readRDS(args[1])
parallel_choice = args[2]
cores           = args[3]

#####################
### tree function ###
#####################

calcTree <- function(current_taxid, nodes){
  tree <- tryCatch(
    {
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
    },
    error = function(cond){
      empty_tree = c(NA)
      message("-----------------------------------")
      message(paste0("error: tax id ", current_taxid, " is not present in nodes, so the function returns an empy tree."))
      return(empty_tree)
    }
  )
  return(tree)
}


writeLines(" ")
writeLines("--- tree function saved ---")

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

writeLines("--- score function saved ---")

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

writeLines("--- heatmap function saved ---")
writeLines(" ")

#####################
## global function ##
#####################

parseTaxAssignment <- function(kmersinfo, current_taxid, nodes=nodes_rds_object){
  message(paste0("analysing tax id ", current_taxid))
  tax_tree = calcTree(current_taxid, nodes)
  tax_score = evalKrakenQual(kmersinfo, tax_tree)
  tax_heatmap = heatmapList(kmersinfo, tax_tree)
  results = list(tax_tree, tax_score, tax_heatmap)
  return(results)
}

writeLines("--- global function saved ---")

if (parallel_choice == "parallel"){


##########################
## run with future_map2 ##
##########################

    writeLines(" ")
    writeLines("-------------------------------")
    writeLines("--- running global function ---")
    writeLines("-------------------------------")
    writeLines(" ")

    plan(multicore, workers = cores)
    heatmap_data = dataset %>% 
      mutate(
        heatmap_var = future_map2(`k-mers info`, `taxID`, parseTaxAssignment)
      ) 

    writeLines(" ")
    writeLines("------------------------------------------")
    writeLines("--- global function executed correctly ---")
    writeLines("------------------------------------------")
    writeLines(" ")

} else {

    ###################
    ## run with map2 ##
    ###################

    writeLines(" ")
    writeLines("-------------------------------")
    writeLines("--- running global function ---")
    writeLines("-------------------------------")
    writeLines(" ")

    heatmap_data = dataset %>% 
      mutate(
        heatmap_var = map2(`k-mers info`, `taxID`, parseTaxAssignment)
      ) 

    writeLines(" ")
    writeLines("------------------------------------------")
    writeLines("--- global function executed correctly ---")
    writeLines("------------------------------------------")
    writeLines(" ")


}

######################
## heatmap assembly ##
######################

heatmap_table = heatmap_data %>% 
  select(heatmap_var) %>% 
  unnest_wider(heatmap_var, names_sep = "_") %>% 
  as.matrix()

heatmap_list = as.matrix(heatmap_table[,3])
row.names(heatmap_list) <- heatmap_data$`read_name`

matrix_fill = lapply(heatmap_list, "length<-", max(lengths(heatmap_list)))
names(matrix_fill) = heatmap_data$`read_name`

heatmap_matrix = t(bind_rows(matrix_fill, .id = NULL ))
class(heatmap_matrix) <- "numeric"

saveRDS(heatmap_matrix, file = "heatmap_matrix.rds")

# pdf("heatmap.pdf")
# heatmap(heatmap_matrix, Colv = NA, Rowv = NA, scale = "column")
# dev.off()


writeLines(" ")
writeLines("-----------------------------------")
writeLines("--- heatmap generated correctly ---")
writeLines("-----------------------------------")

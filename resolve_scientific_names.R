#!/usr/bin/env Rscript

writeLines("--- loading library")
library(taxonomizr)
library(tidyverse)

sqldb = "/home/lescailab/local_REFS/taxonomy_sql/accessionTaxa.sql"

getGenusFromSpecies <- function(taxaid){
  taxonomy = getTaxonomy(taxaid, sqldb)
  genus = as_tibble(taxonomy)$genus[[1]]
  return(genus)
}

getLowestRankName <- function(taxaid, result="name"){
  taxonomy = as_tibble(getTaxonomy(taxaid, sqldb))
  if(!is.na(taxonomy$species[[1]])){
    scientific_name = taxonomy$species[[1]]
    rank_level = "species"
  } else if (!is.na(taxonomy$genus[[1]])){
    scientific_name = taxonomy$genus[[1]]
    rank_level = "genus"
  } else if (!is.na(taxonomy$family[[1]])){
    scientific_name = taxonomy$family[[1]]
    rank_level = "family"
  } else if (!is.na(taxonomy$order[[1]])){
    scientific_name = taxonomy$order[[1]]
    rank_level = "order"
  } else if (!is.na(taxonomy$class[[1]])){
    scientific_name = taxonomy$class[[1]]
    rank_level = "class"
  } else if (!is.na(taxonomy$phylum[[1]])){
    scientific_name = taxonomy$phylum[[1]]
    rank_level = "phylum"
  } else if (!is.na(taxonomy$superkingdom[[1]])){
    scientific_name = taxonomy$superkingdom[[1]]
    rank_level = "superkingdom"
  } else {
    scientific_name = NA
    rank_level = NA
  }
  if(result=="name"){
    return(scientific_name)
  }
  else {
    return(rank_level)
  }
}

writeLines("--- reading the data")
collated_single = readRDS("/home/lescailab/COLLABS/malacrida/glossina_hgtseq_run02/classified_reads_collated_single.rds")
collated_both = readRDS("/home/lescailab/COLLABS/malacrida/glossina_hgtseq_run02/classified_reads_collated_both.rds")

writeLines("--- resolving single reads")
collated_single_resolved = tibble(
  taxID = unique(collated_single$taxID),
  sciname = unlist(lapply(unique(collated_single$taxID), getLowestRankName)),
  rank = unlist(lapply(unique(collated_single$taxID), getLowestRankName, "rank")),
  genus = unlist(lapply(unique(collated_single$taxID), getGenusFromSpecies))
)

writeLines("--- resolving both reads")
collated_both_resolved = tibble(
  taxID = unique(collated_both$taxID),
  sciname = unlist(lapply(unique(collated_both$taxID), getLowestRankName)),
  rank = unlist(lapply(unique(collated_both$taxID), getLowestRankName, "rank")),
  genus = unlist(lapply(unique(collated_both$taxID), getGenusFromSpecies))
)


#####################################################
## creating original read-based resolved datasets ###
#####################################################


writeLines("--- joining resolved and original data for single")
collated_single = collated_single %>%
  left_join(collated_single_resolved, by = "taxID")
writeLines("--- saving results")
saveRDS(collated_single, file = "collated_single_resolved.rds")

writeLines("--- joining resolved and original data for both")
collated_both = collated_both %>%
  left_join(collated_both_resolved, by = "taxID")
writeLines("--- saving results")
saveRDS(collated_both, file = "collated_both_resolved.rds")

writeLines("--- JOB COMPLETED")
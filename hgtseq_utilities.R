library(tidyverse)

evalKrakenQual <- function(kmersinfo){
  kmersTab <- as_tibble(
    matrix(ncol = 2, byrow = TRUE,
           data = unlist(strsplit(
             unlist(strsplit(kmersinfo, " ")), 
             ":"))),
    .name_repair = "minimal"
  )
  names(kmersTab) <- c("taxid", "kmers_found")
  kmersTab$kmers_found = as.numeric(kmersTab$kmers_found)
  kmersTab = kmersTab %>%
    group_by(taxid) %>%
    summarise(tot_kmers = sum(kmers_found))
  ## sum kmers for scoring includes kmers not classified, i.e. 0
  sum_all_kmers = sum(kmersTab$tot_kmers)
  ## then we can filter them out to calculate max == assigned taxid
  kmersTab = kmersTab %>%
    filter(taxid != "0")
  max_taxid = paste0(kmersTab$taxid[which(kmersTab$tot_kmers == max(kmersTab$tot_kmers))],
                     collapse = ",")
  max_kmers = max(kmersTab$tot_kmers)
  taxid_score = max(kmersTab$tot_kmers)/sum_all_kmers
  taxid_classified_score = max(kmersTab$tot_kmers)/sum(kmersTab$tot_kmers)
  parsed = tibble(
    max_taxid = max_taxid,
    max_kmers = max_kmers,
    sum_all_kmers = sum_all_kmers,
    sum_classified_kmers = sum(kmersTab$tot_kmers),
    taxid_score = taxid_score,
    taxid_classified_score = taxid_classified_score
  )
  return(parsed)
}


#### example run
# scored_data = cbind(joint_data,
#                     do.call(rbind,
#                             lapply(joint_data$`k-mers info`, evalKrakenQual)))




#### ANNOTATION FROM CHIPSEQ PACKAGE
# library(ChIPpeakAnno)
# library(org.Hs.eg.db)
# annoData <- readRDS(url("https://github.com/lescai-teaching/datasets_reference_only/raw/main/annotations/chipseq_annoData.RData"))
# integration_noextra <- annotatePeakInBatch(integration_noextra, 
#                                            AnnotationData=annoData, 
#                                            output="nearestBiDirectionalPromoters",
#                                            bindingRegion=c(-2000, 500))
# integration_noextra <- addGeneIDs(integration_noextra,
#                                   "org.Hs.eg.db",
#                                   IDs2Add = c("symbol"))


blacklist <- read_tsv(url("https://raw.githubusercontent.com/lescailab/hgt_blacklist/master/genera_blacklist.txt"), col_names = "genus") %>% pull(genus)
human_flags <- read_tsv(url("https://raw.githubusercontent.com/lescailab/hgt_blacklist/master/human_related_genera_blacklist.txt"), col_names = "genus") %>% pull(genus)

blacklistdata <- tibble(
  genus = blacklist,
  flag = ifelse(blacklist %in% human_flags, "flagged", "blacklisted")
)

getBlacklistStatus <- function(genus_name){
  if (genus_name %in% blacklistdata$genus){
    flag = blacklistdata %>%
      filter(genus == genus_name) %>%
      pull(flag)
  } else {
    flag = "not_blacklisted"
  }
  return(flag)
}



########################################
### CREATE GRANGES WITH CORRECT SEQINFO
### FROM HGTSEQ INTEGRATION TIBBLE
########################################


createIntegrationGR <- function(data,seqinfo){
  
  ## load package if not present
  require(GenomicRanges)
  require(tidyverse)
  
  ## create granges object
  integration <- GRanges(
    data$chr,
    IRanges(
      start = data$position,
      end = data$position
    )
  )
  mcols(integration)<-data
  
  ## remove extra chromosomes
  
  integration_noextra <- integration[
    seqnames(integration) %in% seqnames(seqinfo)
  ]
  
  ## match levels of GRanges object with genome
  matchinglevels = seqlevels(seqinfo)[seqlevels(seqinfo) %in% seqlevels(integration_noextra)]
  integration_noextra = keepSeqlevels(integration_noextra, matchinglevels)
  ## match seqlengths as well
  seqlengths(integration_noextra) = seqlengths(seqinfo)[seqlevels(seqinfo) %in% seqlevels(integration_noextra)]
  
  
  
  return(integration_noextra)
  
}


###############################################
## create genome from sequence dictionary
###############################################

createGenomeGR <- function(dictionary, has_chrX=TRUE, has_chrY=TRUE, autosomal_chrs=c(1:22), genome_name=NULL, chr_prefix=NULL){
  
  require(GenomicRanges)
  require(tidyverse)
  
  dictred = read_tsv(dictionary, skip = 1, col_names = c("class", "sequence", "length", "m5", "as", "url", "species"))
  
  dictred$sequence <- gsub("SN:", "", dictred$sequence)
  dictred$length <- gsub("LN:", "", dictred$length)
  
  newSeqInfoData <- Seqinfo(
    seqnames = dictred$sequence,
    seqlengths = as.numeric(dictred$length),
    isCircular = NA,
    genome = genome_name
  )
  
  standard_chromosomes = paste0(chr_prefix, as.character(autosomal_chrs))
  if(has_chrX){
    standard_chromosomes = c(standard_chromosomes, paste0(chr_prefix, "X"))
  }
  if(has_chrY){
    standard_chromosomes = c(standard_chromosomes, paste0(chr_prefix, "Y"))
  }
  
  newSeqInfoData_noextra <- newSeqInfoData[standard_chromosomes]
  
  genomeGR <- GRanges(
    seqnames(newSeqInfoData_noextra),
    IRanges(
      start = 1,
      end = unname(seqlengths(newSeqInfoData_noextra))
    )
  )
  seqinfo(genomeGR) <- newSeqInfoData_noextra
  
  return(genomeGR)
}









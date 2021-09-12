# Convert R data into Python inputs
require(reshape2)
library(dplyr)

## Set up parameters
inoutpath <- "/home/xiaopengxu/Desktop/data-covid-review/2021-05-31/"

subpath <- function(fname){
  return(file.path(inoutpath,fname))
}

## Retrieve compdata_ext from webpage_data, save to csv
load(subpath("webpage/webpage_data.Rdata"))
write.csv(compdata_ext,subpath("compdata_ext.csv"))

## Get references betweeen papers, save to csv
s_scholar_all <- readRDS(subpath("s_scholar_compdata.RDS"))
dup_index <- readRDS(subpath("dup_index.RDS"))

s_scholar <- s_scholar_all[-dup_index]
compdata_ext_ref <- compdata_ext

ref <- sapply(s_scholar, function(a)if(!any(is.na(a))) paste(a$references$doi, collapse=', ' ) else '')
p_doi <- sapply(s_scholar, function(a)if(!any(is.na(a))) paste(a$doi, collapse=', ' )  else '')

for(i in 1:length(compdata_ext_ref$title)) {
  compdata_ext_ref$ref[i] = ref[i]
  compdata_ext_ref$p_doi[i] = p_doi[i]
}

write.csv(compdata_ext_ref, subpath("compdata_ext_ref.csv"))

## Get raw data from PubMed and arxiv for label verification
raw_pubmed <- readRDS(subpath("cache/pubmed.RDS"))
raw_arxiv <- readRDS(subpath("cache/arxiv.RDS"))
raw_biomedRxiv <- readRDS(subpath("cache/biomedRxiv.RDS"))
#raw_biomedRxiv_details <- readRDS(subpath("cache/biomedRxiv_details.RDS"))

### biomedRxiv convertion
data <- raw_biomedRxiv$rels
data$id <- rownames(data) 
data$rel_date = as.Date(data$rel_date)
data$rel_authors <- vapply(data$rel_authors, paste, collapse = ", ", character(1L))

biomedRxiv <- data.frame(data)
write.csv(biomedRxiv, subpath("biomedRxiv.csv"))

### aRxiv convertion
write.csv(raw_arxiv, subpath("aRxiv.csv"))

### Pubmed convertion
getPub <- function(pubmed, member)
  do.call(c, lapply(pubmed, member))

pubmed <- data.frame(
  title = getPub(raw_pubmed, ArticleTitle),
  abstract = getPub(raw_pubmed, AbstractText),
  journal = getPub(raw_pubmed, ISOAbbreviation),
  doi = getPub(raw_pubmed, DOI),
  pubYear = getPub(raw_pubmed, YearPubDate),
  pubMonth = getPub(raw_pubmed, MonthPubDate),
  pubDay = getPub(raw_pubmed, DayPubDate),
  acptYear = getPub(raw_pubmed, YearAccepted),
  acptMonth = getPub(raw_pubmed, MonthAccepted),
  acptDay = getPub(raw_pubmed, DayAccepted),
  stringsAsFactors = F
)

pubmed_ok <- pubmed %>% mutate(abstract = as.character(abstract))
write.csv(pubmed_ok, subpath("PubMed.csv"))

## Author and Article metrics convertion
CSCoV_scores <- readRDS(subpath("CSCoV.scores.RDS"))
write.csv(CSCoV_scores, subpath("CSCoV_scores.csv"))
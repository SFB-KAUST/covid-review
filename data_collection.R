
## OPTIONS

# path where downlaoded data are stored
inoutpath <- "datanew"
logfile <- "data_collection.log"
# path to local file containing CORD-19 database
# (see https://www.semanticscholar.org/cord19)
cord_path <- "./cord-19_2021-08-19.tar.gz"

# when TRUE, data are downloaded. When FALSE, they are loaded from previous
# cached downloads. Since RISmed library has recently broken, Pubmed papers
# are aloways loaded from a local CORD-19 file. However they are only
# preprocessed if fetchPubmed is TRUE.
fetchPubmed <- F
fetchBioMed <- F
fetchArxiv <- F
fetchDownloads <- F
fetchScholarArticles <- F
fetchScholarAuthors <- F
errToFile <- F # redirect stderror to .err file


## CHECK OPTIONS

subpath <- function(fname){
  return(file.path(inoutpath,fname))
}

checkpath <- function(p){
  if(!file.exists(p))
    stop(paste("Path does not exist:", p))
}

checkpath(inoutpath)

if(!file.exists(subpath("cache")))
  dir.create(subpath("cache"))
if(!file.exists(subpath("figs")))
  dir.create(subpath("figs"))

checkpath(subpath("cache/author_metrics"))


## START

# function to print messages and log them
log <- function(txt){
  txt <- paste0("[", Sys.time(), "] ", txt)
  write(txt, file=subpath(logfile), append=TRUE)
  print(txt)
}

if(errToFile)
  sink(file(subpath(paste0(logfile,".err")), open = "a"),
       type = "message", append = T) else sink(NULL, type="message")

log("Initializing...")

library(dplyr)
library(ggplot2)
library(RISmed)
library(reshape2) # dcast
library(ggpubr) # ggarrange
library(ggstatsplot) # combine_plots


getPub <- function(pubmed, member)
  do.call(c, lapply(pubmed, member))

# this is to put all sources in the same format
preprocess_data <- function(arxiv, biorxiv, pubmed) {
  if(!is.null(pubmed)) {
    pubmed_ok <- data.frame(
      title = getPub(pubmed, ArticleTitle),
      abstract = getPub(pubmed, AbstractText),
      journal = getPub(pubmed, ISOAbbreviation),
      DOI = getPub(pubmed, PMID),
      date = as.Date(paste(
        getPub(pubmed, YearPubmed),
        getPub(pubmed, MonthPubmed),
        getPub(pubmed, DayPubmed), sep="-")),
      collection = "pubmed",
      published = NA,
      nauthors = sapply(getPub(pubmed, Author), nrow),
      stringsAsFactors = F
    )
  } else { # using CORD as RISmed library doesn't work anymore
    if(fetchPubmed) { # preprocessing CORD-19 data
      preprocess_cord <- function() # just to create a temporary scope
      {
        log(paste("Unzipping CORD-19:", cord_path))
        cord_version <- gsub("cord-19_|.tar.gz", "", basename(cord_path))
        untar(cord_path, exdir = file.path(dirname(cord_path), "cord"),
              files = file.path(cord_version, "metadata.csv"))
        log(paste("Loading CORD:", file.path(cord_version, "metadata.csv")))
        cord19 = read.csv(file.path(dirname(cord_path), "cord", cord_version, "metadata.csv"))
        intitle <- grep("covid-19|sars-cov-2", tolower(cord19$title))
        inabstract <- grep("covid-19|sars-cov-2", tolower(cord19$abstract))
        cord19_filtered <- cord19[union(intitle, inabstract), ]
        cord19_new <- cord19_filtered[, c("title", "abstract", "journal", "doi", "publish_time", "source_x")]
        colnames(cord19_new) <- c("title", "abstract", "journal", "DOI", "date", "collection")
        cord19_new$nauthors = sapply(strsplit(as.character(cord19_filtered$authors), ";"), length)
        cord19_new$abstract[cord19_new$abstract==""] <- cord19_new$title[cord19_new$abstract==""]
        cord19_new$abstract <- gsub("&lt|&quot", "", cord19_new$abstract)
        log(paste("Loaded", nrow(cord19_new), "entries"))
        PMids <- union(grep("Medline", cord19_new$collection), grep("PMC", cord19_new$collection))
        pubmed_ok <- cord19_new[PMids,]
        pubmed_ok$collection <- "pubmed"
        log(paste("Extracted", nrow(pubmed_ok), "entries"))
        pubmed_ok$published <- NA
        pubmed_ok <- pubmed_ok[, c(1:6,8,7)] # reordering to stay safe
        pubmed_ok$abstract[is.na(pubmed_ok$abstract)] <- ""
        pubmed_ok$title[is.na(pubmed_ok$title)] <- ""
        pubmed_ok$date <- as.Date(pubmed_ok$date)
        log("Storing cord_pubmed")
        saveRDS(pubmed_ok, subpath("cache/cord_pubmed.RDS"))
        return(pubmed_ok)
      }
      pubmed_ok <- preprocess_cord()
    } else pubmed_ok <- readRDS(subpath("cache/cord_pubmed.RDS"))
  }
  
  pubmed_ok <- pubmed_ok %>% mutate(abstract = as.character(abstract))
  
  arxiv$journal_ref[arxiv$journal_ref == ""] <- NA
  arxiv_ok <- data.frame(
    title = arxiv$title,
    abstract = arxiv$abstract,
    journal = arxiv$journal_ref,
    DOI = substr(arxiv$id, 1, 10),
    date = as.Date(arxiv$submitted),
    collection = "arxiv",
    published = arxiv$journal_ref,
    nauthors = x <- sapply(arxiv$authors, function(x) length(strsplit(x,"|",fixed=T)[[1]])),
    stringsAsFactors = F
  )
  
  biorxiv <- biorxiv$rel
  biorxiv_ok <- data.frame(
    title = biorxiv$rel_title,
    abstract = biorxiv$rel_abs,
    journal = NA,
    DOI = biorxiv$rel_doi,
    date = as.Date(biorxiv$rel_date),
    collection = biorxiv$rel_site,
    published = biorxiv$published,
    nauthors = sapply(biorxiv$rel_authors, function(x) length(x$author_name)),
    stringsAsFactors = F
  )
  
  alldata <- rbind(pubmed_ok, arxiv_ok, biorxiv_ok)
  alldata$abstract[alldata$abstract==""] <- alldata$title[alldata$abstract==""]
  alldata$abstract <- gsub("&lt|&quot", "", alldata$abstract)
  return(alldata)
}


# Querying PubMed

if(F) { # using CORD-19 data instead of querying pubmed
  
query <- '((covid-19[Title/Abstract]) OR (SARS-COV-2[Title/Abstract]))'
#the following PMIDs return errors and must be skipped:
pmid_exclude <- c("33690889")

if(fetchPubmed) {
  log("Querying Pubmed...")
  npapers <- EUtilsSummary(query, retmax=1)@count
  log(paste("Items to fetch:", npapers))
  chunksize <- 5000
  pubmed <- list()
    #for(i in 1:ceiling(npapers/chunksize)){
  for(i in 1:ceiling(npapers/chunksize)){
      log(paste0("Querying chunk ", i, " of ", ceiling(npapers/chunksize), "..."))
      goterr <- T
      while(goterr) {
        attempts <- 1
        tryCatch({
          if(i==21) ## weird bug
            qi <- EUtilsSummary(query, retstart=(i-1)*chunksize+1, retmax=chunksize) else
              qi <- EUtilsSummary(query, retstart=(i-1)*chunksize+0, retmax=chunksize)
            goterr <- F
        }, error=function(e) {
          if(attempts > 15)
            stop("Quitting after 15 failed attempts.")
          log("Got an error, retrying...")
          goterr <- T
          attempts <- attempts + 1
        }
        )
        qi@PMID <- setdiff(qi@PMID,pmid_exclude)
        pubmed[[i]] <- EUtilsGet(qi)
      }
    }

  saveRDS(pubmed, subpath("cache/pubmed.RDS"))
  when <- date()
} else {
  log("Querying Pubmed (dry run)...")    
  pubmed <- readRDS(subpath("cache/pubmed.RDS"))
  when <- file.mtime(subpath("cache/pubmed.RDS"))
}

  log(paste("Total entries:",
            sum(sapply(pubmed, function(x) length(x@ArticleId)))))
} else {
  pubmed <- NULL
}




## Querying bioRxiv and medRxiv

if(fetchBioMed) {
  log("Querying bioRxiv and medRxiv...")
  options(timeout= 4000000) 
  library(RCurl)
  library(jsonlite)
  
  ### API documentation is at https://api.biorxiv.org/covid19/help ###
  
  group_url <- "https://api.biorxiv.org/covid19/PAGING"
  biorxivPAGE <- jsonlite::fromJSON(gsub("PAGING","0", group_url))
  biorxivTOT <- biorxivPAGE$messages$total
  
  biorxivL <- list()
  i <- 0
  while(i*30 < biorxivTOT){
    cat(i*30, "of", biorxivTOT, "\r")
    biorxivPAGE <- jsonlite::fromJSON(gsub("PAGING", i*30, group_url))
    biorxivL[[length(biorxivL)+1]] <- biorxivPAGE$collection
    i <- i+1
  }
  
  biorxiv <- do.call(rbind, biorxivL)
  
  #for back compatibility:
  temp <- biorxiv
  biorxiv <- list()
  biorxiv$rels <- temp
  rm(temp)
  biorxiv$rels$rel_site <- tolower(biorxiv$rels$rel_site)
  saveRDS(biorxiv, subpath("cache/biomedRxiv.RDS"))
  
  log("Downloading article details...")
  
  details <- list()
  for(i in 1:nrow(biorxiv$rels))
  {
    cat(i, "of", nrow(biorxiv$rels), "\r")
    details[[i]] <- tryCatch(
      jsonlite::fromJSON(
        paste0("https://api.biorxiv.org/details/",
               biorxiv$rels$rel_site[i], "/",
               biorxiv$rels$rel_doi[i])
      ),
      error = function(e){
        log(paste0("Error at entry: ", i))
        return(NA)
      }
    )
  }
  when <- date()
  saveRDS(details, subpath("cache/biomedRxiv_details.RDS"))
} else {
    log("Querying bioRxiv and medRxiv (dry run)...")    
    biorxiv <- readRDS(subpath("cache/biomedRxiv.RDS"))
    details <- readRDS(subpath("cache/biomedRxiv_details.RDS"))
    when <- file.mtime(subpath("cache/biomedRxiv.RDS"))
}

log(paste("Total entries:", nrow(biorxiv$rels)))  



published <- sapply(details, function(x) {
  if("collection" %in% names(x)){
    if(!is.null(x$collection$published[1]))
      return(x$collection$published[1])
  }
  return(NA)
}
)
biorxiv$rels <- cbind(biorxiv$rels, published = published)





## Querying arXiv

library(aRxiv)
q <- 'abs:sars OR ti:sars OR abs:covid OR ti:covid'

if(fetchArxiv) {
  log("Querying arXiv...")
  qc <- arxiv_count(q)
  log(paste("Will fetch this number of items:", qc))
  arxiv <- NULL
  attempts <- 1
  repeat { #sometimes it quits at partial downloads
    arxiv <- rbind(arxiv, arxiv_search(q, limit = qc, start=nrow(arxiv)))
    if(nrow(arxiv) >= qc)
      break
    attempts <- attempts + 1
    if(attempts > 20)
      stop("arXiv download didn't finish after 20 attemtps")
    log(paste0("Stopped earlier by server at size: ",
               nrow(arxiv), ". Querying more..."))
  }
  saveRDS(arxiv, subpath("cache/arxiv.RDS"))
  when <- date()
} else {
  log("Querying arXiv (dry run)...")
  arxiv <- readRDS(subpath("cache/arxiv.RDS"))
  when <- file.mtime(subpath("cache/arxiv.RDS"))
}

log(paste("Total entries:", nrow(arxiv)))

matchstr <- paste(arxiv$title, arxiv$abstract)
w <- grep("covid-19|sars-cov-2", tolower(matchstr))
arxiv <- arxiv[w,]
arxiv$submitted <- as.Date(arxiv$submitted)




## Joining data
log("Joining data...")
alldata <- preprocess_data(arxiv, biorxiv, pubmed)
saveRDS(alldata, subpath("alldata.RDS"))


tab <- table(alldata$collection)
tab <- c(tab, total=sum(tab))
log(paste("Joined entries:",
          paste(names(tab), tab, collapse=", ")))

showstats <- function(data) {
    log(paste("NA titles:", sum(is.na(data$title))))
    log(paste("NA abstracts:", sum(is.na(data$abstract))))
        log(paste("Duplicated (non-NA) titles:",
                  sum(duplicated(data$title, incomparables=NA))))
        log(paste("Duplicated (non-NA) abstracts:",
                  sum(duplicated(data$abstract, incomparables=NA))))
}
showstats(alldata)


log("Identifying computational articles...")
## Extraction of computational papers

## Keywords to identify computational contents
keywords <- c(
  "bayes",
  "model fit",
  "mathematical prediction",
  "simulation model",
  "model simulation",
  "virtual screening",
  "molecular dynamics simulation",
  "simulation experiment",
  "machine learning",
  "information mining",
  "supervised learning",
  "unsupervised learning",
  "deep learning",
  "computational",
  "artificial intelligence",
  "bioinformatic",
  "neural network",
  "mathematical model",
  "sequencing",
  "image analysis",
  "image processing",
  "data mining",
  "classifier",
  "in silico",
  "in-silico",
  "exom",
  "transcriptom",
  "proteom",
  "metagenom",
  " omics",
  "immunopeptidom",
  "omicscience",
  "metabolom",
  "pharmacogenom",
  "nutrigenom",
  "vaccinom",
  "phylogenom",
  "radiom",
  "interactom",
  "genome-wide",
  "whole-genome",
  "whole genome",
  "genomic epidemiology",
  "genome haplotypes",
  "genome sequences",
  "genomes",
  "genomic structure",
  "lipidom",
  "multiom",
  "immunome",
  "phenomics",
  "glycomics",
  "microbiom",
  "toxicogenom",
  "phosphoproteom",
  "glycoproteom",
  "metatranscriptom",
  "virom",
  "active learning",
  "adversarial network",
  "predictive model",
  "bayesian model",
  "deep model",
  "information retrieval",
  "multi-omics",
  "sequence alignment",
  "deep sequencing",
  "network analysis",
  "deep-learning",
  "knowledge graph",
  "network model",
  "structural model",
  "literature mining",
  "logistic regression",
  "text mining",
  "text-mining",
  "machine intelligence",
  "optimization model",
  "correlation analysis",
  "ensemble learning",
  "transfer learning",
  "prediction model",
  "forecasting model",
  "digital health",
  "lstm",
  "probabilistic",
  "dataset",
  "aurcroc"
)


## matching keywords in paper titles or abstracts
title_abs <- tolower(paste(alldata$title, alldata$abstract))
matched_papers <- grep(paste(keywords, collapse="|"), tolower(title_abs), fixed=F)
compdata <- alldata[matched_papers, ]

#removing papers with wrong publication date (NA or future dates)
compdata <- compdata[-which(compdata$date > as.Date(when) | is.na(compdata$date)),]

#computing reverse matches (which keywords matched each paper)
rev_matches <- sapply(
  title_abs[matched_papers],
  function(paper)
    paste(keywords[sapply(keywords, function(x) grepl(x, paper))], collapse=", "),
  USE.NAMES = F
)

compdata <- cbind(compdata, keywords=rev_matches)
# storing the database of computational studies, including duplicates
saveRDS(compdata, subpath("compdata_dups.RDS"))

tab <- table(compdata$collection)
tab <- c(tab, total=sum(tab))
log(paste("Computational entries:",
          paste(names(tab), tab, collapse=", ")))

showstats(compdata)


### Scraping download statistics from bioRxiv and medRxiv

idx <- compdata$collection=="biorxiv"|compdata$collection=="medrxiv"
dois <- compdata[idx,c("collection","DOI")]

if(fetchDownloads) {
  if(!file.exists(subpath("cache/biomedrxiv_scrape")))
    dir.create(subpath("cache/biomedrxiv_scrape"))
  
    log("Scraping download statistics...")
    cmd <- paste0('wget https://www.SERVER_HERE.org/content/DOI_HEREv1.article-metrics',
                ' -O FNAME_HERE')
  for(i in 1:nrow(dois)) {
    print(paste(i, "of", nrow(dois)))
    fname <- paste0(subpath("cache/biomedrxiv_scrape/"), gsub("[./]", "_", dois$DOI[i]))
    cmd_i <- gsub("SERVER_HERE", dois$collection[i], cmd)
    cmd_i <- gsub("DOI_HERE", dois$DOI[i], cmd_i)
    cmd_i <- gsub("FNAME_HERE", fname, cmd_i)
    system(cmd_i, ignore.stdout = T)
    Sys.sleep(30)
  }
  
  library(rvest)
  biomedstats <- list()
  # timed_biomedstats <- list()
  for(i in 1:length(dois$DOI)) {
    cat(i, "of", nrow(dois), "\r")
    fname <- paste0(subpath("cache/biomedrxiv_scrape/"), gsub("[./]", "_", dois$DOI[i]))
    
    stats_table <- tryCatch({
      html_table(html_nodes(read_html(file(fname)),
                            ".highwire-stats.sticky-enabled"))[[1]]
    }, error = function(e){
      log(paste("error:", i, fname))
      return(matrix(rep(NA,4),1))
    })
    biomedstats[[i]] <- apply(stats_table[,-1,drop=F],2,sum)
    # if(nrow(stats_table)>0){
    #  timed_biomedstats[[i]] <- cbind(stats_table,
    #                                  at_days=paper_age_at_date(compdata,i,stats_table))
    # }
  }
  
  biomedstats <- do.call(rbind, biomedstats)
  rownames(biomedstats) <- dois$DOI
  saveRDS(biomedstats, subpath("biomedrxiv_stats.RDS"))
  #saveRDS(timed_biomedstats, subpath("cache/timed_biomedrxiv_stats.RDS"))
  
} else {
    log("Scraping download statistics (dry run)...")
    biomedstats <- readRDS(subpath("biomedrxiv_stats.RDS"))
}
foundstats <- sum(apply(biomedstats,1,function(x)!all(is.na(x))))
log(paste("Download statistic entries:", foundstats))


## Scraping metrics from Semantic Scholar {#semantic_scrape}

if(fetchScholarArticles) {
  log("Getting article metrics from Semantic Scholar...")
  library(jsonlite)
  res <- list()
  jurl <- "https://api.semanticscholar.org/v1/paper/"
  
  s_scholar <- list()
  for(i in 1:nrow(compdata))
  {
    cat(i, "of", nrow(compdata), "\r")
    id <- compdata$DOI[i]
    if(compdata$collection[i]=="arxiv") {
      id <- paste0("arXiv:", id)
    } else if(compdata$collection[i]=="pubmed") {
      id <- paste0("PMID:", id)
    }
    
    s_scholar[[i]] <- tryCatch(
      jsonlite::fromJSON(paste0(jurl, id)),
      error = function(e){
        log(paste0("error: ", i, ", ", jurl, id))
        return(NA)
      }
    )
    Sys.sleep(5)
  } 
  saveRDS(s_scholar, subpath("s_scholar_compdata.RDS"))
  log(paste("Stored queries results:", nrow(s_scholar)))
} else {
    log("Getting article metrics from Semantic Scholar (dry run)...")
    s_scholar <- readRDS(subpath("s_scholar_compdata.RDS"))
    log(paste("Loaded queries results:", length(s_scholar)))
}





## getting Author metrics from Semantic Scholar

# note: previously downloaded author metrics are skipped, only new authors are
# downloaded. This takes days.

totauth <- sum(unlist(sapply(s_scholar, function(x) if(!is.na(x[1])) nrow(x$authors) else 0)))

if(fetchScholarAuthors) {
  log("Updating author metrics from Semantic Scholar...")
  library(jsonlite)
  jurl <- "https://api.semanticscholar.org/v1/author/"
  
  for(i in 1:length(s_scholar))
  {
    cat(i, "of", length(s_scholar), "\r")
    
    if(is.na(s_scholar[[i]])) next
    
    fname <- paste0(subpath("cache/author_metrics/"),
                    s_scholar[[i]]$paperId, ".RDS")
    
    if(file.exists(fname) || length(s_scholar[[i]]$authors)==0) next
    
    author_metrics_i <- data.frame(
      id = s_scholar[[i]]$authors$authorId,
      npapers = NA, infCit = NA
    )
    
    for(j in 1:nrow(author_metrics_i)) {
      a_data <- tryCatch(
        jsonlite::fromJSON(paste0(jurl, author_metrics_i$id[j])),
        error = function(e){
          log(paste("error:", i, j, jurl))
          return(list(papers=matrix(NA,0,0),
                      influentialCitationCount=NA))
        }
      )
      author_metrics_i[j, c("npapers","infCit")] <-
        c(nrow(a_data$papers), a_data$influentialCitationCount)
      Sys.sleep(3)
    }
    saveRDS(author_metrics_i, fname)
  }
} else {
    log("Updating author metrics from Semantic Scholar (dry run)...")
    }

log("Summarizing author metrics...")

authcount <- 0
fs <- list.files(subpath("cache/author_metrics"), "*.RDS")
authmet <- list()
for(i in 1:length(s_scholar)) {
  authmet[[i]] <- NA
  if(is.na(s_scholar[[i]][1])){next}
  id <- s_scholar[[i]]$paperId
  fname <- paste0(subpath("cache/author_metrics/"), id, ".RDS")
  if(file.exists(fname)){
    data <- readRDS(fname)
    authmet[[i]] <- c(sum=apply(data[,2:3],2,sum),
                      mean=apply(data[,2:3],2,mean),
                      max=apply(data[,2:3],2,max))
    authcount <- authcount + nrow(data)
  }
}
authmet <- do.call(rbind, authmet)
saveRDS(authmet, subpath("authmet.RDS"))
log(paste("Loaded author query results:", authcount))


## Handling duplicates (done afterwards for retro-compatibility issues)
log("Removing duplicates...")
abs_nopunc <- tolower(trimws(gsub('[[:punct:] ]+',' ',compdata$abstract)))
title_nopunc <- tolower(trimws(gsub('[[:punct:] ]+',' ',compdata$title)))
dups <- union(which(duplicated(abs_nopunc) & !is.na(abs_nopunc)),
              which(duplicated(title_nopunc) & !is.na(title_nopunc)))

# removing duplicates: if there's a preprint and a journal, prefer journal
preprints <- c("biorxiv", "medrxiv", "arxiv")
removeme <- vector("numeric")
for(i in seq_along(dups)){
  w <- which(compdata$abstract==compdata$abstract[dups[i]])
  if(length(w)<2) w <- which(compdata$title==compdata$title[dups[i]])
  # detecting preprints, may be both are, in that case will remove the first (newer)
  out <- vector("numeric")
  out <- which(tolower(c(compdata$journal[w[1]],
                         compdata$journal[w[2]])) %in% preprints)
  if(length(out)==0)
    out <- which(tolower(c(compdata$collection[w[1]],
                           compdata$collection[w[2]])) %in% preprints)
  if(length(out)==0)
    out <- 1
  
  if(length(out)>0)
    removeme <- c(removeme, w[out[1]])
}

## Exporting article scores
if(F){
log("Exporting scores...")
normbydate <- function(x, dates=compdata$date)
  #return(x / (as.numeric(as.Date(when)-dates)+1))
  return(x / (as.numeric(max(compdata$date)-dates)+1))


auth.papers <- authmet[,"mean.npapers"]
auth.papers.score <- ecdf(auth.papers)(auth.papers)
auth.citations <- authmet[,"mean.infCit"]
auth.citations.score <- ecdf(auth.citations)(auth.citations)
auth.score <- apply(cbind(auth.papers.score, auth.citations.score), 1,
                    mean, na.rm=T)
art.citations <- sapply(s_scholar, function(x) {
  if("citations" %in% names(x))
    if(!is.null(nrow(x[["citations"]])))
      return(nrow(x[["citations"]]))
  return(0)
})
art.citations.score <- normbydate(art.citations)
art.citations.score <- ecdf(art.citations)(art.citations)
biomed_idx <- compdata$collection=="biorxiv"|compdata$collection=="medrxiv"
# make sure the idx points to articles for which we have DL stats
stopifnot(compdata$DOI[biomed_idx]==rownames(biomedstats))
art.views <- art.views.score <- rep(NA, nrow(compdata))
art.views[biomed_idx] <- biomedstats[,"Abstract"]
art.views.score[biomed_idx] <- normbydate(biomedstats[,"Abstract"],
                                              compdata$date[biomed_idx])
art.views.score[biomed_idx] <-
  ecdf(art.views.score[biomed_idx])(art.views.score[biomed_idx])

allscores <- data.frame(
  auth.papers,
  auth.papers.score,
  auth.citations,
  auth.citations.score,
  art.citations,
  art.citations.score,
  art.views,
  art.views.score
)

allscores <- cbind(
  allscores, final.score=apply(
    allscores[,c("auth.papers.score", "auth.citations.score",
                 "art.citations.score", "art.views.score")],
    1, mean, na.rm=T)
)

allscores <- allscores[-removeme,]
saveRDS(allscores, subpath("CSCoV.scores.RDS"))
} else log("Score export skipped.")

compdata <- compdata[-removeme,]
saveRDS(compdata, subpath("compdata.RDS"))

## this is used by preprint_scoring to export compdata_ext_2 for ML models
saveRDS(removeme, subpath("dup_index.RDS"))

log("All done.")


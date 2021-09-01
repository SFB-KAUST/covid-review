
## OPTIONS

inoutpath <- "datanew"
logfile <- "topic_modeling.log"
# if F loads previous time-hungry results instead of computing them
doheavystuff <- T
errToFile <- T # redirect stderror to .err file

# !NOTE: The following needs to be done once
# ud_model <- udpipe_download_model(language="english")
# saveRDS(ud_model, "data/udpime_model_english.RDS")


## CHECK OPTIONS

subpath <- function(fname){
  return(file.path(inoutpath,fname))
}

compdata <- readRDS(subpath("compdata.RDS"))

checkpath <- function(p){
  if(!file.exists(p))
    stop(paste("Path does not exist:", p))
}

checkpath(inoutpath)
checkpath(subpath("figs"))
if(!file.exists(subpath("LDAdata")))
    dir.create(subpath("LDAdata"))



## START

library(dplyr)
library(ggplot2)
library(ggstatsplot)

log <- function(txt){
  txt <- paste0("[", Sys.time(), "] ", txt)
  write(txt, file=subpath(logfile), append=TRUE)
  print(txt)
}

if(errToFile)
  sink(file(subpath(paste0(logfile,".err")), open = "a"),
       type = "message", append = T) else sink(NULL, type="message")


## Extracting lemmas

library(udpipe)

annotate <- function(abs) {
  log("Annotating tokens...")
  library(foreach)
  library(doParallel)
  md <- readRDS(subpath("LDAdata/udpime_model_english.RDS"))
  #md$file_model = "/home/francesco/git/covid-papers/data/english-ewt-ud-2.4-190531.udpipe"
  md$file_model = subpath("LDAdata/english-ewt-ud-2.4-190531.udpipe")
  ud_model <- udpipe_load_model(md)

  registerDoParallel(26)
  i <- 1; batch <- 20
  res <- foreach(i=1:(ceiling(length(abs)/batch))) %dopar% {
    j <- min((i-1)*batch+1, length(abs))
    k <- min(j+batch-1, length(abs))
    #print(paste(format(Sys.time(), "%H:%M:%S"), j,k,length(abs)))
    as.data.frame(
      udpipe_annotate(ud_model, abs[j:k], doc_id=paste0(j:k))
    )[, c("doc_id", "lemma", "upos")]
  }
  res <- do.call(rbind, res)
  res <- res[order(as.numeric(res$doc_id)),]
  res$doc_id <- paste0("doc", res$doc_id)
  return(res)
}

if(doheavystuff) {
  log("Extracting lemmas...")
  comp_lemmas <- annotate(tolower(compdata$abstract))
  saveRDS(comp_lemmas, subpath("LDAdata/udpipe_compdata_annotations.RDS"))
} else {
  log("Loading extracted lemmas...")
  comp_lemmas <- readRDS(subpath("LDAdata/udpipe_compdata_annotations.RDS"))
}

res <- comp_lemmas
dtf <- subset(res, upos %in% c("NOUN", "ADJ", "VERB"))
dtf <- document_term_frequencies(dtf, "doc_id", term = "lemma")
log(paste("DTF size:", paste(dim(dtf), collapse="x")))

dtm <- document_term_matrix(dtf)
log(paste("DTM size:", paste(dim(dtm), collapse="x")))

dtm <- dtm_remove_sparseterms(dtm, remove_emptydocs = F)
dtm <- dtm_remove_terms(dtm, c("covid", "sar", "cov", "coronavirus",
                               "nlmcategory"), remove_emptydocs = F)
dtm <- dtm_remove_tfidf(dtm, cutoff=0.045,remove_emptydocs = F)
topic_outliers <- c(which(compdata$DOI==32511578), which(apply(dtm,1,sum)==0))

## EXTRACTING TOPICS

library(tidytext)
maketopics <- function(dtm, k, seed=NULL){
  if(is.null(seed))
    seed <- sample(1:1000,1)
  library(topicmodels)
  #m <- LDA(dtm, k, control=list(seed=1))
  m <- LDA(dtm, k, control=list(seed=seed))

  m_topics <- tidy(m, matrix = "beta")

  m_top_terms <- m_topics %>%
    group_by(topic) %>%
    top_n(30, beta) %>%
    ungroup() %>%
    arrange(topic, -beta)

  return(list(model=m, topterms=m_top_terms))
}




if(doheavystuff) {
 
  log("Building topic model...")
 registerDoParallel(26)
  res <- foreach(i=1:26) %dopar% {
   maketopics(dtm[-topic_outliers,], 6, seed=NULL)
  }

  ## Multiple runs give alternative solutions to be manually inspected.
  ## The following arbitrarily chooses the first one.
  
  res <- res[[1]]
  saveRDS(res, subpath("LDAdata/LDA_results.RDS"))
} else {
  res <- readRDS(subpath("LDAdata/LDA_results.RDS"))
  log("Loading topic model...")
}

lda_m <- res$model


log("Assigning topics...")

pred_data <- predict(lda_m, dtm)
compdata_ext <- cbind(compdata, topic=pred_data[,-c(1:5)])
ids <- as.numeric(gsub("doc","",rownames(dtm)))
pred <- data.frame(topic=pred_data$topic, Collection=compdata$collection[ids])
pred$topic <- res$renamer[pred$topic]
colnames(compdata_ext)[grep("topic", colnames(compdata_ext))] <- paste0("topic.",res$renamer)
compdata_ext <- cbind(compdata_ext, topic=pred$topic)


## Fixing papers uncorrectly classified as Imaging
## The model is mislead to classify as "Imaging" those papers that heavily use
## computational models that are typical of image analysis (such as CNNs). This
## is manually fixed in the following.
if(T) {
  imgset <- which(compdata_ext$topic=="Imaging")
  words <- apply(dtm[imgset, ], 1, function(x)names(x[x>0]))
  dontfix <- which(sapply(words, function(x) any(!is.na(
    grep("image|^ray$|^radiol|tomography|^ct",x, ignore.case=T)))))
  fixme <- imgset[-dontfix]
  temp <- compdata_ext[fixme, "topic.Imaging"]
  compdata_ext[fixme, "topic.Imaging"] <- 0
  maxtopic <- res$renamer[apply(compdata_ext[imgset,10:15], 1, which.max)]
  compdata_ext[fixme, "topic.Imaging"] <- temp
  compdata_ext[imgset, "topic"] <- maxtopic
  pred$topic[imgset] <- maxtopic
}


### Mapping articles to 2D embedding

log("Building Article map...")

if(doheavystuff) {
  log("Computing distances...")
  library(philentropy)
  library(umap)
  d <- distance(as.matrix(dtm[-topic_outliers,]), method = "jaccard")
  um <- umap(as.matrix(d))$layout
  saveRDS(um, subpath("LDAdata/topics_umap.RDS"))
} else {
  log("Loading distances...")
  um <- readRDS(subpath("LDAdata/topics_umap.RDS"))
}


## The following additional code exports data for the ML scripts

# data produced in previous
s_scholar <- readRDS(subpath("s_scholar_compdata.RDS"))
biomedstats <- readRDS(subpath("biomedrxiv_stats.RDS"))
compdata <- readRDS(subpath("compdata.RDS"))
scores <- readRDS(subpath("CSCoV.scores.RDS"))
dup_idx <- readRDS(subpath("dup_index.RDS"))

toDate <- function(dates)
  sort(as.Date(paste("28", dates), format="%d %b %Y"))
idx <- compdata$collection=="biorxiv"|compdata$collection=="medrxiv"

# Exporting for ML models
biomedstats <- biomedstats[rownames(biomedstats) %in% compdata$DOI,]
s_citations <- unlist(sapply(s_scholar, function(x) if(!is.na(x)) length(x$citations$doi) else 0))[-dup_idx]
s_influential <- unlist(sapply(s_scholar, function(x) if(!is.na(x)) x$influentialCitationCount else 0))[-dup_idx]
compdata_ext_2 <- cbind(compdata_ext, numcit=s_citations, influcit=s_influential)
compdata_ext_2 <- cbind(numcit=s_citations, influcit=s_influential,
                        "DL_Abstract"=NA, "DL_Full"=NA, "DL_Pdf"=NA)
rownames(compdata_ext_2) <- compdata$DOI
compdata_ext_2[rownames(biomedstats), c("DL_Abstract", "DL_Full", "DL_Pdf")] <- biomedstats
saveRDS(compdata_ext_2, subpath("compdata_ext_2.RDS"))
########################


log("All done.")


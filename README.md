# CSCoV: Computational Studies about COVID-19

This repository includes source codes to build the CSCoV database. The latest build can be directly downloaded from Zenodo at https://zenodo.org/record/5136684.

The code includes two parts:
- Data collection and categorization (in R)
- ML models (in Python).

## Data collection and categorization

- data_collection.R: collects papers from Pubmed, arXiv, bioRxiv, and medRxiv. Also gathers article and author-related information from Semantic Scholar.
- topic_modeling.R: analyzes articles to split them into topics.

## Machine Learning models

ML models will use data from R scripts to provide an additional score for preprint papers.

Order to run this project:
- 0.Doc2vec_embedding.ipynb
- 1.DeepWalk_embeddings.ipynb
- 2.Feature_selection.ipynb
- 3.Visualize_graph.ipynb
- 4.Deep_Learning.ipynb

Input files:
- data/compdata_cit_authmet.csv
- data/compdata_ref_author.csv

Input files to tune DL params:
- data/lda_doc2vec_feats.csv
- data/lda_doc2vec_targs.csv
- data/graph_embedding.csv

Install dependencies using command:
```
pip install -r requirements.txt
```

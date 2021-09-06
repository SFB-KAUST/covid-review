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

Install dependencies:
```shell
pip install -r requirements.txt
```
Order to run scripts:
- 0.LDA_topics.ipynb: evaluate topic numbers of LDA models.
- 1.Doc2vec_embedding.ipynb: retrieve Doc2Vec embeddings.
- 2.DeepWalk_embeddings.ipynb: get DeepWalk embeddings from the citation network.
- 3.Add_pub_metric_info.ipynb: manually verify the labels.
- 4.Deep_Learning.ipynb: train DL models and predict probability.
- 5.Visualize_graph.ipynb: get citation edges, node degrees and LDA topics for citation graph visualizaiton.


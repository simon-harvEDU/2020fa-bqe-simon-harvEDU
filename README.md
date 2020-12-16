# BiologyQueryEnhancer - BQE 
## Abstract
With the Biology Query Enhancer I am developing a tool for Biologists and Bioinformaticians to effectively search and query the two major biological and medical data bases for literature and knowledge: *pubmed* and *biorxiv*. The tool is made to query multiple questions to *both* data bases. The output can be either a csv data frame or a word cloud as a graphical output (see below). The search is quering the abstract and the text. Both tools are inspired by and using code from Biopython and https://github.com/blairbilodeau/arxiv-biorxiv-search. For the graphical output bigrams are used. The tool is ready to use with the programs necessary saved in the pipfile. 

![Alt text](/example_output.png "wordcloud search for p53 and cancer")



## 
example usage

```
pipenv shell 
python ./bqe/run_search.py
python ./bqe/run_wordcloud.py
```

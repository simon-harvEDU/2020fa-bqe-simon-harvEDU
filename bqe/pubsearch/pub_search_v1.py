import os, re
from Bio import Entrez
import time as time
import pandas as pd

# https://biopython.org/docs/1.75/api/Bio.Entrez.html


def search(query):
    Entrez.email = "your.email@example.com"
    handle = Entrez.esearch(
        db="pubmed", sort="relevance", retmax="10", retmode="xml", term=query
    )
    results = Entrez.read(handle)
    return results


def fetch_details(id_list):
    ids = ",".join(id_list)
    Entrez.email = "your.email@example.com"
    handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
    results = Entrez.read(handle)
    return results


def pubmedsearch(
    target_genes, search_terms, fast_search = True
):
    """ get the pubmed count or titles associated with genes and your search term. You can choose between output as count for the number of hits or the abstracts to get back 
    For EXAMPLE: search_terms(your_target_list = ["ABC2", "ABC1"], 
                                your_search_terms = ["fibrosis"], 
                                output = "count") """

    # formar and check input:
    # write an assertion statement here
    #if type(search_terms) != list:
    #    return "please provide the search terms as a list e.g.: ['heart']"

    sleep_time = 1
    if fast_search:
        sleep_time = 0.3

    res_d = {}
    res_abstract = {}

    print(search_terms)
    for term in search_terms:
        print(term)
    for term in search_terms:

        for target in target_genes:

            # define the search term:
            search_term = "{} AND {}".format(target, term)
            print(search_term)

            titles_list = []
            abstract_list = []
            time.sleep(sleep_time)

            # search the data base:
            results = search(search_term)  # query
            id_list = results[
                "IdList"
            ]  # list of UIDs - this will give us the article IDs
            chunk_size = 30  # how much data you want to read in one instance - there is a limit to get server answers

            ## add a timer here to know how many genes from the list were searched yet and how many more to go ##

            if (
                id_list
            ):  # if there something found for the query it should be in the ID list

                for chunk_i in range(
                    0, 60, chunk_size
                ):  # if there are more than 50 papers on the search query, we cap it there, this already shows that there is plenty of literature
                    chunk = id_list[chunk_i : chunk_i + chunk_size]

                    try:
                        # print("searchtermfired")
                        time.sleep(sleep_time)
                        papers = fetch_details(chunk)
                        for articles in papers["PubmedArticle"]:
                            # print(articles['PubmedData'].keys())
                            # print(articles["MedlineCitation"]["Article"]['Abstract']['AbstractText'])

                            abstract_list.append(
                                articles["MedlineCitation"]["Article"]["Abstract"][
                                    "AbstractText"
                                ]
                            )
                            # print(articles["MedlineCitation"]['Abstract'])

                            titles_list.append(
                                articles["MedlineCitation"]["Article"]["ArticleTitle"]
                            )

                        res_d[str(target) + "_" + term + "_title"] = titles_list
                        res_d[str(target) + "_" + term + "_abstract"] = abstract_list

                    except:  # occasionally a chunk might annoy your parser
                        time.sleep(sleep_time)
                        pass
            else:  # if there is nothing in the list
                res_d[str(target) + "_" + term + "_title"] = ["none"]
                res_d[str(target) + "_" + term + "_abstract"] = ["none"]
                
    res_d = pd.DataFrame.from_dict(res_d, orient="index")
    df = res_d.transpose()
    return df

class multi_pubsearch():
    def __init__(self, gene_list, combination_terms, max_records):
        self.max_records = max_records
        self.gene_list = gene_list
        self.combination_terms = combination_terms
        self.df = pd.DataFrame()
        
    def search_pubmed(self):
        self.df = pubmedsearch(target_genes = self.gene_list , search_terms = self.combination_terms)
        
        return self.df






#from pubsearch import pub_search_v1
from pubsearch import pub_search_v1 as ps
#from pub_search_v1 import multi_pubsearch

if __name__ == "__main__":
    res = ps.multi_pubsearch(gene_list = ["ABC1", "P53"] ,  combination_terms = ["fibrosis"], max_records = 2).search_pubmed()
    # out = pub_search.search_terms_beta(
    #     your_target_list=["ABC1", "P53"], your_search_terms=["fibrosis"]
    # )
    print(res)

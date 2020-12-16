# import biorxiv
from biorxiv_search import biorxiv

if __name__ == "__main__":
    print("module start")
    test = biorxiv.multi_medsearch(
        gene_list=["p53", "POSTN"], combination_terms=["fibrosis"], max_records=10
    ).search_biorxiv()
    print(test)

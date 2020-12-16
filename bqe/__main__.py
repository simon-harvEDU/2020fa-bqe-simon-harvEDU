#  import os
# from lits import try_class
# from lits import arx_easy
from luigi import build
from bqe import bqe_search

from bqe.bqe_search import combined_search, graphical_output


if __name__ == "__main__":
    build([combined_search(gene_list= ["p53", "CCND1"], combination_terms = ["apoptosis", "cell cycle"])], local_scheduler=True),
    # luigi --module combined_search graphical_output --gene-list '["p53", "CCND1"]' --combination_terms '["apoptosis", "cell cycle"]' --local-scheduler


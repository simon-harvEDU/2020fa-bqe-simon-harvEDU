#  import os
# from lits import try_class
# from lits import arx_easy
from luigi import build
import bqe_search

from bqe_search import combined_search, graphical_output

build([combined_search(gene_list= ["p53", "CCND1"], combination_terms = ["apoptosis", "cell cycle"])], local_scheduler=True),
build([graphical_output(gene_list= ["p53"], combination_terms = ["apoptosis"])], local_scheduler=True),
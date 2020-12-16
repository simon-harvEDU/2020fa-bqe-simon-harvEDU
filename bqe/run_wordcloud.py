#  import os
# from lits import try_class
# from lits import arx_easy
from luigi import build
import bqe_search

from bqe_search import graphical_output

build([graphical_output(gene_list= ["p53"], combination_terms = ["apoptosis"])], local_scheduler=True),
import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi import ExternalTask, Parameter, Task, run, LocalTarget

import pandas as pd
from biorxiv_search import biorxiv
from pubsearch import pub_search_v1
import datetime

from wordcloud import WordCloud, STOPWORDS 
import matplotlib.pyplot as plt 
from word_cloud_m import word_cloud 
import re

from csci_utils.luigi.target import SuffixPreservingLocalTarget

class luigi_biorxivsearch(Task):
    max_records = luigi.IntParameter(default=2)
    gene_list = luigi.ListParameter(default = ["p53", "POSTN"] )
    combination_terms= luigi.ListParameter(default = "cancer")

    def output(self):
        names_gene_list = "_".join(self.gene_list)
        names_comb_list = "_".join(self.combination_terms)
        return LocalTarget('./data/{}{}_biorxiv.csv'.format(str(names_gene_list),str(names_comb_list)))

    def run(self):
        df = biorxiv.multi_medsearch(gene_list = self.gene_list , combination_terms = self.combination_terms , max_records = self.max_records).search_biorxiv()

        f = self.output().open("w")
        df.to_csv(f, sep="\t", encoding="utf-8", index=None)
        f.close()


class luigi_pubmedsearch(Task):
    max_records = luigi.IntParameter(default=2)
    gene_list = luigi.ListParameter(default = ["p53", "POSTN"] )
    combination_terms= luigi.ListParameter(default = ["cancer"])

    def output(self):
        names_gene_list = "_".join(self.gene_list)
        names_comb_list = "_".join(self.combination_terms)
        return LocalTarget('./data/{}{}_pubmed.csv'.format(str(names_gene_list),str(names_comb_list)))

    def run(self):
        print(self.gene_list, self.combination_terms)
        df = pub_search_v1.multi_pubsearch(gene_list = self.gene_list ,  combination_terms = self.combination_terms , max_records = self.max_records).search_pubmed()
        f = self.output().open("w")
        df.to_csv(f, sep="\t", encoding="utf-8", index=None)
        f.close()


class combined_search(Task):
    gene_list = luigi.ListParameter(default = ["p53", "POSTN"] ) 
    combination_terms= luigi.ListParameter( default = ["cancer"] )
    max_records = luigi.IntParameter(default=2)

    def requires(self):
        return {"biorxiv": luigi_biorxivsearch(gene_list = self.gene_list, combination_terms = self.combination_terms, max_records = self.max_records), "pubmed": luigi_pubmedsearch(gene_list = self.gene_list, combination_terms = self.combination_terms, max_records = self.max_records)}

    def output(self):
        names_gene_list = "_".join(self.gene_list)
        names_comb_list = "_".join(self.combination_terms)
        return LocalTarget('./data/{}{}_combined.csv'.format(str(names_gene_list),str(names_comb_list)))

    def run(self):
        bio_r = self.input()['biorxiv'].open('r')
        df_bio = pd.read_csv(bio_r, sep='\t')
        df_bio.columns = [(str(i+'_bio') )for i in df_bio.columns]

        pubmed_r = self.input()['pubmed'].open('r')
        df_pub = pd.read_csv(pubmed_r, sep='\t')
        df_pub.columns = [(str(i+'_pub') )for i in df_pub.columns]
        
        # merge both df
        df_merged = pd.concat([df_pub, df_bio], axis=1)
        f1 = self.output().open("w")
        df_merged.to_csv(f1, sep='\t', encoding='utf-8', index=None)
        f1.close()


class graphical_output(Task):
    gene_list = luigi.ListParameter(default = ['p53'] ) 
    combination_terms= luigi.ListParameter( default = ['treatment'] )
    max_records = luigi.IntParameter(default=50)
    database = luigi.Parameter(default = 'pubmed')

    def requires(self):
        if self.database == 'pubmed':
            return luigi_pubmedsearch(gene_list = self.gene_list, combination_terms = self.combination_terms)
        elif self.database == 'biorxiv':
            return luigi_biorxivsearch(gene_list = self.gene_list, combination_terms = self.combination_terms)

    def output(self):
        target = SuffixPreservingLocalTarget(
            "./data/{}{}_wordcloud.png".format(self.gene_list[0], self.combination_terms[0]),
            format=luigi.format.Nop,
        )
        return target

    def run(self):
        # df processing
        df = self.input().open('r')
        df = pd.read_csv(df, sep='\t')
        df = df[[col for col in df.columns if 'abstract' in col]]
        c_name = str(self.gene_list[0]) +  '_' + str(self.combination_terms[0]) +  '_abstract'  
        test_words = ' '.join(df[c_name])

        strings_exclude = self.gene_list + self.combination_terms
        
        with self.output().temporary_path() as self.tmp_path:
            word_cloud.mk_wordcloud(test_words, filename_out = self.tmp_path ,strings_exclude = strings_exclude)
from luigi.contrib.external_program import ExternalProgramTask
from luigi import ExternalTask, Parameter, Task, run, LocalTarget

from arx_easy_v1 import biomedrxivsearch


#from lits import arx_easy_v1
import datetime
import try_class 

class biorxiv(Task):
    
    # biorxiv search
    search_term = Parameter(default = "p53")
    combination_term = Parameter(default = "apoptosis")

    # output
    out_path = Parameter(default="data_biorxiv/")

    def output(self):
        return biomedrxivsearch(
            kwd = [self.search_term, self.combination_term]
            )
    
    def run(self):


class pubmed(Task):
    search_term = Parameter(default = "p53")
    combination_term = Parameter(default = "apoptosis")

    out_path = Parameter(default="data_pubmed/")

    def requires(self):
        test =  try_class.my_search(self.search_term)
        a = test.search()
        b = test.fetch_details( id_list = a)

        return b
    
    def run(self):
        f = self.output().open("w")
        f.write("test")
        f.close()
    
    def output(self):
        return LocalTarget("./sometest_pubmed.txt")


    test =  try_class.my_search("p53")
    a = test.search()
    #print(a)
    b = test.fetch_details( id_list = a)


class combined(Task):
    def requires(self):
        return {"pubmed" : pubmed(),
        "biorxiv" : biorxiv()
        }
    
    def run(self):
        f = self.output().open("w")
        print  >> f, "something"
        f.close()
    
    def output(self):
        return LocalTarget("output_test.txt")

if __name__ == "__main__":
    run()

#     def run(self):
#         with self.input().open("r") as f1, self.output().open("w") as f2:
#             f2.write(f1.read())       

# class SavedModel(luigi.ExternalTask):
#     MODEL_ROOT = luigi.Parameter(default="s3://harvedu/pset_4/")
#     SHARED_RELATIVE_PATH = luigi.Parameter(default="saved_models/")
#     MODEL = luigi.Parameter(default="candy.pth")  # Filename of the model
#     binary = luigi.Parameter(default=True)

#     def output(self):

#         if self.binary:
#             return S3Target(
#                 self.MODEL_ROOT + self.SHARED_RELATIVE_PATH + str(self.MODEL),
#                 format=luigi.format.Nop,
#             )
#         else:
#             return S3Target(
#                 self.MODEL_ROOT + self.SHARED_RELATIVE_PATH + str(self.MODEL)
#             )


# class DownloadModel(luigi.Task):
#     S3_ROOT = luigi.Parameter(default="s3://harvedu/pset_4/")
#     SHARED_RELATIVE_PATH = luigi.Parameter(default="saved_models/")
#     S3_FILE = luigi.Parameter(default="candy.pth")

#     LOCAL_ROOT = os.path.abspath("data")

#     binary = luigi.Parameter(default=True)

#     def requires(self):
#         return SavedModel(
#             MODEL_ROOT=self.S3_ROOT,
#             SHARED_RELATIVE_PATH=self.SHARED_RELATIVE_PATH,
#             MODEL=self.S3_FILE,
#         )

#     def output(self):
#         if self.binary:
#             return luigi.LocalTarget(
#                 self.LOCAL_ROOT
#                 + "/"
#                 + str(self.SHARED_RELATIVE_PATH)
#                 + "/"
#                 + str(self.S3_FILE),
#                 format=luigi.format.Nop,
#             )
#         else:
#             return luigi.LocalTarget(
#                 self.LOCAL_ROOT
#                 + "/"
#                 + str(self.SHARED_RELATIVE_PATH)
#                 + "/"
#                 + str(self.S3_FILE)
#             )

#     def run(self):
#         with self.input().open("r") as f1, self.output().open("w") as f2:
#             f2.write(f1.read())


# class DownloadImage(luigi.Task):
#     S3_ROOT = luigi.Parameter(default="s3://harvedu/pset_4/")
#     SHARED_RELATIVE_PATH = luigi.Parameter(default="images/")
#     S3_FILE = luigi.Parameter(default="luigi.png")

#     LOCAL_ROOT = os.path.abspath("data")

#     binary = luigi.Parameter(default=True)

#     def requires(self):
#         return ContentImage(
#             IMAGE_ROOT=self.S3_ROOT,
#             SHARED_RELATIVE_PATH=self.SHARED_RELATIVE_PATH,
#             IMAGE=self.S3_FILE,
#         )

#     def output(self):
#         if self.binary:
#             return luigi.LocalTarget(
#                 self.LOCAL_ROOT
#                 + "/"
#                 + str(self.SHARED_RELATIVE_PATH)
#                 + "/"
#                 + str(self.S3_FILE),
#                 format=luigi.format.Nop,
#             )
#         else:
#             return luigi.LocalTarget(
#                 self.LOCAL_ROOT
#                 + "/"
#                 + str(self.SHARED_RELATIVE_PATH)
#                 + "/"
#                 + str(self.S3_FILE)
#             )

#     def run(self):
#         with self.input().open("r") as f1, self.output().open("w") as f2:
#             f2.write(f1.read())




# # coding: utf-8
# import luigi
# import os
# from luigi import ExternalTask, Parameter, Task, run
# from luigi.contrib.s3 import S3Target
# from csci_utils.luigi.target import SuffixPreservingLocalTarget
# from luigi.contrib.external_program import ExternalProgramTask
# from pset_4.tasks.data import DownloadImage, DownloadModel, SavedModel, ContentImage


# class Stylize(ExternalProgramTask):
#     model_name = luigi.Parameter(default="mosaic.pth")
#     model_path = luigi.Parameter(default="data/saved_models/")

#     image_name = luigi.Parameter(default="luigi.png")
#     image_path = luigi.Parameter(default="data/images/")

#     output_path = luigi.Parameter(default="data/output/")
#     output_img = luigi.Parameter(default="styl_images.png")

#     def requires(self):
#         return {
#             "image": DownloadImage(S3_FILE=self.image_name),
#             "model": DownloadModel(S3_FILE=self.model_name),
#         }

#     def output(self):
#         target = SuffixPreservingLocalTarget(
#             "{}{}{}".format(self.output_path, self.model_name, self.output_img),
#             format=luigi.format.Nop,
#         )
#         return target

#     def program_args(self):
#         return [
#             "python",
#             "-m",
#             "neural_style",
#             "eval",
#             "--content-image",
#             "{}{}".format(self.image_path, self.image_name),
#             "--model",
#             "{}{}".format(self.model_path, self.model_name),
#             "--output-image",
#             self.tmp_path,
#             "--cuda",
#             "0",
#         ]

#     def run(self):
#         with self.output().temporary_path() as self.tmp_path:
#             super().run()

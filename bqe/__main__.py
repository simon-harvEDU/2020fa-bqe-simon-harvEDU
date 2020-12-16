#  import os
# from lits import try_class
# from lits import arx_easy
from lits import arx_easy

# from arx_easy import biomedrxivsearch
# from arx_easy_v2 import biomedrxivsearch
import datetime
from lits import try_class

# import try_class


if __name__ == "__main__":
    records_df3 = arx_easy.biomedrxivsearch(
        start_date=datetime.date.today().replace(day=1),
        end_date=datetime.date.today(),
        journal="biorxiv",
        subjects=[],
        kwd=["domain", "Single-Cell"],
        kwd_type="all",
        abstracts=True,
        athr=[],
        max_records=50,
        max_time=300,
    )
    print(records_df3)

    test = try_class.my_search("p53")
    a = test.search()
    # print(a)
    b = test.fetch_details(id_list=a)
    print(b)

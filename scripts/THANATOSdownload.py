#!/usr/bin/python

import wget
import pandas as pd
from tqdm import tqdm, trange
from concurrent.futures import ThreadPoolExecutor

def get_tuples(data):
    return list(data[['url',"Accession"]].itertuples(index=False, name=None))

# def make_chunks(tup, n=1000):
#     """Yield successive n-sized chunks from lst."""
#     for i in range(0, len(tup), n):
#         yield lst[i:i + n]

def download_file(tup, directory="C://Users//sonja//REPOS//ShortProject2020//scripts//THANATOS//"):
    link, filekey = tup
    filename = directory + filekey + ".html"
    wget.download(link,filename)
#     print(link,filename)

def run_():

    tup = tqdm(get_tuples(data))
    with ThreadPoolExecutor() as exe:
        exe.map(download_file, tup)
        
#         for t in tqdm(tup):
#             res = [exe.submit(download_file,t) for t in tup]

#         for f in tqdm(concurrent.futures.as_completed(res), total=len(tup)):
#             temp = f.result()
# #             f.flush()
# #             os.fsync(f.fileno())



if __name__ == "__main__":
    #extract identifiers for each gene from its gene page
    #load file with links
    data = pd.read_csv('THANATOS_missing.tsv', sep='\t', encoding = 'utf-8')

    run_()

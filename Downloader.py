import os
import pandas as pd
import re
import requests
from TableLoader import TableLoader
import time
from tqdm import tqdm
import xml.etree.ElementTree as ET


class Downloader(object):
    wd = os.getcwd()
    root = None
    table = None
    pubmed_ids = set()
    pmc_ids = set()

    def __init__(self):
        # self.download_table()  # Downloads
        table_loader = TableLoader()
        self.table = table_loader.table
        # self.extract_ids()
        self.load_ids()
        self.delta_ids()
        self.download_content()

    def download_table(self):
        self.get_root()
        self.get_names()
        self.get_orphacodes()
        # self.download_orphanet()  # Downloads
        self.get_queries()
        # self.get_esearches()  # Downloads
        self.get_counts()
        self.save_table()

    def get_root(self):
        tree = ET.parse("en_product7.xml")
        self.root = tree.getroot()

    def get_names(self):
        name_list = list()
        for disorder in self.root.iter('Disorder'):
            name_list.append(disorder[2].text)
        name_list = pd.DataFrame(name_list)
        name_list.columns = ['Disease']
        self.table = name_list

    def get_orphacodes(self):
        orphacode_list = list()
        for disorder in self.root.iter('Disorder'):
            orphacode_list.append(disorder[0].text)
        orphacode_list = pd.DataFrame(orphacode_list)
        orphacode_list.columns = ['OrphaCode']
        self.table = self.table.merge(orphacode_list, left_index=True, right_index=True)

    def download_orphanet(self):
        print("Downloading Orphanet web pages...")
        for index in tqdm(self.table.index):
            orphacode = self.table.iloc[index, 1]
            url = f'https://www.orpha.net/consor/cgi-bin/OC_Exp.php?Expert={orphacode}'
            response = requests.get(url, allow_redirects=True)
            with open(self.wd + f"/orphanet/{orphacode}.html", "wb") as output:
                output.write(response.content)
        print("Done.")

    def get_queries(self):
        print("Getting the PubMed query of each disease...")
        query_list = list()
        for index in tqdm(self.table.index):
            orphacode = self.table.iloc[index, 1]
            url_search = open(self.wd + f"/orphanet/{orphacode}.html", "r")
            url_search = url_search.read()
            if "https://pubmed.ncbi.nlm.nih.gov/?term=" in url_search:
                url_search = re.sub(r'\n', '', url_search)
                url_search = re.sub(r'.+https://pubmed\.ncbi\.nlm\.nih\.gov/\?term=', '', url_search)
                url_search = re.sub(r'\'.+', '', url_search)
                query = url_search
                query_list.append(query)
            else:
                query_list.append(None)

        query_list = pd.DataFrame(query_list)
        query_list.columns = ['Query']
        self.table = self.table.merge(query_list, left_index=True, right_index=True)
        print("Done.")

    def get_esearches(self):
        print("Downloading esearch results...")
        for index in tqdm(self.table.index):
            query = self.table.iloc[index, 2]
            if query is not None:
                orphacode = self.table.iloc[index, 1]
                url_pubmed = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={query}&retmax=100000'
                response_pubmed = requests.get(url_pubmed, allow_redirects=True)
                with open(self.wd + f"/pubmed/{orphacode}-esearch.xml", "wb") as file_pubmed:
                    file_pubmed.write(response_pubmed.content)
                url_pmc = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=PMC&term={query}&retmax=100000'
                response_pmc = requests.get(url_pmc, allow_redirects=True)
                with open(self.wd + f"/pmc/{orphacode}-esearch.xml", "wb") as file_pmc:
                    file_pmc.write(response_pmc.content)
                time.sleep(0.4)
        print("Done.")

    def get_counts(self):
        print("Getting the PubMed article count of each disease...")
        count_list = list()
        for index in tqdm(self.table.index):
            if self.table.iloc[index, 2] is None:
                count_list.append(0)
            else:
                orphacode = self.table.iloc[index, 1]
                tree = ET.parse(self.wd + f"/pubmed/{orphacode}-esearch.xml")
                root = tree.getroot()
                count = root.find('Count').text
                count_list.append(count)

        count_list = pd.DataFrame(count_list)
        count_list.columns = ['Count']
        self.table = self.table.merge(count_list, left_index=True, right_index=True)
        print("Done.")

    def save_table(self):
        self.table.to_csv("rare_diseases.csv", index=False)

    def extract_ids(self):
        print("Extracting IDs from PubMed...")
        for index in tqdm(self.table.index):
            if self.table.iloc[index, 2] is not None:
                orphacode = self.table.iloc[index, 1]
                tree = ET.parse(self.wd + f"/pubmed/{orphacode}-esearch.xml")
                root = tree.getroot()
                article_ids = list()
                for ID in root.iter('Id'):
                    article_ids.append(ID.text)
                    self.pubmed_ids.add(ID.text)
                article_ids = pd.DataFrame(article_ids)
                article_ids.to_csv(self.wd + f"/pubmed/{orphacode}-ids.csv", index=False)
        self.pubmed_ids = list(self.pubmed_ids)
        self.pubmed_ids = pd.DataFrame(self.pubmed_ids)
        self.pubmed_ids.to_csv(self.wd + "/pubmed/all_ids.csv", index=False)
        print("Done.")

        print("Extracting IDs from PubMed Central...")
        for index in tqdm(self.table.index):
            if self.table.iloc[index, 2] is not None:
                orphacode = self.table.iloc[index, 1]
                tree = ET.parse(self.wd + f"/pmc/{orphacode}-esearch.xml")
                root = tree.getroot()
                article_ids = list()
                for ID in root.iter('Id'):
                    article_ids.append(ID.text)
                    self.pmc_ids.add(ID.text)
                article_ids = pd.DataFrame(article_ids)
                article_ids.to_csv(self.wd + f"/pmc/{orphacode}-ids.csv", index=False)
        self.pmc_ids = list(self.pmc_ids)
        self.pmc_ids = pd.DataFrame(self.pmc_ids)
        self.pmc_ids.to_csv(self.wd + "/pmc/all_ids.csv", index=False)
        print("Done")

    def load_ids(self):
        self.pubmed_ids = pd.read_csv(self.wd + "/pubmed/all_ids.csv")
        self.pmc_ids = pd.read_csv(self.wd + "/pmc/all_ids.csv")

    def delta_ids(self):
        self.pubmed_ids = self.pubmed_ids.squeeze()
        self.pubmed_ids = self.pubmed_ids.unique().tolist()
        self.pubmed_ids = set(map(str, self.pubmed_ids))
        for file in os.listdir(self.wd + "/pubmed/"):
            if file.endswith(".txt"):
                ID = file.removesuffix(".txt")
                self.pubmed_ids.remove(ID)
        self.pubmed_ids = list(self.pubmed_ids)
        self.pubmed_ids = pd.DataFrame(self.pubmed_ids)

        self.pmc_ids = self.pmc_ids.squeeze()
        self.pmc_ids = self.pmc_ids.unique().tolist()
        self.pmc_ids = set(map(str, self.pmc_ids))
        for file in os.listdir(self.wd + "/pmc/"):
            if file.endswith(".txt"):
                ID = file.removesuffix(".txt")
                self.pmc_ids.remove(ID)
        self.pmc_ids = list(self.pmc_ids)
        self.pmc_ids = pd.DataFrame(self.pmc_ids)

    def download_content(self):
        print("Downloading abstracts from PubMed...")
        for index in tqdm(self.pubmed_ids.index):
            pubmed_id = self.pubmed_ids.iloc[index, 0]
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pubmed_id}&retmode=xml&api_key=15b8d0248116528f9ed88b7e47796350b108'
            response = requests.get(url, allow_redirects=True)
            tree = ET.ElementTree(ET.fromstring(response.text))
            root = tree.getroot()
            for abstract in root.iter("AbstractText"):
                with open(self.wd + f"/pubmed/{pubmed_id}.txt", "w", encoding="utf8") as output:
                    abstract = " ".join(abstract.itertext())
                    output.write(abstract)
            time.sleep(0.1)
        print("Done.")

        print("Downloading articles from PubMed Central...")
        for index in self.pmc_ids.index:
            pmc_id = self.pmc_ids.iloc[index, 0]
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=PMC&id={pmc_id}&retmode=xml&api_key=15b8d0248116528f9ed88b7e47796350b108'
            response = requests.get(url, allow_redirects=True)
            tree = ET.ElementTree(ET.fromstring(response.text))
            root = tree.getroot()
            for body in root.iter("body"):
                with open(self.wd + f"/pubmed/{pmc_id}.txt", "w", encoding="utf8") as output:
                    body = " ".join(body.itertext())
                    output.write(body)
            time.sleep(0.1)
        print("Done.")

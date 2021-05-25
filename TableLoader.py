import os
import pandas as pd


class TableLoader:
    wd = os.getcwd()
    table = None

    def __init__(self):
        self.table = pd.read_csv(self.wd + "/rare_diseases.csv")
        self.table = self.table.where(pd.notnull(self.table), None)

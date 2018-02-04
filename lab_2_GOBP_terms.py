import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


name_file_to_open = "Cell-Cycle-Set.xlsx"

def open_file(name):



	xl = pd.ExcelFile("Cell-Cycle-Set.xlsx")


	df = xl.parse("Sheet1")


	df = df.dropna()
	df = df.reset_index(drop=True)


	return df

df = open_file(name_file_to_open)


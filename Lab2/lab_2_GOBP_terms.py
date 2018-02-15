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


print (df.head())
GOBP_val_count = {}

# for i in range(df.shape[0]):
# 	print (i)


for row_GOBP in df.GOBP:
	row_GOBP = row_GOBP.split(';')
	for item in row_GOBP:
		if item in GOBP_val_count:
			GOBP_val_count[item] +=1
		else:
			GOBP_val_count[item] = 1

print (GOBP_val_count)



# for row_GOBP in df.GOBP:
# 	row_GOBP = row_GOBP.split(';')
# 	for item in row_GOBP:
# 		GOBP_val_count[item] +=1

# print (GOBP_val_count)


# d = {'key':0, 'rand': 0}

# if 'rand' in d:
# 	print (d)
# 	d['ff'] = 0
# 	print (d)

# if 'ff' in d:
# 	d['ff'] += 1
# 	print (d)


#we need to count
#create a dictionary of strings iteger pair,add every turn that is not in list 
#(this generates a dictionary with all terms and value = 0)
#for every term in row_gobp (iterate across dictionary, if found add 1)
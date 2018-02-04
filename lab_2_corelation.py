import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


name_file_to_open = "Cell-Cycle-Set.xlsx"

def open_file(name):



	xl = pd.ExcelFile("Cell-Cycle-Set.xlsx")

	print(xl.sheet_names)

	df = xl.parse("Sheet1")


	df = df.dropna()
	df = df.reset_index(drop=True)
	print (df.head())

	print (df.shape)
	print (df.GOBP[2])

	return df


def GOXX_row_finder(df, word_finder, GOXX):
	list_of_rows_with_cell_cycle = []

	for row_index in range((df.shape)[0]):

		row_GOXX_all = df[GOXX][row_index].split(';')

		for item in row_GOXX_all:

			if str(item) == word_finder:
				print (row_index,item)
				list_of_rows_with_cell_cycle.append(row_index)



	print (list_of_rows_with_cell_cycle)

	df = df.iloc[list_of_rows_with_cell_cycle]

	print (df.head())
	return df



df = open_file(name_file_to_open)

df_cell_cycle = GOXX_row_finder(df, 'cell cycle' , 'GOBP' )
df_ribosome = GOXX_row_finder(df, 'ribosome', 'GOCC')



#weeker corelation
fig, axis = plt.subplots(1,2, figsize=(11,4))
axis[0].scatter(df_cell_cycle.mean_RNA_G1, df_cell_cycle.mean_protein_G1, c='r')
axis[0].scatter(df_cell_cycle.mean_RNA_S, df_cell_cycle.mean_protein_S, c='g')
axis[0].scatter(df_cell_cycle.mean_RNA_G2, df_cell_cycle.mean_protein_G2, c="b")

axis[1].scatter(df_ribosome.mean_RNA_G1, df_ribosome.mean_protein_G1, c='r')
axis[1].scatter(df_ribosome.mean_RNA_S, df_ribosome.mean_protein_S, c='g')
axis[1].scatter(df_ribosome.mean_RNA_G2, df_ribosome.mean_protein_G2, c="b")



fig, axis = plt.subplots(2,3)
corr_G1_CC = df_cell_cycle[['mean_RNA_G1','mean_protein_G1']].corr()
corr_S_CC = df_cell_cycle[['mean_RNA_S', 'mean_protein_S']].corr()
corr_G2_CC = df_cell_cycle[['mean_RNA_G2', 'mean_protein_G2']].corr()

axis[0][0].matshow(corr_G1_CC)
axis[0][1].matshow(corr_S_CC)
axis[0][2].matshow(corr_G2_CC)

print ('corr_G1_CC: \n', corr_G1_CC)
print ('corr_S_CC: \n', corr_S_CC)
print ('corr_G2_CC: \n', corr_G2_CC)

# stronger corelation


corr_G1_R = df_ribosome[['mean_RNA_G1','mean_protein_G1']].corr()
corr_S_R = df_ribosome[['mean_RNA_S', 'mean_protein_S']].corr()
corr_G2_R = df_ribosome[['mean_RNA_G2', 'mean_protein_G2']].corr()

axis[1][0].matshow(corr_G1_R)
axis[1][1].matshow(corr_S_R)
axis[1][2].matshow(corr_G2_R)

print ('corr_G1_R: \n', corr_G1_R)
print ('corr_S_R: \n', corr_S_R)
print ('corr_G2_R: \n', corr_G2_R)





plt.show()


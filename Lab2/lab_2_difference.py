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


def diff_cell_cycle(df, nameX,nameX2):

	Diff_X = df[nameX] - df[nameX2]
	Diff_X -= Diff_X.mean()
	Diff_X /= Diff_X.var()
	return Diff_X



def GOXX_row_finder(df, word_finder, GOXX):
	list_of_rows_with_cell_cycle = []

	for row_index in range((df.shape)[0]):

		row_GOXX_all = df[GOXX][row_index].split(';')

		for item in row_GOXX_all:

			if str(item) == word_finder:
				print (row_index,item)
				list_of_rows_with_cell_cycle.append(row_index)



	df = df.iloc[list_of_rows_with_cell_cycle]

	print (df.head())
	return df

def plot_difference_given_df(df,Title):

	Diff_G1_S_R = diff_cell_cycle(df,'mean_RNA_S','mean_RNA_G1')
	Diff_S_G2_R = diff_cell_cycle(df,'mean_RNA_G2','mean_RNA_S')
	Diff_G2_G1_R = diff_cell_cycle(df,'mean_RNA_G1','mean_RNA_G2')
	# we want to find difference, and for each take mean away and subtract by variance


	Diff_G1_S_P = diff_cell_cycle(df,'mean_protein_S','mean_protein_G1')
	Diff_S_G2_P = diff_cell_cycle(df,'mean_protein_G2','mean_protein_S')
	Diff_G2_G1_P = diff_cell_cycle(df,'mean_protein_G1','mean_protein_G2')



	fig, axis = plt.subplots(1,1)
	axis.scatter(Diff_G1_S_R,Diff_G1_S_P, color= ['springgreen'], alpha = 0.3)
	axis.scatter(Diff_S_G2_R,Diff_S_G2_P, color= ['red'], alpha = 0.3 )
	axis.scatter(Diff_G2_G1_R,Diff_G2_G1_P, color= ['blue'],  alpha = 0.3)
	axis.set_title(Title)
	axis.set_xlabel('RNA')
	axis.set_ylabel('Protein')

	fig, axis = plt.subplots(1,3, figsize = (15, 3))
	axis[0].plot([0,0], [-5,5], linestyle = ':')
	axis[0].plot([-5,5], [0,0], linestyle = ':')
	axis[0].scatter(Diff_G1_S_R,Diff_G1_S_P, color= ['springgreen'], alpha = 0.3)
	axis[1].plot([0,0], [-5,5], linestyle = ':')
	axis[1].plot([-5,5], [0,0], linestyle = ':')
	axis[1].scatter(Diff_S_G2_R,Diff_S_G2_P, color= ['red'], alpha = 0.3 )
	axis[2].plot([0,0], [-5,5], linestyle = ':')
	axis[2].plot([-5,5], [0,0], linestyle = ':')
	axis[2].scatter(Diff_G2_G1_R,Diff_G2_G1_P, color= ['blue'],  alpha = 0.3)



	green_patch = mpatches.Patch(color='green', label='Diff_G1_S')
	red_patch = mpatches.Patch(color='red', label='Diff_S_G2')
	blue_patch = mpatches.Patch(color='blue', label='Diff_G2_G1')
	plt.legend(handles=[green_patch,red_patch, blue_patch])




plot_ribo_cell_cycle = True
		
word_to_find, GOXX = 'intracellular part', 'GOCC'
df = open_file(name_file_to_open)
import matplotlib.patches as mpatches

if plot_ribo_cell_cycle == True:
	df_cell_cycle = GOXX_row_finder(df, 'cell cycle' , 'GOBP' )
	df_ribosome = GOXX_row_finder(df, 'ribosome', 'GOCC')

	plot_difference_given_df(df, 'Difference all')
	plot_difference_given_df(df_cell_cycle, 'Difference cell cycle')
	plot_difference_given_df(df_ribosome, 'Difference ribosome')

else:
	df_X = GOXX_row_finder(df, word_to_find, GOXX)
	plot_difference_given_df(df_X, 'Difference ' + word_to_find) 

	





print (df.head())

plt.show()

# .str.contains('')



# Y = X - Mx / max(X)-min(X)


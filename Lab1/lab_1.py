import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


xl = pd.ExcelFile("Cell-Cycle-Set.xlsx")

print(xl.sheet_names)

df = xl.parse("Sheet1")

print (df.head())


df = df.dropna()


# print (df.head())

# print (np.array(df.mean_RNA_G1))

fig, axis = plt.subplots(1,1) 
axis.hist(df.mean_RNA_G1)
axis.set_title("RNA_G1 Hist")
axis.set_xlabel('Distribution')
axis.set_ylabel('Frequency')
axis.hist(df.mean_protein_G1)
axis.set_title("Protein_G1 Hist")
axis.set_xlabel('Distribution')
axis.set_ylabel('Frequency')



fig, axis = plt.subplots(1,1) 
axis.hist(df.mean_RNA_S)
axis.set_title("RNA_S Hist")
axis.set_xlabel('Distribution')
axis.set_ylabel('Frequency')
axis.hist(df.mean_protein_S)
axis.set_title("Protein_S Hist")
axis.set_xlabel('Distribution')
axis.set_ylabel('Frequency')

fig, axis = plt.subplots(1,1) 
axis.hist(df.mean_RNA_G2)
axis.set_title("RNA_G2 Hist")
axis.set_xlabel('Distribution')
axis.set_ylabel('Frequency')
axis.hist(df.mean_protein_G2)
axis.set_title("Protein_G2 Hist")
axis.set_xlabel('Distribution')
axis.set_ylabel('Frequency')


fig = plt.figure()
plt.hist(df.mean_RNA_G1)
plt.hist(df.mean_protein_G1)
plt.hist(df.mean_RNA_S)
plt.hist(df.mean_protein_S)
plt.hist(df.mean_RNA_G2)
plt.hist(df.mean_protein_G2)
# plt.set_title("Hist")
# plt.set_xlabel('Distribution')
# plt.set_ylabel('Frequency')
















corr_all = df.corr()
corr_protein = df[['mean_RNA_G1',  'mean_RNA_S',  'mean_RNA_G2']].corr()
corr_RNA = df[['mean_protein_G1' , 'mean_protein_S',  'mean_protein_G2' ]].corr()
corr_RNA = df[['mean_protein_G1' , 'mean_protein_S',  'mean_protein_G2' ]].corr()
corr_G1 = df[['mean_RNA_G1',  'mean_protein_G1']].corr()
corr_S = df[['mean_RNA_S',  'mean_protein_S']].corr()
corr_G2 = df[['mean_RNA_G2',  'mean_protein_G2']].corr()

print ('All: \n',np.array((corr_all)))
print('Protein: \n',np.array(corr_protein))
print ('RNA: \n',np.array(corr_RNA))
print ('G1: \n', corr_G1)
print ('S: \n', corr_S)
print ('G2: \n', corr_G2)


fig, axis = plt.subplots(1,3,figsize=(15,7 ))
axis[0].scatter(df.mean_RNA_G1, df.mean_protein_G1, c="r")
axis[1].scatter(df.mean_RNA_S, df.mean_protein_S, c="g")
axis[2].scatter(df.mean_RNA_G2, df.mean_protein_G2)

m,b = np.polyfit(df.mean_RNA_G1, df.mean_protein_G1, 1)
# m2,m,b = np.polyfit(df.mean_RNA_G1, df.mean_protein_G1, 2)
# print (a)
# print (b)

axis[0].plot(df.mean_RNA_G1,df.mean_RNA_G1*m+b)

# axis[0].plot(df.mean_RNA_G1,df.mean_RNA_G1*m2**2+m*df.mean_RNA_G1+ b)
# print (df.Gene_Name)

plt.show()	
import numpy as np
from string import *
import pandas as pd

#-------------------------------------------------------
#This function returns to values for cae of match or mismatch
def Diagonal(n1,n2,pt):
    if(n1 == n2):
        return pt['MATCH']
    else:
        return pt['MISMATCH']

#------------------------------------------------------------   
#This function gets the optional elements of the aligment matrix and returns the elements for the pointers matrix.
def Pointers(di,ho,ve):

    pointer = max(di,ho,ve) #based on python default maximum(return the first element).
    print (di, ho, ve)

    if(di == pointer):
        return 'D'
    elif(ho == pointer):
        return 'H'
    else:
         return 'V'    


blosum50_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
blosum50 = pd.read_csv('blosum50.txt', names = blosum50_amino_acids, header=None, sep=',')

def blosum50_extract(blosum,letter1, letter2):
	# find the index of the second letter,
	index_of_second_letter = blosum.columns.get_loc(letter2)
	return (blosum[letter1][index_of_second_letter])
	#print (blosum50.columns.get_loc('N'))
	#print (blosum50['N'][2])
#--------------------------------------------------------
#This function creates the aligment and pointers matrices
def NW(s1,s2,match = 4,mismatch = -4, gap = -8):
    penalty = {'MATCH': match, 'MISMATCH': mismatch, 'GAP': gap} #A dictionary for all the penalty valuse.
    length_p1 = len(s1) + 1 #The dimension of the matrix columns.
    length_p2  = len(s2) + 1 #The dimension of the matrix rows.
    al_mat = np.zeros((length_p2 ,length_p1),dtype = int) #Initializes the alighment matrix with zeros.
    p_mat = np.zeros((length_p2 ,length_p1),dtype = str) #Initializes the alighment matrix with zeros.
    #Scans all the first rows element in the matrix and fill it with "gap penalty"
    for i in range(length_p2 ):
        al_mat[i][0] = penalty['GAP'] * i
        p_mat[i][0] = 'V'
    #Scans all the first columns element in the matrix and fill it with "gap penalty"
    for j in range (length_p1):
        al_mat[0][j] = penalty['GAP'] * j
        p_mat [0][j] = 'H'
    #Fill the matrix with the correct values.

    p_mat [0][0] = 0 #Return the first element of the pointer matrix back to 0.
    for i in range(1,length_p2 ):
        for j in range(1,length_p1):
            di = al_mat[i-1][j-1] + blosum50_extract(blosum50,s1[j-1],s2[i-1])  # Diagonal(s1[j-1],s2[i-1],penalty) #The value for match/mismatch -  diagonal.
            ho = al_mat[i][j-1] + penalty['GAP'] #The value for gap - horizontal.(from the left cell)
            ve = al_mat[i-1][j] + penalty['GAP'] #The value for gap - vertical.(from the upper cell)
            al_mat[i][j] = max(di,ho,ve) #Fill the matrix with the maximal value.(based on the python default maximum) # if the value is negative, just put 0.
            p_mat[i][j] = Pointers(di,ho,ve)
    print (np.matrix(al_mat))
    print (np.matrix(p_mat))

A = NW('HEAGAWGHEE', 'PAWHEAE')

#https://www.biostars.org/p/231391/
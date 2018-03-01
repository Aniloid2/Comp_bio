import numpy as np
import pandas as pd
from io import StringIO


blosum50_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
blosum50 = pd.read_csv('blosum50.txt', names = blosum50_amino_acids, header=None, sep=',')

def blosum50_extract(blosum,letter1, letter2):
	# find the index of the second letter,
	index_of_second_letter = blosum.columns.get_loc(letter2)
	return (blosum[letter1][index_of_second_letter])
	#print (blosum50.columns.get_loc('N'))
	#print (blosum50['N'][2])




def define_HW_value(lh, lv,ld, BV):
	# find if we have a match
	value = max(ld + BV,lh, lv )
	direction = ''
	if (value == ld+BV):
		direction = 'D'
	elif(value == lh):
		direction = 'H'
	elif(value == lv):
		direction = 'V'
	else:
		direction = 'ERROR'

	return (value, direction)



blosum50_match_mismatch = blosum50_extract(blosum50, 'N', 'N')

# so we are given two strings, we have to create a empy matrix that will hold the value, and another matrix that will hold the directions
def Need_wun(protein1, protein2):
	gap = -8 # default gap
	length_p1 = len(protein1) + 1
	length_p2 = len(protein2) + 1
	value_matrix = np.zeros((length_p2, length_p1))
	direction_matrix= np.zeros((length_p2, length_p1) ,dtype=str)
	# generators for both leftest and rightest colum
	Horizontal_gap = [0] + [ gap + gap*x for x in range(length_p2-1)]
	Vertical_gap = [0] + [ gap + gap*x for x in range(length_p1-1)]
	value_matrix[0,:] = Vertical_gap
	value_matrix[:,0] = Horizontal_gap
	direction_matrix[:,0] = 'V'
	direction_matrix[0,:] = 'H'
	direction_matrix[0,0] = '0'


	for i in range(1, length_p2):
		for j in range(1,length_p1):

			# for every value we want to find the H + gap,V + gap and D + blosum_extract
			local_v = value_matrix[i-1,j] + gap
			local_h = value_matrix[i, j-1]  + gap
			local_d = value_matrix[i-1,j-1]
			#get two letters from the strings that we need for this iteration
			H_letter = protein2[i-1]
			V_letter = protein1[j-1]
			#print (H_letter, V_letter)
			blosum_value = blosum50_extract(blosum50,H_letter,V_letter)
			#print (blosum_value)
			#print (local_h, local_v, local_d, blosum_value)
			value, direction = define_HW_value(local_h, local_v, local_d, blosum_value)
			# update the current matrix we are on
			value_matrix[i][j] = value
			direction_matrix[i][j] = direction

	return (value_matrix, direction_matrix)







def Back_algo_Need_W(value_matrix, direction_matrix, proteinj, proteini):

	counter_i, counter_j = value_matrix.shape[0] , value_matrix.shape[1]

	AlignmentA, AlignmentB = [],[]
	i = value_matrix.shape[0] 
	j = value_matrix.shape[1] 


	while (i >0 or j >0):
		local_direction = direction_matrix[i -1, j - 1]
		if (i > 0 and j>0 and local_direction == 'D'):
			AlignmentA.append(proteini[i-2])
			AlignmentB.append(proteinj[j-2])
			# D
			i -=1
			j -=1
		elif (i > 0 and local_direction == 'V' ):
			# V
			AlignmentA.append('-')
			AlignmentB.append(proteini[i-2])
			i -=1
		elif (j >0 and local_direction == 'H' ):
			# H 
			AlignmentB.append('-')
			AlignmentA.append(proteinj[j-2])
			j -=1
		elif (local_direction == '0'):
			break


	AlignmentA.reverse()
	AlignmentB.reverse()

	return (AlignmentA, AlignmentB)




value_matrix, direction_matrix = Need_wun('HEAGAWGHEE', 'PAWHEAE')


A, B = Back_algo_Need_W(value_matrix, direction_matrix, 'HEAGAWGHEE', 'PAWHEAE')

print (" Alignment for ('HEAGAWGHEE', 'PAWHEAE'): \n", A, '\n', B)

value_matrix, direction_matrix = Need_wun('PQPTTPVSSFTSGSMLGRTDTALTNTYSAL', 'PSPTMEAVEASTASHPHSTSSYFATTYYHL')

A, B = Back_algo_Need_W(value_matrix, direction_matrix, 'PQPTTPVSSFTSGSMLGRTDTALTNTYSAL', 'PSPTMEAVEASTASHPHSTSSYFATTYYHL')

print (" Alignment for ('PQPTTPVSSFTSGSMLGRTDTALTNTYSAL', 'PSPTMEAVEASTASHPHSTSSYFATTYYHL'): \n", A,'\n', B)


# mat = np.matrix(direction_matrix)
# with open('lab3_Need_wun_seq.txt','wb') as f:
#     for line in mat:
#         np.savetxt(f, line, fmt='%.2s')
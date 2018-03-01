
import numpy as np
from math import log
from matplotlib import pyplot as plt

import xlsxwriter



class State():
	def __init__(self, Emmis, Trans):
		self.Trans = Trans
		self.Emmis = Emmis

	def print_states(self):
		print (self.Trans, self.Emmis)

	def turn_to_logs(self):
		for i,j in self.Emmis.items():
			self.Emmis[i] = log(j)

		for i in range(len(self.Trans)):
			temp = log(self.Trans[i])
			self.Trans[i] = temp





class HMM():
	def __init__(self, sequence, states):
		self.sequence = sequence
		self.states = states

	def viterbi(self):

			# initialising the states to be in log form for no overflow
			for i in range(len(self.states)):
				self.states[i].turn_to_logs()
			
			# initialising V and pointer
			globalV = []
			globalV.append([1] + [0 for i in range(len(self.states))])

			Global_pointer = [[None for i in range(len(self.states))]]	
			
			

			# bulk of viterbi 
			for t in range(1,len(self.sequence)): 
				# a local probability &  pointer which will later be appended
				PL = [0 for i in range(len(self.states))] 
				local_pointer = [0 for i in range(len(self.states))] 
				for i in range(len(self.states)):
					# under the current state find the maximum probability V(t - 1)+Aij which will lead to next state by looking at the previeous state
					# and append it, then time/add by the emissional prob
					localV = []
					for j in range(len(self.states)):
						transmiss = self.states[i].Trans[j] 

						Vj = globalV[t-1][j] 

						localV.append(Vj + transmiss)

					emiss = self.states[i].Emmis[self.sequence[t]]
					# keep track which one of the two was the maximum by assinging an index pointer
					PL[i] = emiss + max(localV)
					local_pointer[i] = localV.index(max(localV))
						

 

				globalV.append(PL)

				Global_pointer.append(local_pointer)




			# find which state had the highest value
			last_entry_max_index = int(globalV[-1].index(max(globalV[-1])))

			Final_seq = []

			# go back wards to find the states which result in maximum probability sequence, backward algo

			for t in range(len(self.sequence)-1,-1, -1):

				Final_seq.append(last_entry_max_index)

				last_entry_max_index = Global_pointer[t][int(last_entry_max_index)]

			# reverse it
			R_Seq = Final_seq[::-1]

			# plot it 
			fig, axis = plt.subplots()
			axis.plot([t for t in range(len(self.sequence) )], Final_seq[::-1], 'x')
			axis.set_title('Sequence Graphed')
			axis.set_xlabel('Sequence No')
			axis.set_ylabel('State belonging')

			# put it in color
			workbook = xlsxwriter.Workbook('Phase_lambda_sequenced.xlsx')
			worksheet = workbook.add_worksheet()

			csv_file =[]
			for i in range(0, len(self.sequence), 10):
				csv_line = []
				for j in range(i, i+10):
					try:

						csv_line.append(self.sequence[j] )

					except:
						pass

				csv_file.append(csv_line)


			for i in range(len(csv_file)):
				for j in range(len(csv_file[i])):

					if R_Seq[j + 10*i] == 0:

						F = workbook.add_format({'font_color': 'red'})
						worksheet.write(i, j, self.sequence[j + 10*i], F)

					else:

						F = workbook.add_format({'font_color': 'blue'})
						worksheet.write(i,j, self.sequence[j + 10*i], F)

			workbook.close()

		

SeQ = list(" 5453525456666664365666635661416626365666621166211311155566351166565663466653642535666662541345464155")


state1 = State({'1': 1/6, '2' : 1/6, '3': 1/6, '4': 1/6, '5': 1/6, '6':1/6}, [9/10, 1/10])
state2 = State({'1': 1/10, '2': 1/10, '3':1/10, '4':1/10, '5':1/10, '6': 1/2}, [1/10, 9/10])



dosonest_casino = HMM((SeQ), (state1, state2))
dosonest_casino.viterbi()

Phase_L = [' '] + [nucleo for nucleo in open('PL.txt').read() if nucleo != '\n' if nucleo != ' ']

state3 = State({'A' : 0.2698, 'T' : 0.3237, 'C' : 0.2080, 'G': 0.1985}, [0.9998, 0.0002])
state4 = State({'A' : 0.2459, 'T': 0.2079 , 'C' : 0.2478, 'G' : 0.2984}, [0.0003, 0.9997])



ritch_region = HMM(Phase_L, (state3, state4))
ritch_region.viterbi()


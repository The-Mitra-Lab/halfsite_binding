"""
This is a toolkit to compute the free energy of interactions between Gal4 and a sequence of DNA

NJM
V6
6/29/17


v2 This version changes the dictionary from 0 to 16 and changes it to a dict of the free energy
v3 Now the input is from a text file, and can be as long as the text file.  Designed to be independed from sequence
v4 Updated to be more similar in style to the pssm class. Text file read in moved to a separate generator file.
   Max, min, and search functions added.  Search combines the calculate and cutOff methods and can calculate with reverse complement
v5 Text file input is now in the form of get_psfem("tf_name") and has been returned to this file.  Condensed matrix files to also inlcude
   the sequence, kd, and recommended threshold.  To load a specific psfem, use get_psfem("tf_name"), provided a psfem text file of the tf is 
   in the TFs folder.
v6 ROC Curve functions from analyze_sig_hits_v2 have been added
"""



from math import log,e
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from sklearn import metrics
import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics

class psfem:
	"""
	Requires the base Kd used for the sequence and the file name of the free energy matrix
	infoFile - file with kd and sequence
	inputFile - file with psfem
	"""
	def __init__(self, sequence, baseKd,psfeMatrix,recommended_threshold = 0.0): 
		"""Constants"""
		self.REF_CONCENTRATION = 1.0
		self.BOLTZMANN = 1.380*pow(10,-23)
		self.T = 303
		self.consensus_sequence = Seq(sequence,IUPAC.unambiguous_dna)
		self.baseKd = float(baseKd)
		self.consensus_fe = log(self.baseKd*self.REF_CONCENTRATION) 
		self.recommended_threshold = recommended_threshold

		self.psfem = {"C":[],
				        "T":[],
				        "A":[],
				        "G":[]
				   }
		self.psfem = psfeMatrix

		
		self.length = len(self.psfem["C"])
		
		maximum = self.maxE()
		self.max_fe = maximum[0]
		self.max_sequence = Seq(maximum[1],IUPAC.unambiguous_dna)

		minimum =  self.minE()
		self.min_fe = minimum[0]
		self.min_sequence = Seq(minimum[1],IUPAC.unambiguous_dna)

    #not sure if implemented correctly
	def __getitem__(self, letter):
		return self.psfem[letter]

	"""Functions"""
	#Calculates free energy from a given Kd
	def calcEnergy(self,Kd):
		energy = log(Kd*self.REF_CONCENTRATION)# energy in units of KbT
		return energy

	#Calculates the change in free energy for a specific relative Affinity (used for each changed base)  USED FOR PSFEM MATRIX
	def calcDeltaDeltaG(self,relativeAffinity): #calculate for a single index
		newKd = (1.0/relativeAffinity)*self.baseKd
		return self.calcEnergy(newKd)-self.consensus_fe

	#Calculates the sum of the changes in free energy from all of the changed base pairs
	#Free energy of a sequence
	def calcDeltaG(self,my_seq): #add all of the single indexes
		sumDDG = 0.0
		for x in range(0,self.length):
			cursor = my_seq[x:x+1]
			if(cursor == ''):
				break
			ddg = self.psfem[str(cursor)][x]
			sumDDG = sumDDG + ddg
		return sumDDG

	"""From pssm program"""
	def reverse_complement(self): 
		values = {} 
		values["A"] = self.psfem["T"][::-1] 
		values["T"] = self.psfem["A"][::-1] 
		values["G"] = self.psfem["C"][::-1] 
		values["C"] = self.psfem["G"][::-1]

		return self.__class__(str((self.consensus_sequence).reverse_complement()), self.baseKd, values, self.recommended_threshold) 
	
	#Returns a list of final free energies for each sequence starting from 0 to n-length of base sequence
	def calculate(self, sequence):
		FinalGList = []
		FinalG = 0
		for x in range(0,len(sequence)-self.length+1):
			if(len(sequence) >= self.length):
				inputSeq = str(sequence[x:x+self.length])
			else:
				print("Missed")
				break
			finalDDG = self.calcDeltaG(inputSeq)
			FinalG = self.consensus_fe + finalDDG
			FinalGList.append(FinalG)
		if(len(FinalGList)==1):
			return FinalGList[0]
		else:
			return FinalGList
	#									   couldn't do recommended_threshold, idk why
	def search(self, sequence, threshold = -5.0, both_strands = True):

		n = len(sequence)
		m = self.length
		if(both_strands):
			rc = self.reverse_complement()

		orientation = 1
		for x in range(0,n-m+1):
			s = sequence[x:x+m]
			fe = self.calculate(s)
			if(fe < threshold):
				startPos = x+1#1 indexed
				orientation = 1
				endPos = x+m#returns final pos in sequence
				yield(startPos,endPos,fe,orientation)
			if(both_strands):
				score = rc.calculate(s)
				if (score < threshold):
					orientation = -1
					startPos = x+m #1 indexed  based on main strand
					endPos = x+1 #returns final pos in sequence
					yield(startPos,endPos,fe,orientation)





	"""Helper Methods"""
	#loops through array and finds all final free energies below a certain threshold.
	@staticmethod
	def cutOff(array, threshold):
		returnArray = []
		for x in range(0,len(array)):
			if(array[x] < threshold):
				returnArray.append(x)
		return returnArray

	def maxE(self):
		y = 0
		s = ""
		for x in range(0,len(self.psfem['C'])):
			maxVal = self.psfem["C"][x]
			maxIndex = "C"
			if(self.psfem["C"][x] > maxVal):
				maxVal = self.psfem["C"][x]
				maxIndex = "C"
			if(self.psfem["T"][x] > maxVal):
				maxVal = self.psfem["T"][x]
				maxIndex = "T"
			if(self.psfem["A"][x] > maxVal):
				maxVal = self.psfem["A"][x]
				maxIndex = "A"
			if(self.psfem["G"][x] > maxVal):
				maxVal = self.psfem["G"][x]
				maxIndex = "G"

			y = y + maxVal
			s = s + maxIndex
		y = y + self.consensus_fe
		return [y,s]

	def minE(self):
		y = 0
		s = ""
		for x in range(0,len(self.psfem['C'])):
			minVal = self.psfem["C"][x]
			minIndex = "C"
			if(self.psfem["C"][x] < minVal):
				minVal = self.psfem["C"][x]
				minIndex = "C"
			if(self.psfem["T"][x] < minVal):
				minVal = self.psfem["T"][x]
				minIndex = "T"
			if(self.psfem["A"][x] < minVal):
				minVal = self.psfem["A"][x]
				minIndex = "A"
			if(self.psfem["G"][x] < minVal):
				minVal = self.psfem["G"][x]
				minIndex = "G"

			y = y + minVal
			s = s + minIndex
		y = y + self.consensus_fe
		return [y,s]

"""read and get based off of pssm class functions"""
def read_in_psfem_file(filename):
	psfeMatrix = {}
	with open(filename) as fin:
		"""For info"""
		words =[]
		for line in fin:
			words.append(line.split("\t"))

		for x in range(0,len(words)):
			for y in range(0,len(words[x])):
				if('\n' in words[x][y]):
					words[x][y] = words[x][y][0:len(words[x][y])-1]
		sequence = words[0][1]
		s = words[1][1]
		baseKd = float(s)

		recommended_threshold = float(words[2][1])

		"""For values"""
		table = words[3:len(words)]
		for x in range(0,4):
			c = str(table[0][x])
			psfeMatrix.update({c:[]})
			for y in range(1,len(table)):
				psfeMatrix[str(table[0][x])].append(float(table[y][x]))
	return psfem(sequence,baseKd,psfeMatrix,recommended_threshold)

def get_psfem(tf_name,directory = "./PWM/"):
	file_list = os.listdir(directory)
	pattern = tf_name
	for file_name in file_list:
		if re.search(pattern,file_name,re.I):#case insensitive
			if file_name[0] != ".": #make sure not editing ._ files
				return read_in_psfem_file(directory+file_name)



"""Added from analyze_sig_hits_v2 and edited to fit psfem"""

def add_seq_to_frame(sighits_frame,yeast_genomic_sequence_file='./S288C_reference_sequence_R61-1-1_20080605.fsa',bases_into_gene=150):
	"""This function takes a significant hits files as input and returns a frame with the significant hits
	file as a frame plus the sequence of the intergenic region as an added column.  The sequence is a 
	SeqRecord object"""
	
	#open significant hits file and get top (or bottom) r entries.  
	

	#load the yeast genome into memory
	yeast_genome = {}
	count = 1
	for seq_record in SeqIO.parse(yeast_genomic_sequence_file,"fasta"):
		yeast_genome[count] = seq_record
		yeast_genome[count].alphabet = IUPAC.unambiguous_dna
		count = count + 1

	seq_frame = pd.DataFrame(index = sighits_frame.index,columns = ['Sequence'])

	#loop through the intergenic region names and put the corresponding sequences in a list of sequence records
	#name = name + chormosome and start and end.  Sequence = grabbed sequence

	for idx,row in sighits_frame.iterrows():
		ig_sequence = yeast_genome[int(sighits_frame.ix[idx,'Chr'])].seq[int(sighits_frame.ix[idx,'Start'])-1-bases_into_gene:int(sighits_frame.ix[idx,'End'])+bases_into_gene]
		ig_sequence.alphabet = IUPAC.unambiguous_dna
		seq_frame.ix[idx]["Sequence"] = ig_sequence
	sighits_frame = sighits_frame.join(seq_frame)
	return sighits_frame

#Nik - Made edit to function to allow for harbison p value.  For harbison to work, poisson must be false.
def make_comparison_frame(file1,file2,name1 = 'TF1',name2 = 'TF2',poisson=True,bs = True,harbison = False):
	frame1 = pd.read_csv(file1, delimiter = "\t",index_col=0)
	frame2 = pd.read_csv(file2, delimiter = "\t",index_col=0)
	frame1seq = add_seq_to_frame(frame1)
	if bs == True:
		comparison_frame = pd.DataFrame(index = frame1seq.index,columns = [name1+' TPH BS',name1+' pvalue',name2+' TPH BS',name2+' pvalue','Sequence','Log2FC','Left Common Name','Right Common Name'])
		comparison_frame[name1+' TPH BS'] = frame1seq['TPH BS']
		comparison_frame[name2+' TPH BS'] = frame2['TPH BS']
	else:
		comparison_frame = pd.DataFrame(index = frame1seq.index,columns = [name1+' TPH',name1+' pvalue',name2+' TPH',name2+' pvalue','Sequence','Log2FC','Left Common Name','Right Common Name'])
		comparison_frame[name1+' TPH'] = frame1seq['TPH']
		comparison_frame[name2+' TPH'] = frame2['TPH']
	comparison_frame['Sequence'] = frame1seq['Sequence']
	comparison_frame['Left Common Name'] = frame1seq['Left Common Name']
	comparison_frame['Right Common Name'] = frame1seq['Right Common Name']
	if poisson:
		comparison_frame[name1+' pvalue'] = frame1seq['Poisson pvalue']
		comparison_frame[name2+' pvalue'] = frame2['Poisson pvalue']
	else:
		if harbison:
			comparison_frame[name1+' pvalue'] = frame1seq['Harbison IG pvalue']
			comparison_frame[name2+' pvalue'] = frame2['Harbison IG pvalue']
		else:
			comparison_frame[name1+' pvalue'] = frame1seq['CHG pvalue']
			comparison_frame[name2+' pvalue'] = frame2['CHG pvalue']
	if bs == True:
		temp1 = frame1seq['TPH BS']
		temp2 = frame2['TPH BS']
	else:
		temp1 = frame1seq['TPH']
		temp2 = frame2['TPH']
	#set TPH values of zero to 1 so I can take the log
	temp1.loc[temp1 == 0] = 1
	temp2.loc[temp2 == 0] = 1
	comparison_frame['Log2FC'] = np.log2(temp2/temp1)
	return comparison_frame

#making an empty list
#df['empty_list'] = np.empty((len(df), 0)).tolist()

def add_psfem_info(inframe,tfname,cutoff,remove_overlaps=True):
	#This function takes a significant hits frame with sequence information add 
	#adds annotation about the presence of binding sites for a tf
	#it adds 3 columns.  One with a list of coordinates for each hit
	#one with the scores, and one with the orientations

	hit_starts_megalist = []
	hit_score_megalist = []
	hit_orientation_megalist = []
	#get the psfem
	tf_psfem = get_psfem(tfname)

	for idx,row in inframe.iterrows():
		hit_starts_list = []
		hit_coords_list = []
		hit_orientation_list = []
		hit_score_list = []
		#loop through all hits in a sequence
		for startpos,endpos,score,orientation in tf_psfem.search(row['Sequence'],cutoff): #test to remove warning
			#loop through all previous hits and check for overlap
			#if remove_overlaps: #Reason remove_overlaps = False doesn't work
				overlap_flag = 0
				for ind2,val in enumerate(hit_coords_list):
					if range(max(min(startpos,endpos),min(val[0],val[1])),min(max(startpos,endpos),max(val[0],val[1]))+1):
						#replace if overlap and a better hit
						if score > hit_score_list[ind2]:
							hit_starts_list[ind2] = startpos
							hit_coords_list[ind2] = [startpos,endpos]
							hit_orientation_list[ind2] = orientation
							hit_score_list[ind2] = score
							overlap_flag = 1
						#do nothing if not a better hit
						else:
							overlap_flag = 1
				#if no overlap add to list of hits
				if not overlap_flag:
					hit_starts_list.append(startpos)
					hit_coords_list.append([startpos,endpos])
					hit_orientation_list.append(orientation)
					hit_score_list.append(score)

		#put hits in megalist
		hit_starts_megalist.append(hit_starts_list)
		hit_score_megalist.append(hit_score_list)
		hit_orientation_megalist.append(hit_orientation_list)

	inframe[tfname+' pos PSFEM'] = hit_starts_megalist
	inframe[tfname+' score PSFEM'] = hit_score_megalist
	inframe[tfname+' orientaion PSFEM'] = hit_orientation_megalist

	return inframe		

def getScores(inframe,tfname,cutoff):
	outFrame = pd.DataFrame(index = inframe.index,columns = ['Left Common Name','Right Common Name'])
	scores = []
	tf_psfem = get_psfem(tfname)
	for idx,row in inframe.iterrows():
		scoreMin = getSingleScore(row['Sequence'],tfname)
		scores.append(scoreMin)
	outFrame['Left Common Name'] = inframe['Left Common Name']
	outFrame['Right Common Name'] = inframe['Right Common Name']
	outFrame['Score ' + tfname] = scores
	return outFrame

def getSingleScore(sequence,tf):
	tf_psfem = get_psfem(tf)
	scores = [x[2] for x in tf_psfem.search(sequence,tf_psfem.max_fe)]
	if scores:
		y = min(scores)
	else:
		y = tf_psfem.max_fe
	return y

def getDoubleScore(sequence,tf,tf2):
	tf_psfem = get_psfem(tf)
	tf_psfem2 = get_psfem(tf2)
	scores = [x[2] for x in tf_psfem.search(sequence,tf_psfem.max_fe)]
	scores2 = [x[2] for x in tf_psfem2.search(sequence,tf_psfem2.max_fe)]
	if scores:
		y = min(scores)
	else:
		y = tf_psfem.max_fe
	if scores2:
		if(y > min(scores2)):
			y = min(scores2)
	else:
		if(y > tf_psfem2.max_fe):
			y = tf_psfem2.max_fe
	return y


#From analyze_sig_hits
def plotROCCurve(tpr,fpr, color = 'r',label = 'ROC ', mode = "", randomline=True):

    plt.figure(figsize=(6, 6), dpi=100)
    plt.xlabel("FPR", fontsize=16)
    plt.ylabel("TPR", fontsize=16)
    plt.title("ROC Curve " + str(mode), fontsize=16)
    plt.plot(fpr, tpr, color=color, linewidth=2, label=label)

    if randomline:
        x = [0.0, 1.0]
        plt.plot(x, x, linestyle='dashed', color='red', linewidth=2, label='random')

    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.tick_params(axis='both', which='major', labelsize=16)
	#plt.tick_params(axis='both', which='minor', labelsize=8)
    plt.legend(fontsize=14, loc='best')
    plt.tight_layout()
    plt.show()



def get_tpr_fpr_lists(positives,negatives,scoring_function,lowerbetter = True):
	flist = [positives, negatives]
	total = pd.concat(flist)
	total["Score"] = total["Sequence"].apply(scoring_function)
	total = total.sort_values(by=['Score'],ascending = lowerbetter)
	foundpositives = 0
	foundnegatives = 0
	num_positives = len(positives)
	num_negatives = len(negatives)
	tpr_list = []
	fpr_list = []
	for idx,row in total.iterrows():
		if idx in positives.index:
			foundpositives += 1.0
		else:
			foundnegatives += 1.0
		tpr_list.append(foundpositives/float(num_positives))
		fpr_list.append(foundnegatives/float(num_negatives))
	return tpr_list,fpr_list
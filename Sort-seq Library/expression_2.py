#!/usr/bin/env python 

import os
import sys
import pandas as pd
from Bio import SeqIO
import csv



reference_file = sys.argv[1]
r1_file = sys.argv[2]
r2_file = sys.argv[3]
output_file = sys.argv[4]


def read_bc_file(ref_file):
	reader=csv.reader(open(ref_file,"r"),delimiter = ',')
	d = {}
	for row in reader:
		seq, bc= row
		d[bc]={"seq":0}
	return d


def ReverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def read_reference_file(ref_file):
	reader=csv.reader(open(ref_file,"r"),delimiter = ',')
	d = {}
	for row in reader:
		seq, bc= row
		d[bc]=seq
	return d


def parse_paired_reads_file(ref_filename,read1_file,read2_file):



	barcode_dict=read_bc_file(ref_filename)
	
	handle1 = open(read1_file,"rU")
	handle2 = open(read2_file,"rU")
	read2_iter = SeqIO.parse(handle2,"fastq")
	
	ref_file = read_reference_file(ref_filename)


	count =0
	for record1 in SeqIO.parse(handle1,"fastq"): 
		record2 = next(read2_iter)
		

		full = (str(record1.seq)+ReverseComplement(record2.seq)).upper()
		#print(full)
		
		count+=1
		if count%10000==1:
			print (count)
		
	

		
		for key in barcode_dict.keys():
		
			if "CTATAGGGCGAATTGGGTA" in full[0:20]:
				
				if key in full[180:]:
					barcode_dict[key]["seq"]+=1
					# match_count=0
					# if ref_file[key][20:26] in full[18:28]:
					# 	match_count+=1
					# if ref_file[key][28:32] in full[26:34]:
					# 	match_count+=1
					# if ref_file[key][34:40] in full[32:42]:
					# 	match_count+=1
					# if ref_file[key][44:50] in full[42:52]:
					# 	match_count+=1
					# if ref_file[key][54:60] in full[52:62]:
					# 	match_count+=1
					# if ref_file[key][64:70] in full[62:72]:
					# 	match_count+=1
					# if ref_file[key][74:80] in full[72:82]:
					# 	match_count+=1
					# if ref_file[key][84:90] in full[82:92]: 
					# 	match_count+=1
					# if ref_file[key][94:100] in full[92:102]:
					# 	match_count+=1

					# if ref_file[key][104:110] in full[102:112]: 
					# 	match_count+=1
					# if ref_file[key][114:120] in full[112:122]: 
					# 	match_count+=1
					# if ref_file[key][124:130] in full[122:132]: 
					# 	match_count+=1
					# if ref_file[key][134:140] in full[132:144]:
					# 	match_count+=1
					# if ref_file[key][154:160] in full[152:164]:
					# 	match_count+=1
					# if ref_file[key][164:170] in full[162:174]:
					# 	match_count+=1
					# if ref_file[key][174:180] in full[172:183]:
					# 	match_count+=1
					# if ref_file[key][184:190] in full[182:192]:
					# 	match_count+=1
					
	

					# if match_count>14:
					# 	barcode_dict[key]["match"]+=1
					# else:
					# 	#print (full)
					# 	barcode_dict[key]["wrong"].append(full[0:150]+full[220:300])
				
				


	

	return barcode_dict


ref = open (reference_file,"r")
table=[]

for row in csv.reader(ref):
	table.append(row)
ref.close()

datadict = parse_paired_reads_file(reference_file, r1_file, r2_file)









new=[]
for row in table:
	barcode = row[1]

	if barcode=="barcode":
		continue

	seq = datadict.get(barcode)['seq']


	new.append([row[0],barcode,seq])





with open(output_file,"wb") as csv_file:
	writer=csv.writer(csv_file)
	for row in new:
		writer.writerow(row)










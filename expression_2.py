#!/usr/bin/env python

import os
import sys
import pandas as pd
from Bio import SeqIO
import csv

#code written by JL 


reference_file = sys.argv[1]
r1_file = sys.argv[2]
r2_file = sys.argv[3]
output_file = sys.argv[4]


def read_bc_file(ref_file):
	""" Function reads reference_file which contains a column with every lib seq and a second column with the barcode"""
	reader=csv.reader(open(ref_file,"r"),delimiter = ',')
	d = {}
	for row in reader:
		seq, bc= row
		d[bc]={"seq":0}
	return d


def ReverseComplement(seq):
	"""Function returns the sequence reverse complement """
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
	"""This function determines the read count for each sequence"""



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

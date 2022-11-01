"""
This module is designed to display ROC curves with one or two lines for tf comparison or any ROC curve with plotROC.

NJM
5/29/18
"""

from KdtoEv2 import psfem
from KdtoEv2 import get_psfem
from KdtoEv2 import make_comparison_frame,add_psfem_info, getScores,get_tpr_fpr_lists, plotROCCurve, getSingleScore, getDoubleScore
from analyze_sig_hits_v2 import pssm, add_pwm_info,make_comparison_frame, get_pssm, get_max_pwm_score
from math import log,e
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
import os

hits = []

"""Used to make ROC curve of one or two TFs"""
def disp_tf_ROC(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4', tf2 ='', PSFEM = False):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	if(PSFEM):
		tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getSingleScore(x,str(tf)))
		if tf2 != "":
			tpr2,fpr2 = get_tpr_fpr_lists(positives,negatives,lambda x: getSingleScore(x,str(tf2)))
			auc2 = metrics.auc(fpr2,tpr2)
		else:
			tpr2 = 0
			fpr2 = 0
			auc2 = 0
	else:
		tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: get_max_pwm_score(x,str(tf)),False) #Lower better false bc pwms
		if tf2 != "":
			tpr2,fpr2 = get_tpr_fpr_lists(positives,negatives,lambda x: get_max_pwm_score(x,str(tf2)),False)
			auc2 = metrics.auc(fpr2,tpr2)
		else:
			tpr2 = 0
			fpr2 = 0
			auc2 = 0
	auc = metrics.auc(fpr,tpr)
	plotROC(tpr,fpr, tf, auc, tf2, tpr2, fpr2, auc2)



def calc_tf_auc(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4', PSFEM = False, halfSite = False):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD', bs = False)
	#cf.to_csv('testCF.txt', sep = '\t')
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	if(PSFEM):
		tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getSingleScore(x,str(tf)))
	elif(halfSite):
		tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: get_sum_pwm_score(x,str(tf)),False) #Lower better false bc pwms
	else:
		tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: get_max_pwm_score(x,str(tf)),False) #Lower better false bc pwms
	auc = metrics.auc(fpr,tpr)
	return auc


def get_sum_pwm_score(sequence,tf):  #i changed the order, so flip the lambda.
	tf_pssm = get_pssm(tf)
	vector = [x[2] for x in tf_pssm.search(sequence,tf_pssm.min())]
	if vector:
		y = sum(vector)
	else:
		y = tf_pssm.min()
	return y



"""Utility function to plot any ROC curve w up to two lines"""
def plotROC(tpr,fpr, name, auc, name2 = "", tpr2 = 0,fpr2 = 0, auc2 = 0, color = 'r',label = 'ROC',randomline=True):
	"""From analyze sig hits"""
	plt.figure(figsize=(6, 6), dpi=100)
	plt.xlabel("FPR", fontsize=16)
	plt.ylabel("TPR", fontsize=16)
	plt.title("ROC Curve", fontsize=16)
	plt.plot(fpr, tpr, color='navy', linewidth=2, label= str(name) + ' AUC = ' + str(auc))
	if name2 != "":
		plt.plot(fpr2, tpr2, color='darkorange', linewidth=2, label= str(name2) + ' AUC = ' + str(auc2))
	x = [0.0, 1.0]
	plt.plot(x, x, linestyle='dashed', color='red', linewidth=2, label='random')

	plt.xlim(0.0, 1.0)
	plt.ylim(0.0, 1.0)
	plt.tick_params(axis='both', which='major', labelsize=16)
	#plt.tick_params(axis='both', which='minor', labelsize=8)
	plt.legend(fontsize=14, loc='best')
	plt.tight_layout()
	plt.show()

def longList(inputList = []):#Leave blank for all
	if len(inputList) == 0:
		path = './recommended/PWM/'
		for filename in os.listdir(path):
			if(filename[0] != '.'): #avoids hidden files
				n = filename.find('.')
				x = filename[n+1:] #removes first name to make sure getPSSM doesn't break
				inputList.append(str(x))
		print inputList
	finalFrame = pd.DataFrame(columns = ["TF", "AUC"])
	for x in inputList:
		print x
		auc = calc_tf_auc(.0001,.1,x)
		itterFrame = pd.DataFrame([[x,auc]], columns = ["TF", "AUC"])
		finalFrame = finalFrame.append(itterFrame)
	return finalFrame





def HS_ROC(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4'):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	#cf = removeCannon(cf)
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: HScount(x),False)
	auc = metrics.auc(fpr,tpr)
	return auc




def removeCannon(df, cutoff = 13.01, tf = 'Gal4'):
	df["Score"] = df["Sequence"].apply(lambda x: get_max_pwm_score(x,tf))
	newdf = df[df["Score"]<cutoff]
	return newdf


def HS_ROCmaxN(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4', n = 200):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getMaxN(x, n),False)
	auc = metrics.auc(fpr,tpr)
	return auc

def HScount(inSeq, perKB = False):
	dna_cursor = Seq(str(inSeq), IUPAC.unambiguous_dna)
	HS = dna_cursor.count("CGG")
	HS = HS + (dna_cursor.reverse_complement()).count("CGG")
	if perKB:
		HS = HS/(len(dna_cursor)/float(1000))
	return HS

def getMaxN(inSeq, n = 200):
	seq = Seq(str(inSeq), IUPAC.unambiguous_dna)
	totals = []
	for x in range(0,len(seq)-n-1):
		inputSeq = seq[x:x+n]
		totals.append(inputSeq.count("CGG"))
	seq = seq.reverse_complement()
	for x in range(0,len(seq)-n-1):
		inputSeq = seq[x:x+n]
		totals.append(inputSeq.count("CGG"))
	if totals:
		return max(totals)
	else:
		return 0



def getSpacing(seq, upper = 15, lower = 0):
	n = 0
	while(seq.find("CGG")>=0):
		i = seq.find("CGG")
		#print seq
		seq = seq[i+3:]
		#print seq
		l = seq.find("CCG")
		if(l >= lower and l <= upper):
			n = n+1
	#print n
	return n

def getSpacingNik(seq, upper = 15, lower = 0):
	n = 0
	while len(seq) >= 6:
		i = seq.find("CGG")
		j = seq.find("CCG")
		if(i != -1 and j != -1):
			if j-i <= (upper+3) and j-i >= (lower+3):
				n = n+1
			seq = seq[i+3:] #removes CGG
			if(seq.find("CGG") != -1):
				seq = seq[seq.find("CGG"):] #move to next CGG
			else:
				return n
		else:
			return n
	return n

def getSpacingAndUp(seq, lower = 0):
	n = 0
	while(seq.find("CGG")>=0):
		i = seq.find("CGG")
		seq = seq[i+3:]
		l = seq.find("CCG")
		if(l >= lower):
			n = n+1
	return n

def getSpacingAndUpNik(seq, lower = 0):
	n = 0
	while len(seq) >= 6:
		i = seq.find("CGG")
		j = seq.find("CCG")
		if(i != -1 and j != -1):
			if j-i >= (lower+3):
				n = n+1
			seq = seq[i+3:] #removes CGG
			if(seq.find("CGG") != -1):
				seq = seq[seq.find("CGG"):] #move to next CGG
			else:
				return n
		else:
			return n
	return n
def spacingROCandUp(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4', lower = 15):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	#cf = removeCannon(cf)
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getSpacingAndUp(x,lower),False)
	auc = metrics.auc(fpr,tpr)
	return auc

def spacingROC(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4', upper = 15, lower = 0):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	#cf = removeCannon(cf)
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getSpacing(x,upper,lower),False)
	auc = metrics.auc(fpr,tpr)
	return auc



def spacingMax400(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4', upper = 15, lower = 0):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	cf = removeCannon(cf)
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getSpacingMax400(x,upper, lower),False)
	auc = metrics.auc(fpr,tpr)
	return auc

def getSpacingMax400(inSeq, upper = 15, lower = 0):
	seq = Seq(str(inSeq), IUPAC.unambiguous_dna)
	totals = []
	for x in range(0,len(seq)-399):
		inputSeq = seq[x:x+400]
		totals.append(getSpacing(inputSeq, upper, lower))
	if totals:
		if max(totals)>0:
			hits.append(max(totals))
		return max(totals)
	else:
		return 0


def spacingMax400AndUp(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4', lower = 0):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	cf = removeCannon(cf)
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getSpacingMax400AndUp(x,lower),False)
	auc = metrics.auc(fpr,tpr)
	return auc
def getSpacingMax400AndUp(inSeq, upper = 15, lower = 0):
	seq = Seq(str(inSeq), IUPAC.unambiguous_dna)
	totals = []
	for x in range(0,len(seq)-399):
		inputSeq = seq[x:x+400]
		totals.append(getSpacingAndUp(inputSeq, lower))
	if totals:
		return max(totals)
	else:
		return 0




def spacingMax400Double(lowerCutoff = .00001, upperCutoff = .1, tf = 'Gal4'):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getSpacingMax400Double(x),False)
	auc = metrics.auc(fpr,tpr)
	return auc

def getSpacingMax400Double(inSeq):
	seq = Seq(str(inSeq), IUPAC.unambiguous_dna)
	totals = []
	for x in range(0,len(seq)-399):
		inputSeq = seq[x:x+400]
		totals.append(getSpacing(inputSeq, 19, 18)+getSpacing(inputSeq,11,10))
	if totals:
		if max(totals)>0:
			hits.append(max(totals))
		return max(totals)
	else:
		return 0





def PWMMax400(lowerCutoff = .00001, upperCutoff = .1, tf = 'HalfSite'):
	cf = make_comparison_frame('sig_prom_Gal4_FL_HarbIG_Final.txt','sig_prom_Gal4_DBD_HarbIG_Final.txt','FL','DBD')
	positives = cf[cf["FL pvalue"]<lowerCutoff]
	negatives = cf[cf["FL pvalue"]>upperCutoff]
	tpr,fpr = get_tpr_fpr_lists(positives,negatives,lambda x: getPWMMax400(x, tf),False)
	auc = metrics.auc(fpr,tpr)
	return auc

def getPWMMax400(inSeq,tf):
	seq = Seq(str(inSeq), IUPAC.unambiguous_dna)
	totals = []
	for x in range(0,len(seq)-399):
		inputSeq = seq[x:x+400]
		totals.append(get_sum_pwm_score(inputSeq, tf))
	if totals:
		if max(totals)>0:
			hits.append(max(totals))
		return max(totals)
	else:
		return 0



"""
print str(spacingMax400(upper = 19,lower = 18))
print sum(hits)
hits1 = hits
hits = []
print str(spacingMax400(upper = 19,lower = 19))
print sum(hits)
hits2 = hits
hits = []
print str(spacingMax400(upper = 18,lower = 18))
print sum(hits)

print hits1
print hits2
print hits
"""
#print calc_tf_auc()
#print PWMMax400(tf = 'HalfSite')
#print PWMMax400(tf = 'Gal4')
#print calc_tf_auc(tf = 'Gal4')
"""
print str(spacingMax400(upper = 6,lower = 2))
print str(spacingMax400(upper = 7,lower = 3))
print str(spacingMax400(upper = 8,lower = 4))
print str(spacingMax400(upper = 9,lower = 5))
print str(spacingMax400(upper = 10,lower = 6))
print str(spacingMax400(upper = 11,lower = 7))
print str(spacingMax400(upper = 12,lower = 8))
print str(spacingMax400(upper = 13,lower = 9))
print str(spacingMax400(upper = 14,lower = 10))
print str(spacingMax400(upper = 15,lower = 11))
print str(spacingMax400(upper = 16,lower = 12))
print str(spacingMax400(upper = 17,lower = 13))
print str(spacingMax400(upper = 18,lower = 14))
print str(spacingMax400(upper = 19,lower = 15))
print str(spacingMax400(upper = 20,lower = 16))
print str(spacingMax400(upper = 21,lower = 17))
print str(spacingMax400(upper = 22,lower = 18))
print str(spacingMax400(upper = 23,lower = 19))
print str(spacingMax400(upper = 24,lower = 20))
print str(spacingMax400(upper = 25,lower = 21))
"""

print str(spacingROC(upper = 11, lower = 10))

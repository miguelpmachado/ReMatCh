import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys
import csv
import numpy
import shlex, subprocess,ftplib
import os.path

def createCheckFile(bamSortedPath, sequenceMedObject):

	check_fileName = bamSortedPath.replace('_sorted', '')

	with open(check_fileName + "_mappingCheck.tab", 'w') as mappingCheckFile:
		mappingCheckFile.write('#Sequence\tDuplication\tIndel\tRawCoverage\tAlternativeQualityScore\tCoverage\tMultipleAllele\tmeanRawCoverageFragmentMapped\tmeanRawCoverageRegionInterest\n')
		for sequence in sequenceMedObject:
			mappingCheckFile.write(sequence + '\t' + str(round(sequenceMedObject[sequence][8], 3)) + '\t' + str(round(sequenceMedObject[sequence][9],3)) + '\t' + str(round(sequenceMedObject[sequence][10],3)) + '\t' +str(sequenceMedObject[sequence][4]) + '\t' +str(sequenceMedObject[sequence][5]) + '\t' +str(sequenceMedObject[sequence][6]) + '\t' + str(round(sequenceMedObject[sequence][1], 3)) + '\t' + str(round(sequenceMedObject[sequence][11], 3)) + '\n')


def filter_vcf(pathToVcf, extraSeq, sequenceMedObject):

	with open(pathToVcf, 'r') as vcfFile:
		with open(pathToVcf + '.temp', 'w') as vcfTemp:

			for line in csv.reader(vcfFile, delimiter="\t"):
				if not line[0].startswith("#"):
					if int(line[1]) > extraSeq and int(line[1]) <= len(sequenceMedObject[line[0]][3]) - extraSeq:
						vcfTemp.write('\t'.join([str(x) for x in line]) + '\n')
				else:
					vcfTemp.write('\t'.join([str(x) for x in line]) + '\n')

	os.remove(pathToVcf)
	os.rename(pathToVcf + '.temp', pathToVcf)



def changeFastaHeadersAndTrimm(FastasequencesFile,TrimmExtraSeq,referencePath):
	
	headersArray=[]
	
	with open(referencePath, 'r') as seqFileRef:
		for line in seqFileRef:
			if '>' in line:
				headersArray.append(line)
			
	with open(FastasequencesFile, 'r') as seqFile:
		with open(FastasequencesFile+".temp", 'w') as tempFile:
			tempStr=''
			headercounter=0
			for line in seqFile:
				
				if '>' in line:
					
					if TrimmExtraSeq!=0 and len(tempStr)>0:
						tempStr=tempStr[TrimmExtraSeq:len(tempStr)-TrimmExtraSeq]
					if headercounter>0:
						tempFile.write(tempStr+"\n")
						
					tempStr=''
					
					tempFile.write(headersArray[headercounter])
					headercounter+=1
				
				else:
					tempStr+=line.replace('\n', '').replace('\r', '')
				
			if TrimmExtraSeq!=0 and len(tempStr)>0:
				tempStr=tempStr[TrimmExtraSeq:len(tempStr)-TrimmExtraSeq]	
			tempFile.write(tempStr)
			
	os.remove(FastasequencesFile)
	os.rename(FastasequencesFile+".temp", FastasequencesFile)

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

class Logger(object):
	def __init__(self, out_directory):
		self.logfile = os.path.join(out_directory, "run.log")
		if os.path.isfile(self.logfile):
			print "Logfile already exists! It will be overwritten..." + "\n"
		self.terminal = sys.stdout
		self.log = open(self.logfile, "w")
	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)
		self.log.flush()
	def flush(self):
		pass



def createCheckFile(bamSortedPath, sequenceMedObject):

	check_fileName = bamSortedPath.replace('_sorted', '')

	with open(check_fileName + "_mappingCheck.tab", 'w') as mappingCheckFile:
		mappingCheckFile.write('#Sequence\tDuplication\tIndel\tCoverage\tLowSNPsQualityScore\tSNPCoverage\tSNPMultipleAllele\tmeanSequenceCoverage\tmeanTrimmedSequenceCoverage\n')
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

def removeFromArray(toClear):
	for i in toClear:
		os.system('rm ' + i)

def checkPrograms():

	print 'Checking dependencies...'
	which_program = ['which', '']
	programs = {'bedtools':['>=','2.22'], 'java':['>=', '1.8'], 'samtools':['=', '1.2'], 'bcftools':['=', '1.2'],'bowtie2':['>=','2.2.6']}

	for program in programs:
		which_program[1] = program
		proc = subprocess.Popen(which_program, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout,stderr = proc.communicate()
		if proc.returncode != 0:
			sys.exit(program + ' not found in PATH.')
		else:
			if program =='java':
				check_version = [stdout.strip('\n'), '-version']
			else:
				check_version = [stdout.strip('\n'), '--version']

			print program + ' found at: ' + check_version[0]
			print check_version
			proc = subprocess.Popen(check_version, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdout,stderr = proc.communicate() 
			print stdout
			version_line=stdout.split('\n')[0].split(' ')[-1]
			print version_line
			version_line=version_line.replace('"', '')
			print version_line
			if 'v' in version_line:
				version_line=version_line.split('v')[1]
			elif 'V' in version_line:
				version_line=version_line.split('V')[1]
			if programs[program][0] == '>=':
				if float('.'.join(version_line.split('.')[1:2])) < float('.'.join(programs[program][1].split('.')[1:2])):
					sys.exit('ReMatCh requires ' + program + ' with version ' + programs[program][1] + ' or above.')
			else:
				if version_line != programs[program][1]:
					sys.exit('ReMatCh requires ' + program + ' with version ' + programs[program][1] + '.')










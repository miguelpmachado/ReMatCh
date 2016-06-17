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
import time

class Logger(object):
	def __init__(self, out_directory):
		time_str = time.strftime("%Y%m%d-%H%M%S")
		self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
		if os.path.isfile(self.logfile):
			print '\n' + 'Logfile already exists! It will be overwritten...'
		if not os.path.isdir(out_directory):
			os.makedirs(out_directory)
		self.terminal = sys.stdout
		self.log = open(self.logfile, "w")
	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)
		self.log.flush()
	def flush(self):
		pass
	def getLogFile(self):
		return self.logfile

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

def checkPrograms(programs_version_dictionary):
	print '\n' + 'Checking dependencies...'
	programs = programs_version_dictionary
	which_program = ['which', '']
	listMissings = []
	for program in programs:
		which_program[1] = program
		proc = subprocess.Popen(which_program, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout,stderr = proc.communicate()
		if proc.returncode != 0:
			listMissings.append(program + ' not found in PATH.')
		else:
			if programs[program][0] == None:
				print program + ' (impossible to determine programme version) found at: ' + stdout.splitlines()[0]
			else:
				check_version = [stdout.splitlines()[0], programs[program][0]]
				proc = subprocess.Popen(check_version, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				stdout,stderr = proc.communicate()
				if stdout == '':
					stdout = stderr
				version_line=stdout.splitlines()[0].split(' ')[-1]
				replace_characters = ['"', 'v', 'V', '+']
				for i in replace_characters:
					version_line=version_line.replace(i, '')
				print program + ' (' + version_line + ') found at: ' + check_version[0]
				if programs[program][1] == '>=':
					program_found_version = version_line.split('.')
					program_version_required = programs[program][2].split('.')
					if float('.'.join(program_found_version[0:2])) < float('.'.join(program_version_required[0:2])):
						listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
					elif float('.'.join(program_found_version[0:2])) == float('.'.join(program_version_required[0:2])):
						if len(program_version_required) == 3:
							if len(program_found_version) == 2:
								program_found_version.append(0)
							if program_found_version[2] < program_version_required[2]:
								listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
				else:
					if version_line != programs[program][2]:
						listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
	if len(listMissings) > 0:
		sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(listMissings))


def removeIndexes(referencePath):

	referenceFileName, extension = os.path.splitext(referencePath)

	os.remove(referenceFileName+"_picard_out.txt")
	os.remove(referenceFileName + ".dict")
	if os.path.isfile(referencePath+'.fai'):
		os.remove(referencePath+'.fai')
	os.system('rm ' + referenceFileName + ".*.bt2")
	os.remove(referenceFileName+"_bowtiBuildLog.txt")

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

from remach_utils import changeFastaHeadersAndTrimm
from remach_utils import filter_vcf
from remach_utils import createCheckFile


def convertToBAM(samPath):

	filename, samfile_extension = os.path.splitext(samPath)

	os.system("samtools view -buh -o " + filename +'_temp.bam' + " " + samPath)
	os.system("rm " + samPath)
	os.system("samtools sort -o " + filename +'.bam' + " " + filename +'_temp.bam')
	os.system("rm "+ filename +'_temp.bam')
	os.system("samtools index " + filename +'.bam')

	return (filename +'')


def rawCoverage(bamSortedPath):
	os.system("bedtools genomecov -d -ibam " + bamSortedPath + ".bam > " + bamSortedPath+".tab")


def checkCoverage(outputPath, coverageThreshold,extraSeq, logFile):

	sequenceMedObject = {}
	#sequenceAndIndex = {}
		
	with open(outputPath+'.tab') as tsv:
		prevName = '';
		countlines = 0
		arrayOfcoverageValues = []
		arrayOfpositionValues = []
		sequenceNames = []
		countSequences = 0
		
		for line in csv.reader(tsv, delimiter="\t"):
			countlines += 1

			if prevName != line[0] and countlines != 1:
				countSequences += 1
				sequenceNames.append(prevName)
				arrayOfcoverageValues.append(int(line[2]))
				arrayOfpositionValues.append(int(line[1]))
				sequenceMedObject[prevName] = [prevName, numpy.average(arrayOfcoverageValues), numpy.std(arrayOfcoverageValues), arrayOfcoverageValues, False, False, False,arrayOfpositionValues]
				arrayOfcoverageValues = []
				arrayOfpositionValues = []
				prevName = line[0]
			else:
				if countlines != 1:
					arrayOfcoverageValues.append(int(line[2]))
					arrayOfpositionValues.append(int(line[1]))
				if countlines == 1:
					prevName = line[0]


		sequenceNames.append(prevName)
		sequenceMedObject[prevName] = [prevName, numpy.average(arrayOfcoverageValues), numpy.std(arrayOfcoverageValues), arrayOfcoverageValues, False, False, False,arrayOfpositionValues]
		countSequences += 1


	for sequence in sequenceMedObject:
		countLowCoverage = 0
		countIndel = 0
		countDuplication = 0
		sequenceLength = len(sequenceMedObject[sequence][3])

		index=0
		for coverageAtPosition in sequenceMedObject[sequence][3]:
			
			if (sequenceMedObject[sequence][7])[index] > extraSeq and (sequenceMedObject[sequence][7])[index] <= sequenceLength-extraSeq:
			
				if coverageAtPosition >= 1.25 * sequenceMedObject[sequence][1]:
					countDuplication += 1
				if coverageAtPosition < 0.1 * sequenceMedObject[sequence][1]:
					countIndel += 1
				if coverageAtPosition < int(coverageThreshold):
					countLowCoverage += 1
			index+=1
		

		sequenceMedObject[sequence].append(float(countDuplication)/sequenceLength)
		sequenceMedObject[sequence].append(float(countIndel)/float(sequenceLength))
		sequenceMedObject[sequence].append(float(countLowCoverage)/float(sequenceLength))


	return sequenceNames, sequenceMedObject



    		
def alleleCalling(bamSortedPath, referencePath, sequenceNames, gatkPath, sampleID, qualityThreshold, coverageThreshold, multipleAlleles, sequenceMedObject,extraSeq, logFile):

	ploidytempFile = bamSortedPath+'_temp_ploi.tab'
	sequencesFile = bamSortedPath + '_sequences.fasta'
	filteredsequencesFile = bamSortedPath + '_sequences_filtered.fasta'
	
	with open(ploidytempFile, 'w') as tempFile:
		tempFile.write(sampleID + '\t' + str(1))

	os.system("samtools mpileup --no-BAQ --fasta-ref " + referencePath + " --uncompressed -t DP,DPR,DV " + bamSortedPath + ".bam | bcftools call --multiallelic-caller --variants-only --samples-file " + ploidytempFile + " --output-type v --output " + bamSortedPath + ".vcf")
	
	print "running bcf"
	logFile.write("running bcf" + '\n')
	
	os.system("bcftools filter --SnpGap 3 --IndelGap 10 --include 'QUAL>=" + str(qualityThreshold) + " && FORMAT/DV>=" + str(coverageThreshold) + " && (FORMAT/DV)/(FORMAT/DP)>=" + str(multipleAlleles) + "' --output-type v --output " + bamSortedPath + "_filtered.vcf " + bamSortedPath + ".vcf")
	os.system("bcftools filter --SnpGap 3 --IndelGap 10 --include 'TYPE=\"snp\" && QUAL>=" + str(qualityThreshold) + " && FORMAT/DV>=" + str(coverageThreshold) + " && (FORMAT/DV)/(FORMAT/DP)>=" + str(multipleAlleles) + "' --output-type v --output " + bamSortedPath + "_filtered_without_indels.vcf " + bamSortedPath + ".vcf")
	

	filter_vcf(bamSortedPath + "_filtered.vcf", extraSeq, sequenceMedObject)
	filter_vcf(bamSortedPath + "_filtered_without_indels.vcf", extraSeq, sequenceMedObject)


	with open(bamSortedPath + ".vcf", 'r') as vcfFile:

		for line in csv.reader(vcfFile, delimiter="\t"):
			if not line[0].startswith("#"):
				quality = line[5]

				if int(line[1]) > extraSeq and int(line[1]) <= len(sequenceMedObject[line[0]][3]) - extraSeq:

					if quality < float(qualityThreshold):
						sequenceMedObject[line[0]][4] = True
					
					qualityCoverage = float(line[9].split(':')[2])
					
					if qualityCoverage < float(coverageThreshold):
						sequenceMedObject[line[0]][5] = True
				
					coverageByAllele = line[9].split(':')[3].split(',')
					alternativeAlleles = line[4].split(',')

					coverageByAllele = [ int(x) for x in coverageByAllele ]

					dominantSNP = coverageByAllele.index(max(coverageByAllele))
					coverageAllele = coverageByAllele[dominantSNP]

					if coverageAllele/qualityCoverage < float(multipleAlleles):
						sequenceMedObject[line[0]][6] = True



	print "running gatk"
	logFile.write("running gatk" + '\n')
	
	os.system("java -jar " + gatkPath + " -T FastaAlternateReferenceMaker -R "+ referencePath +" -o "+ filteredsequencesFile +" -V "+ bamSortedPath + "_filtered.vcf 2> "+sequencesFile+"Log_gatk.txt")
	os.system("java -jar " + gatkPath + " -T FastaAlternateReferenceMaker -R "+ referencePath +" -o "+ bamSortedPath + "_sequences_filtered_without_indels.fasta -V "+ bamSortedPath + "_filtered_without_indels.vcf 2>> "+sequencesFile+"Log_gatk.txt")
	os.system("java -jar " + gatkPath + " -T FastaAlternateReferenceMaker -R "+ referencePath +" -o "+ sequencesFile +" -V "+ bamSortedPath + ".vcf 2>> "+sequencesFile+"Log_gatk.txt")

	os.system("rm " + ploidytempFile)

	
	changeFastaHeadersAndTrimm(filteredsequencesFile, extraSeq,referencePath)
	changeFastaHeadersAndTrimm(bamSortedPath + "_sequences_filtered_without_indels.fasta", extraSeq,referencePath)
	changeFastaHeadersAndTrimm(sequencesFile, 0,referencePath)


	createCheckFile(bamSortedPath, sequenceMedObject)







        

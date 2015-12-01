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



def convertToBAM(samPath):

	filename, samfile_extension = os.path.splitext(samPath)

	os.system("samtools view -buh -o " + filename +'_temp.bam' + " " + samPath)
	os.system("rm " + samPath)
	os.system("samtools sort " + filename +'_temp.bam' + " " + filename)
	os.system("rm "+ filename +'_temp.bam')
	os.system("samtools index " + filename +'.bam')

	return (filename +'')


def rawCoverage(bamSortedPath):

	os.system("bedtools genomecov -d -ibam " + bamSortedPath + ".bam > " + bamSortedPath+".tab")


def changeFastaHeaders(FastasequencesFile,TrimmExtraSeq,sequenceAndIndex):
	with open(FastasequencesFile, 'r') as seqFile:
		with open(FastasequencesFile+".temp", 'w') as tempFile:
			tempStr=''
			for line in seqFile:
				
				if '>' in line:
					
					if TrimmExtraSeq!=0 and len(tempStr)>0:
						tempStr=tempStr[TrimmExtraSeq:len(tempStr)-TrimmExtraSeq-1]
					
					tempFile.write(tempStr)
					tempStr='\n'
					number = line.split('>')[1].strip("\n").strip("\r")
					lineToUse = '>' + sequenceAndIndex[number] + '\n'
					tempFile.write(lineToUse)
				else:
					tempStr+=line.replace('\n', '').replace('\r', '')
					#tempFile.write(line)
				
			tempFile.write(tempStr)
			
	os.remove(FastasequencesFile)
	os.rename(FastasequencesFile+".temp", FastasequencesFile)



def checkCoverage(outputPath, coverageThreshold,extraSeq):

	sequenceMedObject = {}
	sequenceAndIndex = {}
		
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
				sequenceAndIndex[str(countSequences)] = prevName
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
		sequenceAndIndex[str(countSequences)] = prevName

	#print sequenceAndIndex
	#with open(outputPath+'_coverageCheck.tab', 'w') as coverageCheckFile:

	for sequence in sequenceMedObject:
		countLowCoverage = 0
		countIndel = 0
		countDuplication = 0
		sequenceLength = len(sequenceMedObject[sequence][3])

		print "%s: %s bp; %s meanCoverage." % (sequence, sequenceLength, sequenceMedObject[sequence][1])
		
		index=0
		for coverageAtPosition in sequenceMedObject[sequence][3]:
			
			if (sequenceMedObject[sequence][7])[index]> extraSeq and (sequenceMedObject[sequence][7])[index]<= sequenceLength-extraSeq:
			
				if coverageAtPosition >= 1.25 * sequenceMedObject[sequence][1]:
					countDuplication += 1
				elif coverageAtPosition < 0.1 * sequenceMedObject[sequence][1]:
					countIndel += 1
				elif coverageAtPosition < int(coverageThreshold):
					countLowCoverage += 1
			index+=1
		
		sequenceMedObject[sequence].append(str(float(countDuplication)/sequenceLength))
		sequenceMedObject[sequence].append(str(float(countIndel)/float(sequenceLength)))
		sequenceMedObject[sequence].append(str(float(countLowCoverage)/float(sequenceLength)))
			#coverageCheckFile.write(sequence + '\t' + str(float(countDuplication)/sequenceLength) + '\t' + str(float(countcountIndel)/float(sequenceLength)) + '\t' + str(float(countLowCoverage)/float(sequenceLength))+"\n")

	return sequenceNames, sequenceMedObject, sequenceAndIndex
    		
def alleleCalling(bamSortedPath, referencePath, sequenceNames, gatkPath, sampleID, qualityThreshold, coverageThreshold, multipleAlleles, sequenceMedObject, sequenceAndIndex,threadnumb,extraSeq):

	ploidytempFile = bamSortedPath+'_temp_ploi.tab'

	sequencesFile = bamSortedPath + '_sequences.fasta'

	filteredsequencesFile = bamSortedPath + '_sequences_filtered.fasta'
	
	
	

	with open(ploidytempFile, 'w') as tempFile:
		tempFile.write(sampleID + '\t' + str(1))

	os.system("samtools mpileup --no-BAQ --fasta-ref " + referencePath + " --uncompressed -t DP,DPR,DV " + bamSortedPath + ".bam | bcftools call --multiallelic-caller --variants-only --samples-file " + ploidytempFile + " --output-type v --output " + bamSortedPath + ".vcf")
	
	#os.system("bcftools filter --include 'QUAL>=" + str(qualityThreshold) + " && FORMAT/DP>=" + str(coverageThreshold) + "' --output-type v --output " + bamSortedPath + "_filtered.vcf " + bamSortedPath + ".vcf")
	os.system("bcftools filter --SnpGap 3 --IndelGap 10 --include 'QUAL>=" + str(qualityThreshold) + " && FORMAT/DV>=" + str(coverageThreshold) + " && (FORMAT/DV)/(FORMAT/DP)>=" + str(multipleAlleles) + "' --output-type v --output " + bamSortedPath + "_filtered.vcf " + bamSortedPath + ".vcf")
	os.system("bcftools filter --SnpGap 3 --IndelGap 10 --include 'TYPE=\"snp\" && QUAL>=" + str(qualityThreshold) + " && FORMAT/DV>=" + str(coverageThreshold) + " && (FORMAT/DV)/(FORMAT/DP)>=" + str(multipleAlleles) + "' --output-type v --output " + bamSortedPath + "_filtered_without_indels.vcf " + bamSortedPath + ".vcf")

	os.system("java -jar " + gatkPath + " -T FastaAlternateReferenceMaker -R "+ referencePath +" -o "+ filteredsequencesFile +" -V "+ bamSortedPath + "_filtered.vcf")
	os.system("java -jar " + gatkPath + " -T FastaAlternateReferenceMaker -R "+ referencePath +" -o "+ bamSortedPath + "_sequences_filtered_without_indels.fasta -V "+ bamSortedPath + "_filtered.vcf")

	os.system("java -jar " + gatkPath + " -T FastaAlternateReferenceMaker -R "+ referencePath +" -o "+ sequencesFile +" -V "+ bamSortedPath + ".vcf")

	#os.system("rm " + testFile)
	os.system("rm " + ploidytempFile)

	
	changeFastaHeaders(filteredsequencesFile, extraSeq,sequenceAndIndex)
	changeFastaHeaders(bamSortedPath + "_sequences_filtered_without_indels.fasta", extraSeq,sequenceAndIndex)
	changeFastaHeaders(sequencesFile, extraSeq,sequenceAndIndex)


	with open(bamSortedPath + ".vcf", 'r') as vcfFile:
		for line in csv.reader(vcfFile, delimiter="\t"):
			if not line[0].startswith("#"):
				quality = line[5]
				if quality < float(qualityThreshold):
					sequenceMedObject[line[0]][4] = True
				
				'''if line[7].startswith("INDEL"):
					
					deepCoverage = float(line[9].split(':')[1].split('=')[1])
				
					if deepCoverage < float(coverageThreshold):
						sequenceMedObject[line[0]][5] = True'''
				#else:
				qualityCoverage = float(line[9].split(':')[2])
				
				if qualityCoverage < float(coverageThreshold):
					sequenceMedObject[line[0]][5] = True
			
				coverageByAllele = line[9].split(':')[3].split(',')
				alternativeAlleles = line[4].split(',')

				coverageByAllele = [ int(x) for x in coverageByAllele ]

				dominantSNP = coverageByAllele.index(max(coverageByAllele))
				coverageAllele = coverageByAllele[dominantSNP]
				#alternativeSNP = alternativeAlleles[dominantSNP-1] ??

				if coverageAllele/qualityCoverage < float(multipleAlleles):
					sequenceMedObject[line[0]][6] = True

	check_fileName = bamSortedPath.replace('_sorted', '')

	with open(check_fileName + "_mappingCheck.tab", 'w') as mappingCheckFile:
		mappingCheckFile.write('#Sequence\tDuplication\tIndel\tRawCoverage\tAlternativeQualityScore\tCoverage\tMultipleAllele\n')
		for sequence in sequenceMedObject:
			mappingCheckFile.write(sequence + '\t' + str(sequenceMedObject[sequence][8]) + '\t' + str(sequenceMedObject[sequence][9]) + '\t' + str(sequenceMedObject[sequence][10]) + '\t' +str(sequenceMedObject[sequence][4]) + '\t' +str(sequenceMedObject[sequence][5]) + '\t' +str(sequenceMedObject[sequence][6]) + '\n')









	#bcftools filter --SnpGap 3 --IndelGap 10 --include 'TYPE="snp" && QUAL>=10 && MIN(FORMAT/DP)>=10 && MAX(FORMAT/DP)<=920' --output-type z --output CC23_comOutgroup.filtered.gap_snp_qual10_MINformatDP10_MAXformatDP920.vcf.gz CC23_comOutgroup.bcf &&
	#bcftools index CC23_comOutgroup.filtered.gap_snp_qual10_MINformatDP10_MAXformatDP920.vcf.gz &&
	#bcftools view --output-type v CC23_comOutgroup.filtered.gap_snp_qual10_MINformatDP10_MAXformatDP920.vcf.gz | grep --invert-match '#' | wc -l

	#sampleArray = []'''



'''bcftools convert --regions 'GAS_emm81_H293_completeGenome|gki_region|1136072_1137770:601-1098' --output-type v --output gas_$subgroup.gki.vcf gas_$subgroup.bcf.gz
bcftools convert --regions 'GAS_emm81_H293_completeGenome|gtr_region|1114559_1116209:601-1050' --output-type v --output gas_$subgroup.gtr.vcf gas_$subgroup.bcf.gz
bcftools convert --regions 'GAS_emm81_H293_completeGenome|murI_region|314208_315846:601-1038' --output-type v --output gas_$subgroup.murI.vcf gas_$subgroup.bcf.gz
bcftools convert --regions 'GAS_emm81_H293_completeGenome|mutS_region|1660766_1662371:601-1005' --output-type v --output gas_$subgroup.mutS.vcf gas_$subgroup.bcf.gz
bcftools convert --regions 'GAS_emm81_H293_completeGenome|recP_region|1269857_1271516:601-1059' --output-type v --output gas_$subgroup.recP.vcf gas_$subgroup.bcf.gz
bcftools convert --regions 'GAS_emm81_H293_completeGenome|xpt_region|842025_843675:601-1050' --output-type v --output gas_$subgroup.xpt.vcf gas_$subgroup.bcf.gz
bcftools convert --regions 'GAS_emm81_H293_completeGenome|yqiL_region|134461_136095:601-1034' --output-type v --output gas_$subgroup.yqiL.vcf gas_$subgroup.bcf.gz
done &&'''


        

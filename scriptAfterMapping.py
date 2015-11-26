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

def download(dirs2,target_dir2,ref2,success2,f2,link2):
	#new folder for each reference with reference id name
	subprocess.call(['mkdir', target_dir2+"/"+ref2])
				
	#get fasta file for each read file name
	for item in dirs2:
					
		f2.cwd(link2)
		final_target_dir=target_dir2+"/"+ref2 +"/"+ item
		file = open(final_target_dir, 'wb')
		print "Downloading: %s" % item

		f2.retrbinary('RETR %s' % item, file.write)
		file.close()
		print "Downloaded %s" % item
		success2+=1		

	return success2

def download_ERR(ERR_id,target_dir):
	
	if not os.path.exists(target_dir):
		os.makedirs(target_dir)
	
	f=ftplib.FTP('ftp.sra.ebi.ac.uk')
	f.login()
	#for each reference id get the read files
	ref=ERR_id
	failed=0
	success=0

	try:
				
		firstid=ref[0:6]
		#get the read files name from the reference id
		link='/vol1/fastq/'+firstid+"/"+ref
		f.cwd(link)
		dirs=f.nlst();
				
	except:
		try:
			firstid=ref[0:6]
			#get the read files name from the reference id
			link='/vol1/fastq/'+firstid+"/00"+ref[-1]+"/"+ref
			f.cwd(link)
			dirs=f.nlst();
					
											
		except Exception, e:
			failed +=1
			print "Bad ID: " + ref
		else:
			success=download(dirs,target_dir,ref,success,f,link)	
			
	else:
				
		success=download(dirs,target_dir,ref,success,f,link)
	
	f.quit()	
	print "Successfully downloaded %s files and %s ID references were wrong" % (success,failed)	




def downloadFiles(bamPath, ERR_id):

	filename, samfile_extension = os.path.splitext(bamPath)
	dirname = os.path.dirname(filename)
	download_ERR(ERR_id,dirname)







def convertToBAM(samPath):

	filename, samfile_extension = os.path.splitext(samPath)

	os.system("samtools view -buh -o " + filename +'.bam' + " " + samPath)
	#os.system("rm " + samPath)
	os.system("samtools sort " + filename +'.bam' + " " + filename +'_sorted')
	os.system("rm "+ filename +'.bam')
	os.system("samtools index " + filename +'_sorted.bam')

	return (filename +'_sorted')


def rawCoverage(bamSortedPath):

	os.system("bedtools genomecov -d -ibam " + bamSortedPath + ".bam > " + bamSortedPath+".tab")


def checkCoverage(outputPath, coverageThreshold):

	sequenceMedObject = {}
		
	with open(outputPath+'.tab') as tsv:
		prevName = '';
		countlines = 0
		arrayOfpositionValues = []
		sequenceNames = []

		for line in csv.reader(tsv, delimiter="\t"):
			countlines += 1
			if prevName != line[0] and countlines != 1:
				sequenceNames.append(prevName)
				sequenceMedObject[prevName] = [prevName, numpy.average(arrayOfpositionValues), numpy.std(arrayOfpositionValues), arrayOfpositionValues, False, False, False]
				arrayOfpositionValues = []
				prevName = line[0]
			else:
				arrayOfpositionValues.append(int(line[2]))
				if countlines == 1:
					prevName = line[0]

		sequenceNames.append(prevName)
		sequenceMedObject[prevName] = [prevName, numpy.average(arrayOfpositionValues), numpy.std(arrayOfpositionValues), arrayOfpositionValues, False, False, False]

	#with open(outputPath+'_coverageCheck.tab', 'w') as coverageCheckFile:

	for sequence in sequenceMedObject:
		countLowCoverage = 0
		countIndel = 0
		countDuplication = 0
		sequenceLength = len(sequenceMedObject[sequence][3])

		for coverageAtPosition in sequenceMedObject[sequence][3]:
			if coverageAtPosition >= 1.25 * sequenceMedObject[sequence][1]:
				countDuplication += 1
			elif coverageAtPosition < 0.1 * sequenceMedObject[sequence][1]:
				countIndel += 1
			elif coverageAtPosition < int(coverageThreshold):
				countLowCoverage += 1
		
		sequenceMedObject[sequence].append(str(float(countDuplication)/sequenceLength))
		sequenceMedObject[sequence].append(str(float(countIndel)/float(sequenceLength)))
		sequenceMedObject[sequence].append(str(float(countLowCoverage)/float(sequenceLength)))
			#coverageCheckFile.write(sequence + '\t' + str(float(countDuplication)/sequenceLength) + '\t' + str(float(countcountIndel)/float(sequenceLength)) + '\t' + str(float(countLowCoverage)/float(sequenceLength))+"\n")

	return sequenceNames, sequenceMedObject
    		
def alleleCalling(bamSortedPath, referencePath, sequenceNames, gatkPath, sampleID, qualityThreshold, coverageThreshold, multipleAlleles, sequenceMedObject):

	ploidytempFile = bamSortedPath+'_temp_ploi.tab'

	testFile = bamSortedPath + 'test.fasta'

	with open(ploidytempFile, 'w') as tempFile:
		tempFile.write(sampleID + '\t' + str(1))

	os.system("samtools mpileup --no-BAQ --fasta-ref " + referencePath + " --uncompressed -t DP,DPR " + bamSortedPath + ".bam | bcftools call --consensus-caller --gvcf 2 --samples-file "+ ploidytempFile + " --output-type v --output " + bamSortedPath + ".vcf")
	os.system("java -jar " + gatkPath + " -T FastaAlternateReferenceMaker -R "+ referencePath +" -o "+ testFile +" -V "+ bamSortedPath + ".vcf")


	with open(bamSortedPath + ".vcf", 'r') as vcfFile:
		for line in csv.reader(vcfFile, delimiter="\t"):
			if not line[0].startswith("#"):
				if not line[7].startswith('END'):
					quality = line[5]
					if quality < float(qualityThreshold):
						sequenceMedObject[line[0]][4] = True
					
					deepCoverage = float(line[7].split(';')[0].split('=')[1])
					
					if deepCoverage < float(coverageThreshold):
						sequenceMedObject[line[0]][5] = True
					
					coverageByAllele = line[9].split(':')[3].split(',')
					alternativeAlleles = line[4].split(',')

					coverageByAllele = [ int(x) for x in coverageByAllele ]

					dominantSNP = coverageByAllele.index(max(coverageByAllele))
					coverageAllele = coverageByAllele[dominantSNP]
					alternativeSNP = alternativeAlleles[dominantSNP-1]

					if coverageAllele/deepCoverage < float(multipleAlleles):
						sequenceMedObject[line[0]][6] = True

	with open(bamSortedPath + "_mappingCheck.tab", 'w') as mappingCheckFile:
		mappingCheckFile.write('#Sequence\tDuplication\tIndel\tRawCoverage\tAlternativeQualityScore\tCoverage\tMultipleAllele\n')
		for sequence in sequenceMedObject:
			mappingCheckFile.write(sequence + '\t' + str(sequenceMedObject[sequence][7]) + '\t' + str(sequenceMedObject[sequence][8]) + '\t' + str(sequenceMedObject[sequence][9]) + '\t' +str(sequenceMedObject[sequence][4]) + '\t' +str(sequenceMedObject[sequence][5]) + '\t' +str(sequenceMedObject[sequence][6]) + '\n')









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


        
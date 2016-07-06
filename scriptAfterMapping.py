import os
import csv
import numpy
import os.path

import rematch_utils


def convertToBAM(samPath, toClear):

	filename, samfile_extension = os.path.splitext(samPath)

	os.system("samtools view -buh -o " + filename + '_temp.bam' + " " + samPath)
	os.system("rm " + samPath)
	os.system("samtools sort " + filename + "_temp.bam " + filename)
	os.system("rm " + filename + '_temp.bam')
	os.system("samtools index " + filename + '.bam')
	toClear.append(filename + '.bam*')

	return (filename + '')


def rawCoverage(bamSortedPath, toClear):
	command = ['bedtools', 'genomecov', '-d', '-ibam', str(bamSortedPath + '.bam'), '>', str(bamSortedPath + '.tab')]
	run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(command, True, None)
	toClear.append(bamSortedPath + ".tab")


def checkCoverage(outputPath, coverageThreshold, extraSeq, toClear):

	sequenceMedObject = {}

	with open(outputPath + '.tab') as tsv:
		prevName = ''
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
				sequenceMedObject[prevName] = [prevName, numpy.average(arrayOfcoverageValues), numpy.std(arrayOfcoverageValues), arrayOfcoverageValues, False, False, False, arrayOfpositionValues]
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
		sequenceMedObject[prevName] = [prevName, numpy.average(arrayOfcoverageValues), numpy.std(arrayOfcoverageValues), arrayOfcoverageValues, False, False, False, arrayOfpositionValues]
		countSequences += 1

		for sequence in sequenceMedObject:
				countLowCoverage = 0
				countIndel = 0
				countDuplication = 0
				sequenceLength = len(sequenceMedObject[sequence][3])
				sumCoverage = 0

				index = 0
				for coverageAtPosition in sequenceMedObject[sequence][3]:

						if (sequenceMedObject[sequence][7])[index] > extraSeq and (sequenceMedObject[sequence][7])[index] <= (sequenceLength - extraSeq):
								sumCoverage = sumCoverage + coverageAtPosition
								if coverageAtPosition >= 1.25 * sequenceMedObject[sequence][1]:
										countDuplication += 1
								if coverageAtPosition < 0.1 * sequenceMedObject[sequence][1]:
										countIndel += 1
								if coverageAtPosition < int(coverageThreshold):
										countLowCoverage += 1
						index += 1

				sequenceMedObject[sequence].append(float(countDuplication) / float(sequenceLength - (2 * extraSeq)))
				sequenceMedObject[sequence].append(float(countIndel) / float(sequenceLength - (2 * extraSeq)))
				sequenceMedObject[sequence].append(float(countLowCoverage) / float(sequenceLength - (2 * extraSeq)))
				sequenceMedObject[sequence].append(float(sumCoverage) / float(sequenceLength - (2 * extraSeq)))

		return sequenceNames, sequenceMedObject


def alleleCalling(bamSortedPath, referencePath, sequenceNames, gatkPath, sampleID, qualityThreshold, coverageThreshold, multipleAlleles, sequenceMedObject, extraSeq, toClear):

	ploidytempFile = bamSortedPath + '_temp_ploi.tab'
	sequencesFile = bamSortedPath + '_sequences.fasta'

	with open(ploidytempFile, 'w') as tempFile:
		tempFile.write(sampleID + '\t' + str(1))

	print "running bcf"
	command = ['mpileup', '--no-BAQ', '--fasta-ref', referencePath, '--uncompressed', '-t', 'DP,DPR,DV', str(bamSortedPath + '.bam'), '|', 'bcftools', 'call', '--multiallelic-caller', '--variants-only', '--samples-file', ploidytempFile, '--output-type', 'v', '--output', str(bamSortedPath + '.vcf')]
	run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(command, True, None)
	toClear.append(bamSortedPath + ".vcf")
	toClear.append(bamSortedPath + ".vcf.idx")

	rematch_utils.filter_vcf(bamSortedPath + ".vcf", extraSeq, sequenceMedObject)

	with open(bamSortedPath + ".vcf", 'r') as vcfFile:

		for line in csv.reader(vcfFile, delimiter="\t"):
			if not line[0].startswith("#"):
				quality = line[5]

				if quality < float(qualityThreshold):
					sequenceMedObject[line[0]][4] = True

				qualityCoverage = float(line[9].split(':')[2])

				if qualityCoverage < float(coverageThreshold):
					sequenceMedObject[line[0]][5] = True

				coverageByAllele = line[9].split(':')[3].split(',')

				coverageByAllele = [int(x) for x in coverageByAllele]

				dominantSNP = coverageByAllele.index(max(coverageByAllele))
				coverageAllele = coverageByAllele[dominantSNP]

				if coverageAllele / qualityCoverage < float(multipleAlleles):
					sequenceMedObject[line[0]][6] = True

	print "running gatk"
	command = ['java', '-jar', gatkPath, '-T', 'FastaAlternateReferenceMaker', '-R', referencePath, '-o', sequencesFile, '-V', str(bamSortedPath + '.vcf'), '2>>', str(sequencesFile + '_log_gatk.txt')]
	run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(command, True, None)
	toClear.append(sequencesFile + '_log_gatk.txt')

	os.remove(ploidytempFile)

	rematch_utils.changeFastaHeadersAndTrimm(sequencesFile, extraSeq, referencePath)

	rematch_utils.createCheckFile(bamSortedPath, sequenceMedObject)

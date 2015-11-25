import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys

def convertToBAM(samPath):

	filename, samfile_extension = os.path.splitext(samPath)

	os.system("samtools view -buh -o " + filename +'.bam' + " " + samPath)
	os.system("rm " + samPath)
	os.system("samtools sort " + filename +'.bam' + " " + filename +'_sorted')
	os.system("rm "+ filename +'.bam')
	os.system("samtools index " + filename +'_sorted.bam')


def rawCoverage(bamSortedPath, outputPath):

	os.system("bedtools genomecov -d -ibam " + bamSortedPath + " > " + outputPath)


def checkCoverage(outputPath, coverageThreshold):

	import csv

	with open(outputPath) as tsv:
		prevSampleName = '';
		for line in csv.reader(tsv, delimiter="\t"):
			if int(line[2]) < int(coverageThreshold):
				print line[0] + ' low coverage at position: ' + line[1] + ' with value ' + line[2] + ' '

    		

        
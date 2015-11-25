import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys

def convertToBAM(samPath, bamPath, bamSortedPath):

	os.system("samtools view -buh -o " + bamPath + " " + samPath)
	os.system("rm " + samPath)
	os.system("samtools sort " + bamPath + " " + bamSortedPath)
	os.system("rm "+ bamPath)
	os.system("samtools index " + bamSortedPath)


def rawCoverage(bamSortedPath, outputPath):

	os.system("bedtools genomecov -d -ibam " + bamSortedPath + " > " + outputPath)


def checkCoverage(outputPath, coverageThreshold):

	import csv

	with open(outputPath) as tsv:
		prevSampleName = '';
		for line in csv.reader(tsv, delimiter="\t"):
			if int(line[2]) < int(coverageThreshold):
				print line[0] + ' low coverage at position: ' + line[1] + ' with value ' + line[2] + ' '

    		

        
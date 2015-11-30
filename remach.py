import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys

from downloadAndBowtie import downloadAndBowtie
from scriptAfterMapping import checkCoverage
from scriptAfterMapping import convertToBAM
from scriptAfterMapping import rawCoverage
from scriptAfterMapping import alleleCalling


def main():

	parser = argparse.ArgumentParser(description="This program performs a download and read mapping of a given RUN_ID against a reference sequence and returns a report based on coverage.")

	#parser.add_argument('-s', nargs='?', type=str, help="sam file path", required=True)
	parser.add_argument('-r', nargs='?', type=str, help='Path for the reference sequence', required=True)
	parser.add_argument('-t', nargs='?', type=str, help='Output directory path', required=True)
	parser.add_argument('-gatk', nargs='?', type=str, help='Path for the Genome Analysis Toolkit jar file', required=True)
	parser.add_argument('-cov', nargs='?', type=str, help='Minimum coverage', required=True)
	parser.add_argument('-qual', nargs='?', type=str, help='Minimum mapping quality', required=True)
	parser.add_argument('-mul', nargs='?', type=str, help='Multiple alleles', required=True)
	parser.add_argument('-l', nargs='?', type=str, help='list with ids of the sequencing run', required=True)


	args = parser.parse_args()

	runReMaCh(args)

def runReMaCh(args):

	if not os.path.isdir(args.t):
		os.mkdir(args.t)
		print str(args.t) + ' directory created!'

	
	with open(args.l, 'r') as run_ids:

		count_runs = 0
		buildBowtie = True
		for run_id in run_ids:

			count_runs += 1

			if count_runs > 1:
				buildBowtie = False

			run_id = run_id.strip()

			samFilePath = downloadAndBowtie(args.r, run_id, args.t, buildBowtie)
	
			sortedPath = convertToBAM(samFilePath)
			rawCoverage(sortedPath)
			sequenceNames, sequenceMedObject = checkCoverage(sortedPath, args.cov)
			alleleCalling(sortedPath, args.r, sequenceNames, args.gatk, run_id, args.qual, args.cov, args.mul, sequenceMedObject)


if __name__ == "__main__":
	main()
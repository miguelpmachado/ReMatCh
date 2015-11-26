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

	parser = argparse.ArgumentParser(description="This program calls for alleles")

	#parser.add_argument('-s', nargs='?', type=str, help="sam file path", required=True)
	parser.add_argument('-r', nargs='?', type=str, help='reference path', required=True)
	parser.add_argument('-t', nargs='?', type=str, help='targetDir', required=True)
	parser.add_argument('-gatk', nargs='?', type=str, help='gatk jar path', required=True)
	parser.add_argument('-cov', nargs='?', type=str, help='coverage', required=True)
	parser.add_argument('-qual', nargs='?', type=str, help='mapping quality', required=True)
	parser.add_argument('-mul', nargs='?', type=str, help='multiple alleles', required=True)
	parser.add_argument('-sr', nargs='?', type=str, help='run ID', required=True)


	args = parser.parse_args()

	runTest(args)

def runTest(args):
	samFilePath = downloadAndBowtie(args.r, args.sr, args.t)
	
	sortedPath = convertToBAM(samFilePath)
	rawCoverage(sortedPath)
	sequenceNames, sequenceMedObject = checkCoverage(sortedPath, args.cov)
	alleleCalling(sortedPath, args.r, sequenceNames, args.gatk, args.sr, args.qual, args.cov, args.mul, sequenceMedObject)


if __name__ == "__main__":
	main()
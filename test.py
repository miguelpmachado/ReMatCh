import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys

from scriptAfterMapping import checkCoverage
from scriptAfterMapping import convertToBAM
from scriptAfterMapping import rawCoverage
from scriptAfterMapping import alleleCalling


def main():

	parser = argparse.ArgumentParser(description="This program calls for alleles")
	parser.add_argument('-s', nargs='?', type=str, help="sam file path", required=True)
	parser.add_argument('-r', nargs='?', type=str, help='reference path', required=True)
	parser.add_argument('-cov', nargs='?', type=str, help='cov path', required=True)


	args = parser.parse_args()

	runTest(args)

def runTest(args):
	sortedPath = convertToBAM(args.s)
	sequenceNames = rawCoverage(sortedPath)
	checkCoverage(sortedPath, args.cov)
	alleleCalling(sortedPath, referencePath, sequenceNames)


if __name__ == "__main__":
	main()
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


def main():

	parser = argparse.ArgumentParser(description="This program calls for alleles")
	parser.add_argument('-s', nargs='?', type=str, help="sam file path", required=True)
	parser.add_argument('-b', nargs='?', type=str, help="bam file path", required=True)
	parser.add_argument('-bs', nargs='?', type=str, help='bam sorted path', required=True)
	parser.add_argument('-t', nargs='?', type=str, help='coverage tab file', required=True)
	parser.add_argument('-cov', nargs='?', type=str, help='cov path', required=True)


	args = parser.parse_args()

	runTest(args)

def runTest(args):
	convertToBAM(args.s, args.b, args.bs)
	rawCoverage(args.bs, args.t)
	checkCoverage(args.t, args.cov)


if __name__ == "__main__":
	main()
import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys
import glob

from downloadAndBowtie import downloadAndBowtie
from scriptAfterMapping import checkCoverage
from scriptAfterMapping import convertToBAM
from scriptAfterMapping import rawCoverage
from scriptAfterMapping import alleleCalling
from SeqFromWebTaxon import GetSequencesFromTaxon

def main():

	parser = argparse.ArgumentParser(description="This program performs a download and read mapping of a given RUN_ID against a reference sequence and returns a report based on coverage.")

	#parser.add_argument('-s', nargs='?', type=str, help="sam file path", required=True)
	parser.add_argument('-r', nargs='?', type=str, help='Path for the reference sequence', required=True)
	parser.add_argument('-t', nargs='?', type=str, help='Output directory path', required=True)
	parser.add_argument('-gatk', nargs='?', type=str, help='Path for the Genome Analysis Toolkit jar file', required=True)
	parser.add_argument('-picard', nargs='?', type=str, help='Path for the Picard jar file', required=True)
	parser.add_argument('-cov', nargs='?', type=int, help='Minimum coverage', required=True)
	parser.add_argument('-qual', nargs='?', type=float, help='Minimum mapping quality', required=True)
	parser.add_argument('-mul', nargs='?', type=float, help='Multiple alleles', required=True)
	parser.add_argument('-threads', nargs='?', type=int, help='Number of threads', required=False, default= 1)
	parser.add_argument('-rmFastq', nargs='?', type=bool, help='Remove fastq files after the analysis', required=False, default = False)
	parser.add_argument('-l', nargs='?', type=str, help='Path to a list with ids of the sequencing run', required=True)
	parser.add_argument('-tax', nargs='?', type=str, help='Name taxon to download sequences', required=False)
	parser.add_argument('-allplat', nargs='?', type=bool, help='Use all platforms', required=False, default = False)

	args = parser.parse_args()

	runReMaCh(args)

def runReMaCh(args):

	if not os.path.isdir(args.t):
		os.mkdir(args.t)
		print str(args.t) + ' directory created!'
	
	if args.tax:
		GetSequencesFromTaxon(args.tax,args.l,True)
	
	if args.allplat == False:
		platform="Illumina"
	else:
		platform=''
	
	
	with open(args.l, 'r') as run_ids:

		count_runs = 0
		buildBowtie = True
		firstLine = True
		for run_id in run_ids:
			run=False
			if args.tax and firstLine == True:
				firstLine = False
				continue
			
			elif args.tax and platform:
					run_info=run_id.split("\t")
					run_id=run_info[0]
					run_plat=run_info[1]
					
					if platform in run_plat and not "Analyzer" in run_plat:
						run=True
			
			elif args.tax:
				run_info=run_id.split("\t")
				run_id=run_info[0]
				run_plat=run_info[1]
						
				run=True
			
			else:
				run=True
				
			if run==True:	

				startTime = datetime.now()	
				
				count_runs += 1

				if count_runs > 1:
					buildBowtie = False

				run_id = run_id.strip()

				samFilePath = downloadAndBowtie(args.r, run_id, args.t, buildBowtie, args.picard, args.threads)
		
				sortedPath = convertToBAM(samFilePath)
				rawCoverage(sortedPath)
				sequenceNames, sequenceMedObject, sequenceAndIndex = checkCoverage(sortedPath, args.cov)
				alleleCalling(sortedPath, args.r, sequenceNames, args.gatk, run_id, args.qual, args.cov, args.mul, sequenceMedObject, sequenceAndIndex)

				gzSizes = 0

				filesToRemove = glob.glob(os.path.join(args.t, run_id) + '/*.fastq.gz')

				for files in filesToRemove:
					gzSizes += float(os.path.getsize(files))

				if args.rmFastq == True:
					for i in filesToRemove:
						os.remove(i)

				run_time = str(datetime.now() - startTime)

				with open(os.path.join(args.t, run_id, run_id + '_runtime.txt'), 'w') as runTimeFile:
					runTimeFile.write("#runTime\tfileSize\n")
					runTimeFile.write(str(run_time) + '\t' + str(gzSizes) + '\n')


if __name__ == "__main__":
	main()

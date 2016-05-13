#!/usr/bin/env python

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
from rematch_utils import removeFromArray

def main():

	parser = argparse.ArgumentParser(prog='rematch.py', description="ReMatCh is an application which combines a set of bioinformatic tools for reads mapping against a reference, finds the allelic variants and produces a consensus sequence. It also allows direct sample download from ENA database to be used in the analysis.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	requiredNamed = parser.add_argument_group('required arguments')
	#parser.add_argument('-s', nargs='?', type=str, help="sam file path", required=True)
	requiredNamed.add_argument('-r', nargs=1, type=str, metavar=('/path/reference.fasta'), help='Path for the reference sequence', required=True)
	requiredNamed.add_argument('-d', '--workdir', nargs=1, type=str, metavar=('/path/to/workdir'), help='Working directory. Downloaded files will be stored here, but it can also already contain folders with fastq files. Results will be stored here.', required=True)
	requiredNamed.add_argument('-gatk', nargs=1, type=str, metavar=('/path/to/gatk.jar'), help='Path for the Genome Analysis Toolkit jar file', required=True)
	requiredNamed.add_argument('-picard', nargs=1, metavar=('/path/to/picard'), type=str, help='Path for Picard', required=True)
	requiredNamed.add_argument('-l', nargs=1, metavar=('/path/to/idenfifiersList.txt'), type=str, help='Path to a list with ids to run. IDs can be ENA run accession numbers for download or direcory names where fastqs are stored in --workdir. Run accession numbers retrieved from ENA using -tax will be stored here.' , required=True)
	parser.add_argument('-cov', '--min-coverage', nargs='?', metavar=('N'), type=int, help='Minimum coverage depth required for base calling and SNP calling.', default = 10, required=False)
	parser.add_argument('-qual', '--min-quality', nargs='?', metavar=('N'), type=int, help='Minimum mapping quality for SNP calling', default = 10, required=False)
	parser.add_argument('-mul', nargs='?', metavar=('0.0 - 1.0'), type=float, help='Minimum read mapping frequency to be considered as ', default = 0.75, required=False)
	parser.add_argument('-threads', nargs='?', metavar=('N'), type=int, help='Number of threads used to run bowtie2', required=False, default= 1)
	parser.add_argument('-tax', nargs='?', metavar=('Streptococcus pneumoniae'), type=str, help='Name taxon to download sequences. Results will be stored in /path/to/idenfifiersList.txt', required=False)
	parser.add_argument('-xtraSeq', nargs='?', type=int, help='For trimming extra sequence lenght 5\' and 3\' ', required=False, default = 0)
	parser.add_argument('-bowtieBuild', help='Run build bowtie', action='store_true')
	parser.add_argument('-clean', help='Clean intermediate files produced by the application (.bam, .vcf, index files, coverage file)', action='store_true')
	parser.add_argument('-rmFastq', help='Remove fastq files after the analysis', action='store_true')
	parser.add_argument('-allplat', help='Use all platforms. By default, only Illumina runs are used', action='store_true')

	args = parser.parse_args()

	runReMaCh(args)

def runReMaCh(args):

	toClear = []

	if not os.path.isdir(args.t):
		os.mkdir(args.t)
		print str(args.t) + ' directory created!'
	
	if args.tax:
		GetSequencesFromTaxon(args.tax,args.l,True)
	
	if not args.allplat:
		platform="Illumina"
	else:
		platform=''
	
	ids_with_problems = open(os.path.join(args.t,'ids_with_problems.txt'), 'w')
	logFile = open(os.path.join(args.t,'log_file.txt'), 'a')

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

				if count_runs > 1 or not args.bowtieBuild:
					buildBowtie = False

				run_id = run_id.strip()


				print "\n######\ndownloading and bowtieying\n######\n"
				logFile.write("\n######\ndownloading and bowtieying\n######\n")
				samFilePath, singOrPaired, numFilesDownloaded = downloadAndBowtie(args.r, run_id, args.t, buildBowtie, args.picard, args.threads, logFile, toClear)
				print "\n######\ndownloaded and bowtied\n######\n"
				logFile.write("\n######\ndownloaded and bowtied\n######\n")
				
				if numFilesDownloaded == 0:
					ids_with_problems.write(run_id + '\n')
					pass

				if not samFilePath==False:
				
					sortedPath = convertToBAM(samFilePath, toClear)

					rawCoverage(sortedPath, toClear)
					print "\n######\nChecking coverage\n######\n"
					logFile.write("\n######\nChecking coverage\n######\n")
					sequenceNames, sequenceMedObject = checkCoverage(sortedPath, args.cov,args.xtraSeq, logFile, toClear)
					print "\n######\nChecked coverage goint to perform allele call\n######\n"
					logFile.write("\n######\nChecked coverage goint to perform allele call\n######\n")
					alleleCalling(sortedPath, args.r, sequenceNames, args.gatk, run_id, args.qual, args.cov, args.mul, sequenceMedObject,args.xtraSeq, logFile, toClear)
					print "\n######\nallele called everything\n######\n"	
					logFile.write("\n######\nallele called everything\n######\n")
					
					
					gzSizes = 0

					filesToRemove = glob.glob(os.path.join(args.t, run_id) + '/*.fastq.gz')

					for files in filesToRemove:
						gzSizes += float(os.path.getsize(files))

					if args.rmFastq:
						for i in filesToRemove:
							os.remove(i)

					run_time = str(datetime.now() - startTime)

					with open(os.path.join(args.t, run_id, run_id + '_runtime.txt'), 'w') as runTimeFile:
						runTimeFile.write("#runTime\tfileSize\tlibraryLayout\n")
						runTimeFile.write(str(run_time) + '\t' + str(gzSizes) +"\t"+singOrPaired+ '\n')
				else:
					print run_id+" not downloaded sucessfully"
					logFile.write(run_id+" not downloaded sucessfully" + '\n')
					pass

		ids_with_problems.close()

	if args.clean:
		removeFromArray(toClear)


if __name__ == "__main__":
	main()

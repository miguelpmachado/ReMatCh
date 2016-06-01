#!/usr/bin/env python

import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys
import glob
import time

from downloadAndBowtie import downloadAndBowtie
from scriptAfterMapping import checkCoverage
from scriptAfterMapping import convertToBAM
from scriptAfterMapping import rawCoverage
from scriptAfterMapping import alleleCalling
from SeqFromWebTaxon import GetSequencesFromTaxon
from rematch_utils import removeFromArray
from rematch_utils import checkPrograms
from mergeResults import mergeResults
from rematch_utils import Logger
from rematch_utils import removeIndexes


def main():

	parser = argparse.ArgumentParser(prog='rematch.py', description="ReMatCh is an application which combines a set of bioinformatic tools for reads mapping against a reference, finds the allelic variants and produces a consensus sequence. It also allows direct sample download from ENA database to be used in the analysis.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	requiredNamed = parser.add_argument_group('required arguments')
	
	#parser.add_argument('-s', nargs='?', type=str, help="sam file path", required=True)
	requiredNamed.add_argument('-r', nargs="?", type=str, metavar=('/path/reference.fasta'), help='Path for the reference sequence', required=False)
	requiredNamed.add_argument('-d', '--workdir', nargs="?", type=str, metavar=('/path/to/workdir'), help='Working directory. Downloaded files will be stored here under sampleID/fastq/, but it can also already contain folders with fastq files. Results will be stored here.', required=False)
	requiredNamed.add_argument('-gatk', nargs="?", type=str, metavar=('/path/to/gatk.jar'), help='Path for the Genome Analysis Toolkit jar file', required=False)
	requiredNamed.add_argument('-picard', nargs="?", metavar=('/path/to/picard'), type=str, help='Path for Picard', required=False)
	requiredNamed.add_argument('-l', nargs="?", metavar=('/path/to/identifiersList.txt'), type=str, help='Path to a list with ids to run. IDs can be ENA run accession numbers for download or directory names where fastqs are stored in --workdir. Run accession numbers retrieved from ENA using -tax will be stored here.' , required=False)
	parser.add_argument('-cov', '--minCoverage', nargs='?', metavar=('N'), type=int, help='Minimum coverage depth required for base calling and SNP calling.', default = 10, required=False)
	parser.add_argument('-qual', '--minQuality', nargs='?', metavar=('N'), type=int, help='Minimum mapping quality for SNP calling', default = 10, required=False)
	parser.add_argument('-mul', '--multipleAlleles', nargs='?', metavar=('0.0 - 1.0'), type=float, help='Minimum reads frequency (confidence) of dominant nucleotide to consider absence of multiple alleles at a given SNP position.', default = 0.75, required=False)
	parser.add_argument('-threads', nargs='?', metavar=('N'), type=int, help='Number of threads used to run bowtie2', required=False, default= 1)
	parser.add_argument('-tax', nargs='?', metavar=('"Streptococcus pneumoniae"'), type=str, help='Name taxon to download sequences. Results will be stored in /path/to/idenfifiersList.txt', required=False)
	parser.add_argument('-xtraSeq', nargs='?', type=int, help='For trimming extra sequence lenght 5\' and 3\' ', required=False, default = 0)
	parser.add_argument('--asperaKey', nargs=1, type=str, metavar='/path/to/asperaweb_id_dsa.openssh', help='Tells ReMatCh to download run files from ENA using Aspera Connect. The path to Private-key file asperaweb_id_dsa.openssh normaly found in ~/.aspera/connect/etc/asperaweb_id_dsa.openssh needs to be provided.', required=False, default=None) ## mpmachado ##
	parser.add_argument('-bowtieBuild', help='Run build bowtie', action='store_true')
	parser.add_argument('-clean', help='Clean intermediate files produced by the application (.bam, .vcf, index files, coverage file)', action='store_true')
	parser.add_argument('-rmFastq', help='Remove fastq files after the analysis', action='store_true')
	parser.add_argument('-allplat', help='Use all platforms. By default, only Illumina runs are used', action='store_true')

	#Merge results
	mergedResults = parser.add_argument_group('merge results arguments. To be used after ReMatCh run')
	mergedResults.add_argument('--mergeResults', nargs=1, metavar=('/path/to/workdir'), type=str, help='Merge all ReMatCh results available at --workdir. Option to be used alone or with --sequenceCoverage.', required=False)
	mergedResults.add_argument('--sequenceCoverage', nargs=1, metavar=('0.0 - 1.0'), type=float, help='Minimum sequence length to consider the gene to be present. This is a relative measure. To be used with --mergeResults', default=0.8, required=False)

	args = parser.parse_args()

	if args.mergeResults:
		print 'Running mergeResults with --mergeResults ' + str(args.mergeResults[0]) + ' and --sequenceCoverage '+ str(args.sequenceCoverage)
		mergeResults(args.mergeResults[0], args.sequenceCoverage)
	else:
		if not args.r or not args.workdir or not args.gatk or not args.picard or not args.l:
			parser.error('You must pass all the required arguments. For more information type -h')
		else:
			if not os.path.isdir(args.workdir):
				os.mkdir(args.workdir)
				print str(args.workdir) + ' directory created!'
			
			start_time = time.time()
			
			sys.stdout = Logger(args.workdir)
			checkPrograms(args)
			
			# Print arguments passed and shell PATH variable
			print '\n' + 'COMMAND:'
			print sys.executable + ' ' + os.path.abspath(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:])
			print '\n' + 'PATH variable:'
			print os.environ['PATH']
			
			runReMaCh(args)
			
			end_time = time.time()
			time_taken = end_time - start_time
			hours, rest = divmod(time_taken,3600)
			minutes, seconds = divmod(rest, 60)
			print ">>> Runtime :" + str(hours) + "h:" + str(minutes) + "m:" + str(round(seconds, 2)) + "s" + "\n"

def runReMaCh(args):

	toClear = []

	if args.tax:
		GetSequencesFromTaxon(args.tax,args.l,True)
	
	if not args.allplat:
		platform="Illumina"
	else:
		platform=''
	
	ids_with_problems = open(os.path.join(args.workdir,'ids_with_problems.txt'), 'w')
	ids_no_problems = open(os.path.join(args.workdir,'ids_no_problems.txt'), 'w')
	#logFile = open(os.path.join(args.workdir,'log_file.txt'), 'a')

	with open(args.l, 'r') as run_ids:

		count_runs = 0
		buildBowtie = True
		firstLine = True
		for run_id in run_ids:
			toClear = []
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

				print "\nRunning ID: " + str(run_id)
				samFilePath, singOrPaired, filesDownloaded = downloadAndBowtie(args.r, run_id, args.workdir, buildBowtie, args.picard, args.threads, toClear, args.asperaKey) # mpmachado #
				#print "\n######\ndownloaded and bowtied\n######\n"
				#logFile.write("\n######\ndownloaded and bowtied\n######\n")
				
				if not samFilePath==False:
				
					sortedPath = convertToBAM(samFilePath, toClear)

					rawCoverage(sortedPath, toClear)
					print "Checking coverage..."
					#logFile.write("\nChecking coverage...")
					sequenceNames, sequenceMedObject = checkCoverage(sortedPath, args.minCoverage,args.xtraSeq, toClear)
					print "Performing Allele Call..."
					#logFile.write("\nPerforming Allele Call...")
					alleleCalling(sortedPath, args.r, sequenceNames, args.gatk, run_id, args.minQuality, args.minCoverage, args.multipleAlleles, sequenceMedObject,args.xtraSeq, toClear)
					print str(run_id) + " DONE" 
					#logFile.write(str(run_id) + " DONE")
					
					
					gzSizes = 0

					#filesToRemove = glob.glob(os.path.join(args.workdir, run_id, 'downloads'))

					for files in filesDownloaded:
						gzSizes += float(os.path.getsize(os.path.join(args.workdir, run_id, 'fastq',files)))

					if args.rmFastq:
						os.system('rm -r ' + os.path.join(args.workdir, run_id, 'fastq'))
					
					ids_no_problems.write(run_id + "\n")
					ids_no_problems.flush()

					run_time = str(datetime.now() - startTime)

					with open(os.path.join(args.workdir, run_id, run_id + '_runtime.txt'), 'w') as runTimeFile:
						runTimeFile.write("#runTime\tfileSize\tlibraryLayout\n")
						runTimeFile.write(str(run_time) + '\t' + str(gzSizes) +"\t"+singOrPaired+ '\n')
				else:
					ids_with_problems.write(run_id + "\n")
					ids_with_problems.flush()
					if args.rmFastq:
						try:
							shutil.rmtree(os.path.join(args.workdir, run_id, 'fastq'))
						except Exception as e:
							print e
					print run_id + ' - An error has occurred.'
			
			if args.clean:
				removeFromArray(toClear)
		
		ids_with_problems.close()
		ids_no_problems.close()

		if args.clean:
			removeIndexes(args.r)




if __name__ == "__main__":
	main()

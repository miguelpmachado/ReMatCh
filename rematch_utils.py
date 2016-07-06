import os
import sys
import csv
import shlex, subprocess
import os.path
import time

import argparse
from rematch import runReMatCh
from mergeResults import runMergeResults

from threading import Timer


def parseArguments(version):
	parser = argparse.ArgumentParser(prog='rematch.py', description="ReMatCh is an application which combines a set of bioinformatic tools for reads mapping against a reference, finds the allelic variants and produces a consensus sequence. It also allows direct sample download from ENA database to be used in the analysis.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers()

	parser_rematch = subparsers.add_parser('ReMatCh')

	rematch_required = parser_rematch.add_argument_group('ReMatCh required options')
	rematch_required.add_argument('-r', '--reference', nargs=1, type=argparse.FileType('r'), metavar=('/path/reference.fasta'), help='Path for the reference sequence', required=True)
	rematch_required.add_argument('-d', '--workdir', nargs=1, type=str, metavar=('/path/to/workdir'), help='Working directory. Downloaded files will be stored here under sampleID/fastq/, but it can also already contain folders with fastq files. Results will be stored here.', required=True)
	rematch_required.add_argument('--gatk', nargs=1, type=str, metavar=('/path/to/gatk.jar'), help='Path for the Genome Analysis Toolkit jar file', required=True)
	rematch_required.add_argument('--picard', nargs=1, metavar=('/path/to/picard'), type=str, help='Path for Picard', required=True)
	rematch_required.add_argument('-l', nargs="?", metavar=('/path/to/identifiersList.txt'), type=str, help='Path to a list with ids to run. IDs can be ENA run accession numbers for download or directory names where fastqs are stored in --workdir. Run accession numbers retrieved from ENA using -tax will be stored here.', required=True)

	rematch_optional = parser_rematch.add_argument_group('ReMatCh optional options')
	rematch_optional.add_argument('-cov', '--minCoverage', nargs='?', metavar=('N'), type=int, help='Minimum coverage depth required for base calling and SNP calling.', default=10, required=False)
	rematch_optional.add_argument('-qual', '--minQuality', nargs='?', metavar=('N'), type=int, help='Minimum mapping quality for SNP calling', default=10, required=False)
	rematch_optional.add_argument('-mul', '--multipleAlleles', nargs='?', metavar=('0.0 - 1.0'), type=float, help='Minimum reads frequency (confidence) of dominant nucleotide to consider absence of multiple alleles at a given SNP position.', default=0.75, required=False)
	rematch_optional.add_argument('-threads', nargs='?', metavar=('N'), type=int, help='Number of threads used to run bowtie2', required=False, default=1)
	rematch_optional.add_argument('--tax', nargs=1, metavar=('"Streptococcus pneumoniae"'), type=str, help='Name taxon to download sequences. Results will be stored in /path/to/idenfifiersList.txt', required=False, default=[None])
	rematch_optional.add_argument('-xtraSeq', nargs='?', type=int, help='For trimming extra sequence lenght 5\' and 3\' ', required=False, default=0)
	rematch_optional.add_argument('--asperaKey', nargs=1, type=str, metavar='/path/to/asperaweb_id_dsa.openssh', help='Tells ReMatCh to download run files from ENA using Aspera Connect. The path to Private-key file asperaweb_id_dsa.openssh normaly found in ~/.aspera/connect/etc/asperaweb_id_dsa.openssh needs to be provided.', required=False, default=None)
	rematch_optional.add_argument('-bowtieBuild', help='Run build bowtie', action='store_true')
	rematch_optional.add_argument('-clean', help='Clean intermediate files produced by the application (.bam, .vcf, index files, coverage file)', action='store_true')
	rematch_optional.add_argument('-rmFastq', help='Remove fastq files after the analysis', action='store_true')
	rematch_optional.add_argument('-allplat', help='Use all platforms. By default, only Illumina runs are used', action='store_true')
	rematch_optional.add_argument('--useOmicsDataType', nargs=1, type=useOmicsDataType, metavar='GENOMIC,TRANSCRIPTOMIC', help='Tells ReMatCh to analyse these OMICS data type', required=False, default=['ALL'])
	rematch_optional.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	parser_rematch.set_defaults(func=runReMatCh)

	parser_mergeResults = subparsers.add_parser('mergeResults')

	mergeResults_required = parser_mergeResults.add_argument_group('mergeResults required options')
	mergeResults_required.add_argument('--mrWorkdir', nargs=1, metavar=('/path/to/workdir'), type=str, help='Merge all ReMatCh results available at --workdir. Option to be used alone or with --sequenceCoverage.', required=True)

	mergeResults_optional = parser_mergeResults.add_argument_group('mergeResults optional options')
	mergeResults_optional.add_argument('--mrSequenceCoverage', nargs=1, metavar=('0.0 - 1.0'), type=float, help='Minimum sequence length to consider the gene to be present. This is a relative measure. To be used with --mrWorkdir', default=[0.8], required=False)
	mergeResults_optional.add_argument('--mrOutdir', nargs=1, metavar=('/path/to/Results/Outdir/'), type=str, help='Specify a different directory from --workdir to output the merged results (otherwise merged results will be stored in --workdir). To be used with --mrWorkdir', required=False, default=[None])

	parser_mergeResults.set_defaults(func=runMergeResults)

	args = parser.parse_args()

	return args


# For parseArguments
def useOmicsDataType(arguments):
	possible_choises = ['GENOMIC', 'TRANSCRIPTOMIC', 'SYNTHETIC', 'ALL']
	arguments = arguments.split(',')
	for argument in arguments:
		if argument not in possible_choises:
			print 'Choose from ' + ', '.join(possible_choises)
			argparse.ArgumentParser.error('Choose from ' + ', '.join(possible_choises))
	if len(arguments) > 1 and 'ALL' in arguments:
		print 'ALL data type should be given alone'
		argparse.ArgumentParser.error('ALL data type should be given alone')
	else:
		return list(set(arguments))


def check_create_directory(directory):
	if not os.path.isdir(directory):
		os.makedirs(directory)


def start_logger(workdir):
	sys.stdout = Logger(workdir, time.strftime("%Y%m%d-%H%M%S"))
	logfile = sys.stdout.getLogFile()
	return logfile


def scriptVersionGit(version, directory, script_path):
	print 'Version ' + version
	os.chdir(os.path.dirname(script_path))
	command = ['git', 'log', '-1', '--date=local', '--pretty=format:"%h (%H) - Commit by %cn, %cd) : %s"']
	run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None)
	print stdout
	command = ['git', 'remote', 'show', 'origin']
	run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None)
	print stdout
	os.chdir(directory)


def general_information(logfile, version):
	# Check if output directory exists

	print '\n' + '==========> ReMatCh <=========='
	print '\n' + 'Program start: ' + time.ctime()

	# Tells where the logfile will be stored
	print '\n' + 'LOGFILE:'
	print logfile

	# Print command
	print '\n' + 'COMMAND:'
	script_path = os.path.abspath(sys.argv[0])
	print sys.executable + ' ' + script_path + ' ' + ' '.join(sys.argv[1:])

	# Print directory where programme was lunch
	print '\n' + 'PRESENT DIRECTORY :'
	present_directory = os.path.abspath(os.getcwd())
	print present_directory

	# Print program version
	print '\n' + 'VERSION INNUca.py:'
	scriptVersionGit(version, present_directory, script_path)

	# Print PATH variable
	print '\n' + 'PATH variable:'
	print os.environ['PATH']


class Logger(object):
	def __init__(self, out_directory, time_str):
		self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
		self.terminal = sys.stdout
		self.log = open(self.logfile, "w")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)
		self.log.flush()

	def flush(self):
		pass

	def getLogFile(self):
		return self.logfile


def runCommandPopenCommunicate(command, shell_True, timeout_sec_None):
	run_successfully = False
	if isinstance(command, basestring):
		command = shlex.split(command)
	else:
		command = shlex.split(' '.join(command))

	print 'Running: ' + ' '.join(command)
	if shell_True:
		command = ' '.join(command)
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	else:
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if timeout_sec_None is None:
		stdout, stderr = proc.communicate()
	else:
		timer = Timer(timeout_sec_None, proc.kill)
		timer.start()
		stdout, stderr = proc.communicate()
		timer.cancel()

	if proc.returncode == 0:
		run_successfully = True
	else:
		print 'STDOUT'
		print stdout.decode("utf-8")
		print 'STDERR'
		print stderr.decode("utf-8")
	return run_successfully, stdout, stderr


def createCheckFile(bamSortedPath, sequenceMedObject):

	check_fileName = bamSortedPath.replace('_sorted', '')

	with open(check_fileName + "_mappingCheck.tab", 'w') as mappingCheckFile:
		mappingCheckFile.write('#Sequence\tDuplication\tIndel\tCoverage\tLowSNPsQualityScore\tSNPCoverage\tSNPMultipleAllele\tmeanSequenceCoverage\tmeanTrimmedSequenceCoverage\n')
		for sequence in sequenceMedObject:
			mappingCheckFile.write(sequence + '\t' + str(round(sequenceMedObject[sequence][8], 3)) + '\t' + str(round(sequenceMedObject[sequence][9], 3)) + '\t' + str(round(sequenceMedObject[sequence][10], 3)) + '\t' + str(sequenceMedObject[sequence][4]) + '\t' + str(sequenceMedObject[sequence][5]) + '\t' + str(sequenceMedObject[sequence][6]) + '\t' + str(round(sequenceMedObject[sequence][1], 3)) + '\t' + str(round(sequenceMedObject[sequence][11], 3)) + '\n')


def filter_vcf(pathToVcf, extraSeq, sequenceMedObject):

	with open(pathToVcf, 'r') as vcfFile:
		with open(pathToVcf + '.temp', 'w') as vcfTemp:

			for line in csv.reader(vcfFile, delimiter="\t"):
				if not line[0].startswith("#"):
					if int(line[1]) > extraSeq and int(line[1]) <= len(sequenceMedObject[line[0]][3]) - extraSeq:
						vcfTemp.write('\t'.join([str(x) for x in line]) + '\n')
				else:
					vcfTemp.write('\t'.join([str(x) for x in line]) + '\n')

	os.remove(pathToVcf)
	os.rename(pathToVcf + '.temp', pathToVcf)


def changeFastaHeadersAndTrimm(FastasequencesFile, TrimmExtraSeq, referencePath):

	headersArray = []

	with open(referencePath, 'r') as seqFileRef:
		for line in seqFileRef:
			if '>' in line:
				headersArray.append(line)

	with open(FastasequencesFile, 'r') as seqFile:
		with open(FastasequencesFile + ".temp", 'w') as tempFile:
			tempStr = ''
			headercounter = 0
			for line in seqFile:

				if '>' in line:

					if TrimmExtraSeq != 0 and len(tempStr) > 0:
						tempStr = tempStr[TrimmExtraSeq:len(tempStr) - TrimmExtraSeq]
					if headercounter > 0:
						tempFile.write(tempStr + "\n")

					tempStr = ''

					tempFile.write(headersArray[headercounter])
					headercounter += 1

				else:
					tempStr += line.replace('\n', '').replace('\r', '')

			if TrimmExtraSeq != 0 and len(tempStr) > 0:
				tempStr = tempStr[TrimmExtraSeq:len(tempStr) - TrimmExtraSeq]
			tempFile.write(tempStr)

	os.remove(FastasequencesFile)
	os.rename(FastasequencesFile + ".temp", FastasequencesFile)


def removeFromArray(toClear):
	for i in toClear:
		os.system('rm ' + i)


def checkPrograms(args):

	print '\nChecking dependencies...'
	which_program = ['which', '']
	programs = {'bedtools': ['>=', '2.22'], 'java': ['>=', '1.8'], 'samtools': ['=', '1.2'], 'bcftools': ['=', '1.2'], 'bowtie2': ['>=', '2.2.6'], 'ascp': ['>=', '3.6.1']}
	listMissings = []
	for program in programs:
		if program == 'ascp' and not args.asperaKey:
			print 'No aspera key defined, using ftp.'
			continue
		which_program[1] = program
		run_successfully, stdout, stderr = runCommandPopenCommunicate(which_program, False, None)
		if not run_successfully:
			listMissings.append(program + ' not found in PATH.')
		else:
			if program == 'java':
				check_version = [stdout.strip('\n'), '-version']
			else:
				check_version = [stdout.strip('\n'), '--version']

			print program + ' found at: ' + check_version[0]
			run_successfully, stdout, stderr = runCommandPopenCommunicate(check_version, False, None)
			if program == 'java':
				stdout = stderr
			version_line = stdout.split('\n')[0].split(' ')[-1]
			version_line = version_line.replace('"', '')
			if 'v' in version_line:
				version_line = version_line.split('v')[1]
			elif 'V' in version_line:
				version_line = version_line.split('V')[1]
			if programs[program][0] == '>=':
				program_found_version = version_line.split('.')
				program_version_required = programs[program][1].split('.')
				if float('.'.join(program_found_version[0:2])) < float('.'.join(program_version_required[0:2])):
					listMissings.append('ReMatCh requires ' + program + ' with version ' + programs[program][1] + ' or above.')
				elif float('.'.join(program_found_version[0:2])) == float('.'.join(program_version_required[0:2])):
					if len(program_version_required) == 3:
						if len(program_found_version) == 2:
							program_found_version.append(0)
						if program_found_version[2] < program_version_required[2]:
							listMissings.append('ReMatCh requires ' + program + ' with version ' + programs[program][1] + ' or above.')
			else:
				if version_line != programs[program][1]:
					listMissings.append('ReMatCh requires ' + program + ' with version ' + programs[program][1] + '.')

	if len(listMissings) > 0:
		sys.exit('\nErrors:\n' + '\n'.join(listMissings) + '\n\nInstall all dependencies and try again.')


def removeIndexes(referencePath):

	referenceFileName, extension = os.path.splitext(referencePath)

	if os.path.isfile(referenceFileName + '_picard_out.txt'):
		os.remove(referenceFileName + '_picard_out.txt')
	if os.path.isfile(referenceFileName + '.dict'):
		os.remove(referenceFileName + '.dict')
	if os.path.isfile(referencePath + '.fai'):
		os.remove(referencePath + '.fai')
	try:
		os.system('rm ' + referenceFileName + ".*.bt2")
	except:
		print 'It is not possible to remove ".*.bt2" files'

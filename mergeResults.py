import os
import csv
import sys
import rematch_utils


def runMergeResults(args, version):

	# Setting working directory
	if args.mrOutdir[0] is None:
		workdir = os.path.abspath(os.path.join(args.mrWorkdir[0], ''))
	else:
		workdir = os.path.abspath(os.path.join(args.mrOutdir[0], ''))

	# Create working directory
	rematch_utils.check_create_directory(workdir)

	# Start logger
	logfile = rematch_utils.start_logger(workdir)

	# Get general information
	rematch_utils.general_information(logfile, version)

	mergeResults(args.mrWorkdir[0], args.mrSequenceCoverage[0], workdir)


def mergeResults(workdir, sequenceCoverage, outdir):

	# sampleList = os.listdir(workdir)
	dirs = [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d))]

	print '\n' + 'Analysing ' + str(len(dirs)) + ' samples' + '\n'

	sampledict = {}
	consensusdict = {}
	prevNameSeq = ''
	countdirs = 0

	for i in dirs:
		countdirs += 1
		# time.sleep(0.1)
		sys.stderr.write("\rChecking " + str(countdirs) + " of " + str(len(dirs)) + " folders...")
		sys.stderr.flush()
		sampleName = os.path.basename(i)
		mappingFilePath = os.path.join(workdir, sampleName, 'rematch_results', sampleName + '_mappingCheck.tab')
		consensusFilePath = os.path.join(workdir, sampleName, 'rematch_results', sampleName + '_sequences.fasta')
		# parse mappingCheck file
		if os.path.isfile(mappingFilePath):
			with open(mappingFilePath, 'r') as csvfile:
				mappingreader = csv.reader(csvfile, delimiter='\t')
				sampledict[sampleName] = {}
				for line in mappingreader:
					if not line[0].startswith("#"):
						line[0] = line[0].replace(' ', '_').replace('/', '-')
						if 1 - float(line[3]) >= sequenceCoverage:
							if line[6] == 'False':
								sampledict[sampleName][line[0]] = line[8]
							else:
								sampledict[sampleName][line[0]] = 'Mul_Allele (' + line[8] + ')'
						else:
							sampledict[sampleName][line[0]] = 'Absent (' + line[8] + ')'

		# parse consensus sequences
		if os.path.isfile(consensusFilePath):
			with open(consensusFilePath, 'r') as consensusFile:
				countSequences = -1
				for line in consensusFile:
					line = line.splitlines()[0]
					if line.startswith('>'):
						nameseq = line[1:]
						prevNameSeq = nameseq
						nameseq = nameseq.replace(' ', '_').replace('/', '-')
						if nameseq not in consensusdict:
							consensusdict[nameseq] = []
							consensusdict[nameseq].append(['>' + sampleName + '-' + prevNameSeq, ''])
						else:
							consensusdict[nameseq].append(['>' + sampleName + '-' + prevNameSeq, ''])
						countSequences = len(consensusdict[nameseq]) - 1
					else:
						consensusdict[nameseq][countSequences][1] = consensusdict[nameseq][countSequences][1] + line

	print "\nWriting results..."

	if not os.path.exists(os.path.join(outdir, 'merged_results')):
		os.makedirs(os.path.join(outdir, 'merged_results'))

	with open(os.path.join(outdir, 'merged_results', 'mergedResults.tab'), 'w') as results:
		samples = sampledict.keys()
		genes = sampledict[samples[0]].keys()

		header = ['#Samples'] + genes
		results.write('\t'.join(header) + '\n')

		for i in samples:
			row = [i]
			for j in genes:
					row.append(sampledict[i][j])
			results.write('\t'.join(row) + '\n')

	with open(os.path.join(outdir, 'merged_results', 'mergedResults.transposed.tab'), 'w') as results:
		samples = sampledict.keys()
		genes = sampledict[samples[0]].keys()

		header = ['#Genes'] + samples
		results.write('\t'.join(header) + '\n')

		for j in genes:
			row = [j]
			for i in samples:
				row.append(sampledict[i][j])
			results.write('\t'.join(row) + '\n')

	if not os.path.exists(os.path.join(outdir, 'merged_results', 'consensus_sequences')):
		os.makedirs(os.path.join(outdir, 'merged_results', 'consensus_sequences'))

	for x in consensusdict:
		with open(os.path.join(outdir, 'merged_results', 'consensus_sequences', x + '_merged_sequences.fasta'), 'w') as sequenceResults:
			for z in consensusdict[x]:
				sequenceResults.write(z[0] + '\n')
				sequenceResults.write(z[1] + '\n')

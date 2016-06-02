import os
import csv
import sys
import time


def mergeResults(workdir, sequenceCoverage, outdir):

	#sampleList = os.listdir(workdir)
	dirs = [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d))]
	sampledict = {}
	consensusdict = {}
	prevNameSeq = ''
	countdirs = 0

	for i in dirs:
		countdirs+=1
		time.sleep(0.1)
		sys.stdout.write("\rChecking " + str(countdirs) + " of " + str(len(dirs)) + " folders...")
		sys.stdout.flush()
		sampleName = os.path.basename(i)
		mappingFilePath = os.path.join(workdir, sampleName, 'rematch_results', sampleName+'_mappingCheck.tab')
		consensusFilePath = os.path.join(workdir, sampleName, 'rematch_results', sampleName+'_sequences.fasta')
		if os.path.isfile(mappingFilePath): #parse mappingCheck file
			with open(mappingFilePath, 'r') as csvfile:
				mappingreader = csv.reader(csvfile, delimiter='\t')
				sampledict[sampleName] = {}
				for line in mappingreader:
					if not line[0].startswith("#"):
						if 1 - float(line[3]) >= sequenceCoverage:
							if line[6] == 'False':
								sampledict[sampleName][line[0]] = line[8]
							else:
								sampledict[sampleName][line[0]] = 'Mul_Allele (' + line[8] + ')'
						else:
							sampledict[sampleName][line[0]] = 'Absent (' + line[8] + ')'

		if os.path.isfile(consensusFilePath): #parse consensus sequences
			with open(consensusFilePath, 'r') as consensusFile:
				countSequences = -1
				for line in consensusFile:
					line = line.splitlines()[0]
					# if '>' in line:
					if line.startswith('>'):
						nameseq = line[1:]
						prevNameSeq = nameseq
						nameseq = nameseq.replace(' ', '_')
						if nameseq not in consensusdict:
							consensusdict[nameseq] = []
							# consensusdict[nameseq].append(['>' + sampleName + '_' + nameseq, False])
							consensusdict[nameseq].append(['>' + sampleName + '-' + prevNameSeq, ''])
						else:
							# consensusdict[nameseq].append(['>' + sampleName + '_' + nameseq, False])
							consensusdict[nameseq].append(['>' + sampleName + '-' + prevNameSeq, ''])
						# prevNameSeq = nameseq
						countSequences=len(consensusdict[nameseq])-1
					else:
						#print consensusdict
						# consensusdict[prevNameSeq][countSequences][1] = line 
						consensusdict[nameseq][countSequences][1] = consensusdict[nameseq][countSequences][1] + line
						#print consensusdict[prevNameSeq][countSequences]
	print "\nWriting results..."
	
	if not os.path.exists(os.path.join(outdir, 'merged_results')):
		os.makedirs(os.path.join(outdir, 'merged_results'))

	
	with open(os.path.join(outdir, 'merged_results', 'mergedResults.tab'),'w') as results:
		firstLine = True
		header = 'Samples\t'

		for x in sampledict:
			if firstLine:
				for l in sampledict[x]:
					header += l + '\t'
				results.write(header + '\n')
				firstLine = False

			row = x + '\t'
			for k in sampledict[x]:
				row += sampledict[x][k] + '\t'
			results.write(row + '\n')

	if not os.path.exists(os.path.join(outdir, 'merged_results', 'consensus_sequences')):
		os.makedirs(os.path.join(outdir, 'merged_results', 'consensus_sequences'))

	
	for x in consensusdict:
		# with open(os.path.join(workdir, 'merged_results', 'consensus_sequences', x.strip('\n').strip(' ') + '_merged_sequences.fasta'),'w') as sequenceResults:
		with open(os.path.join(outdir, 'merged_results', 'consensus_sequences', x + '_merged_sequences.fasta'),'w') as sequenceResults:
			for z in consensusdict[x]:
				sequenceResults.write(z[0] + '\n')
				sequenceResults.write(z[1] + '\n')

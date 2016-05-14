import os
import csv
import sys
import time


def mergeResults(workdir, sequenceCoverage):

	#sampleList = os.listdir(workdir)
	dirs = [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d))]
	sampledict = {}
	consensusdict = {}
	prevNameSeq = ''
	countdirs = 0

	for i in dirs:
		countdirs+=1
		time.sleep(1)
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
								sampledict[sampleName][line[0]] = 'Mul.Allele'
						else:
							sampledict[sampleName][line[0]] = 'Absent'

		if os.path.isfile(consensusFilePath): #parse consensus sequences
			with open(consensusFilePath, 'r') as consensusFile:
				countSequences = -1
				for line in consensusFile:
					if '>' in line:
						nameseq = line[1:]
						if nameseq not in consensusdict:
							consensusdict[nameseq] = []
							consensusdict[nameseq].append(['>' + sampleName + '_' + nameseq, False])
						else:
							consensusdict[nameseq].append(['>' + sampleName + '_' + nameseq, False])
						prevNameSeq = nameseq
						countSequences=len(consensusdict[nameseq])-1
					else:
						#print consensusdict
						consensusdict[prevNameSeq][countSequences][1] = line 
						#print consensusdict[prevNameSeq][countSequences]
	print "\nWriting results..."
	
	if not os.path.exists(os.path.join(workdir, 'merged_results')):
		os.makedirs(os.path.join(workdir, 'merged_results'))

	
	with open(os.path.join(workdir, 'merged_results', 'mergedResults.tab'),'w') as results:
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

	if not os.path.exists(os.path.join(workdir, 'merged_results', 'consensus_sequences')):
		os.makedirs(os.path.join(workdir, 'merged_results', 'consensus_sequences'))

	
	for x in consensusdict:
		with open(os.path.join(workdir, 'merged_results', 'consensus_sequences', x.strip('\n').strip(' ') + '_merged_sequences.fasta'),'w') as sequenceResults:
			for z in consensusdict[x]:
				sequenceResults.write(z[0].strip('\n').strip(' ') + '\n')
				sequenceResults.write(z[1].strip('\n').strip(' ') + '\n')

	
						






	




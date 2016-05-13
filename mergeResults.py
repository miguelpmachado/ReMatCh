import os
import csv


def mergeResults(workdir):

	sequenceCoverage = 0.8
	#sampleList = os.listdir(workdir)
	dirs = [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d))]
	sampledict = {}
	print dirs

	for i in dirs:
		sampleName = os.path.basename(i)
		print os.path.join(workdir, sampleName, sampleName+'_mappingCheck.tab')
		mappingFilePath = os.path.join(workdir, sampleName, 'rematch_results', sampleName+'_mappingCheck.tab')
		if os.path.isfile(mappingFilePath):
			with open(mappingFilePath, 'r') as csvfile:
				mappingreader = csv.reader(csvfile, delimiter='\t')
				sampledict[sampleName] = {}
				for line in mappingreader:
					if not line[0].startswith("#"):
						if 1 - float(line[3]) >= sequenceCoverage:
							if line[6] == False:
								sampledict[sampleName][line[0]] = line[8]
							else:
								sampledict[sampleName][line[0]] = 'Mul.Allele'
						else:
							sampledict[sampleName][line[0]] = 'Absent'

	with open(os.path.join(workdir,'mergedResults.tab'),'w') as results:
		firstLine = True
		header = 'Samples\t'

		for x in sampledict:
			if firstLine:
				for l in sampledict[x]:
					header += sampledict[x][l] + '\t'
				results.write(header + '\n')
				firstLine = False

			row = x + '\t'
			for k in sampledict[x]:
				row += sampledict[x][k] + '\t'
			results.write(row + '\n')
						






	




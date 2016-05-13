import os
import csv


def mergeResults(workdir):

	#sampleList = os.listdir(workdir)
	dirs = [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d))]
	sampledict = {}
	print dirs

	for i in dirs:
		sampleName = os.path.basename(i)
		print os.path.join(workdir, sampleName, sampleName+'_mappingCheck.tab')
		mappingFilePath = os.path.join(workdir, sampleName, 'rematch_results', sampleName+'_mappingCheck.tab')
		if os.path.isfile(mappingFilePath):
			with open(mappingFilePath, 'rb') as csvfile:
				mappingreader = csv.reader(csvfile, delimiter='\t', quotechar='#')
				sampledict[sampleName] = {}
				for j in mappingreader:
					print j

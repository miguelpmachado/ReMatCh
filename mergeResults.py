import os
import numpy


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
			sampledict[sampleName] = {}
			mappingCheckFile = numpy.loadtxt(mappingFilePath, dtype={'formats': (numpy.str_, numpy.float, numpy.float, numpy.float, numpy.str_, numpy.str_, numpy.str_, numpy.float, numpy.float)}, comments='#', delimiter='\t', converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)
			print mappingCheckFile.shape
			print mappingCheckFile[1][:]
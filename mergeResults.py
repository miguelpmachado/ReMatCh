import os
import numpy as np


def mergeResults(workdir):

	sampleList = os.listdir(workdir)
	sampledict = {}

	for i in sampleList:
		sampleName = os.path.basename(i)
		if os.path.isfile(os.path.join(workdir, sampleName, sampleName+'_mappingCheck.tab')):
			sampledict[sampleName] = {}
			mappingCheckFile = numpy.loadtxt(coveragefile, comments='#', delimiter='\t', converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)
			print mappingCheckFile.shape
			print mappingCheckFile[1][:]
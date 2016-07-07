#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
ReMatCh
<https://github.com/bfrgoncalves/ReMatCh>

Copyright (C) 2016 Bruno Goncalves <bfgoncalves@medicina.ulisboa.pt>

Last modified: July 04, 2016

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import shutil
from datetime import datetime
import time

import downloadAndBowtie
from scriptAfterMapping import checkCoverage
from scriptAfterMapping import convertToBAM
from scriptAfterMapping import rawCoverage
from scriptAfterMapping import alleleCalling
from SeqFromWebTaxon import GetSequencesFromTaxon
import rematch_utils


def main():
	start_time = time.time()

	version = '1.1'

	args = rematch_utils.parseArguments(version)

	args.func(args, version)

	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken, 3600)
	minutes, seconds = divmod(rest, 60)
	print ">>> Runtime :" + str(hours) + "h:" + str(minutes) + "m:" + str(round(seconds, 2)) + "s" + "\n"


def runReMatCh(args, version):
	workdir = args.workdir[0]

	# Create working directory
	rematch_utils.check_create_directory(workdir)

	# Start logger
	logfile = rematch_utils.start_logger(workdir)

	# Get general information
	rematch_utils.general_information(logfile, version)

	rematch_utils.checkPrograms(args)

	toClear = []

	if args.tax[0]:
		GetSequencesFromTaxon(args.tax[0], args.l, True, True)

	if not args.allplat:
		platform = "Illumina"
	else:
		platform = None

	ids_with_problems = open(os.path.join(workdir, 'ids_with_problems.txt'), 'w')
	ids_no_problems = open(os.path.join(workdir, 'ids_no_problems.txt'), 'w')

	with open(args.l, 'r') as run_ids:
		reference_file = args.reference[0].name

		count_runs = 0
		buildBowtie = True
		firstLine = True
		for run_id_line in run_ids:
			toClear = []
			run = False

			if args.tax[0] is not None and firstLine is True:
				firstLine = False
				continue
			elif args.tax[0] is not None and platform is not None:
					run_info = run_id_line.split("\t")
					run_id = run_info[0]
					run_plat = run_info[1]

					if platform in run_plat and 'Analyzer' not in run_plat:
						run = True
			elif args.tax[0]:
				run_info = run_id_line.split("\t")
				run_id = run_info[0]
				run_plat = run_info[1]

				run = True
			else:
				run = True

			if args.tax[0] is not None and run is True:
				if args.useOmicsDataType[0] != 'All':
					run_info = run_id_line.split("\t")
					omics = run_info[3]
					if omics not in args.useOmicsDataType:
						run = False

			if run:

				startTime = datetime.now()

				count_runs += 1

				if count_runs > 1 or not args.bowtieBuild:
					buildBowtie = False

				run_id = run_id.strip()

				print "\nRunning ID: " + str(run_id)
				samFilePath, singOrPaired, filesDownloaded = downloadAndBowtie.downloadAndBowtie(reference_file, run_id, workdir, buildBowtie, args.picard[0], args.threads, toClear, args.asperaKey[0])

				if samFilePath is not False:

					sortedPath = convertToBAM(samFilePath, toClear)

					rawCoverage(sortedPath, toClear)
					print "Checking coverage..."
					sequenceNames, sequenceMedObject = checkCoverage(sortedPath, args.minCoverage, args.xtraSeq, toClear)
					print "Performing Allele Call..."
					alleleCalling(sortedPath, reference_file, sequenceNames, args.gatk[0], run_id, args.minQuality, args.minCoverage, args.multipleAlleles, sequenceMedObject, args.xtraSeq, toClear)
					print str(run_id) + " DONE"

					gzSizes = 0

					for files in filesDownloaded:
						gzSizes += float(os.path.getsize(os.path.join(workdir, run_id, 'fastq', files)))

					if args.rmFastq:
						shutil.rmtree(os.path.join(workdir, run_id, 'fastq'))

					ids_no_problems.write(run_id + "\n")
					ids_no_problems.flush()

					run_time = str(datetime.now() - startTime)

					with open(os.path.join(workdir, run_id, run_id + '_runtime.txt'), 'w') as runTimeFile:
						runTimeFile.write("#runTime\tfileSize\tlibraryLayout\n")
						runTimeFile.write(str(run_time) + '\t' + str(gzSizes) + "\t" + singOrPaired + '\n')
				else:
					ids_with_problems.write(run_id + "\n")
					ids_with_problems.flush()
					if args.rmFastq:
						try:
							shutil.rmtree(os.path.join(workdir, run_id, 'fastq'))
						except Exception as e:
							print e
					print run_id + ' - An error has occurred.'

			if args.clean:
				rematch_utils.removeFromArray(toClear)

		ids_with_problems.close()
		ids_no_problems.close()

		if args.clean:
			rematch_utils.removeIndexes(reference_file)


if __name__ == "__main__":
	main()

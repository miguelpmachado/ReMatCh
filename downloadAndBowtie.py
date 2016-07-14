#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2015 msilva <msilva@msilva>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import ftplib
import os, os.path, glob
import shutil

import rematch_utils

filesExtensions = ['fastq.gz', 'fq.gz']
pairEnd_file = ['_1.f', '_2.f']


def ftpListFiles(ftp, link):
	ftp.cwd(link)
	dirs = ftp.nlst()

	files = []

	for item in dirs:
		files.append(item)

	if len(files) == 0:
		files = None

	return files


def ftpSearchFileTypes(files):
	files_to_download = []
	if files is not None:
		if len(files) == 1:
			for extensions in filesExtensions:
				if extensions in files[0]:
					files_to_download.append(files[0])
					break
		else:
			for file_ena in files:
				for extensions in filesExtensions:
					if extensions in file_ena:
						for pairEnd_file_number in pairEnd_file:
							if pairEnd_file_number in file_ena:
								files_to_download.append(file_ena)
								break
						break
	if len(files_to_download) == 0:
		files_to_download = None

	return files_to_download


def getFilesList(runID):
	run_successfully = False

	partial_tid = runID[0:6]

	files = None

	try:
		f = ftplib.FTP('ftp.sra.ebi.ac.uk', timeout=3600)
		f.login()

		link = '/vol1/fastq/' + partial_tid + '/' + runID

		try:
			files = ftpListFiles(f, link)
			run_successfully = True
		except Exception as e:
			print link
			print e

			link = '/vol1/fastq/' + partial_tid + "/00" + runID[-1] + '/' + runID
			try:
				files = ftpListFiles(f, link)
				run_successfully = True
			except Exception as e:
				print link
				print e

		try:
			f.quit()
		except Exception as e:
			print e
	except Exception as e:
		print e

	print files
	files = ftpSearchFileTypes(files)
	print files

	return run_successfully, files


def download(dirs2, target_dir2, ref2, success2, f2, link2):
	insucess = 0

	files = ftpListFiles(f2, link2)
	files = ftpSearchFileTypes(files)

	for item in files:

		try:
			f2.cwd(link2)
			final_target_dir = target_dir2 + "/" + item
			file = open(final_target_dir, 'wb')
			print "Downloading: %s" % item

			f2.retrbinary('RETR %s' % item, file.write)
			file.close()
			print "Downloaded %s" % item
			success2 += 1
		except Exception as e:
			print e
			insucess += 1

	return success2, insucess


def download_ERR(ERR_id, target_dir):
	ref = ERR_id
	failed = 0
	success = 0
	insucess = 0

	try:
		f = ftplib.FTP('ftp.sra.ebi.ac.uk', timeout=3600)
		f.login()

		try:

			firstid = ref[0:6]
			# get the read files name from the reference id
			link = '/vol1/fastq/' + firstid + "/" + ref
			f.cwd(link)
			dirs = f.nlst()

		except:
			try:
				firstid = ref[0:6]
				# get the read files name from the reference id
				link = '/vol1/fastq/' + firstid + "/00" + ref[-1] + "/" + ref
				f.cwd(link)
				dirs = f.nlst()
			except Exception as e:
				failed += 1
				print "Bad ID: " + ref
			else:
				success, insucess = download(dirs, target_dir, ref, success, f, link)
		else:
			success, insucess = download(dirs, target_dir, ref, success, f, link)
		try:
			f.quit()
		except Exception as e:
			print e
			print 'Insucess: ' + str(insucess)
	except Exception as e:
		print e

	print "Downloaded %s files successfully, %s fail and %s ID references were wrong" % (success, insucess, failed)
	return insucess


def aspera(run_id, asperaKey, outdir, fileToDownload):
	if fileToDownload is None:
		fileToDownload = ''
	else:
		fileToDownload = '/' + fileToDownload

	aspera_command = ['ascp', '-QT', '-l', '300m', '-i', asperaKey, '', outdir]
	aspera_command[6] = str('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/' + run_id[0:6] + '/' + run_id + fileToDownload)
	run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(aspera_command, False, 3600)
	if not run_successfully:
		print 'It was not possible to download! Trying again:'
		aspera_command[6] = str('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/' + run_id[0:6] + '/00' + run_id[-1] + '/' + run_id + fileToDownload)
		run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(aspera_command, False, 3600)

	return run_successfully


# Download using Aspera Connect
def downloadAspera(run_id, outdir, asperaKey, getAllFiles_Boolean, filesToDownload):
	run_successfully = False

	if getAllFiles_Boolean:
		run_successfully = aspera(run_id, asperaKey, outdir, None)
		if run_successfully:
			files = glob.glob1(os.path.join(outdir, run_id), '*')
			for file_downloaded in files:
				shutil.move(os.path.join(outdir, run_id, file_downloaded), outdir)
			shutil.rmtree(os.path.join(outdir, run_id))
	else:
		if filesToDownload is not None:
			runs = []
			print filesToDownload
			for file_ena in filesToDownload:
				print file_ena
				run_successfully = aspera(run_id, asperaKey, outdir, file_ena)
				runs.append(run_successfully)

			if False in runs:
				run_successfully = False
		else:
			run_successfully = True

	return run_successfully


# Search Fastq files (that were downloaded or already provided by the user)
def searchDownloadedFiles(directory):
	for extension in filesExtensions:
		downloadedFiles = glob.glob1(directory, str('*.' + extension))
		if len(downloadedFiles) > 0:
			break
	return downloadedFiles


def downloadAndBowtie(referencePath, run_id, target_dir, buildBowtie, picardJarPath, threads, toClear, asperaKey, removeFastq):

	bowtieBuildFileName, extension = os.path.splitext(referencePath)

	if buildBowtie:
		print "Indexing Reference file..."

		# Picard
		picardFileName, extension = os.path.splitext(referencePath)
		command = ['java', '-jar ', picardJarPath, ' CreateSequenceDictionary', str('R=' + referencePath), str('O=' + picardFileName + '.dict'), '2>', str(picardFileName + '_picard_out.txt')]
		run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(command, True, None)

		# Bowtie build
		command = ["bowtie2-build", referencePath, bowtieBuildFileName]
		run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(command, False, None)

	if not os.path.isdir(target_dir):
		os.makedirs(target_dir)
	dir_sample = os.path.join(target_dir, run_id)
	if not os.path.isdir(dir_sample):
		os.makedirs(dir_sample)

	dir_with_gz = os.path.join(target_dir, run_id, 'fastq')
	if not os.path.isdir(dir_with_gz):
		os.makedirs(dir_with_gz)

	# download ERR
	aspera_run = False
	ftp_down_insuc = 0

	downloadedFiles = searchDownloadedFiles(dir_with_gz)
	download_step_performed = False
	if len(downloadedFiles) < 1:
		print 'Trying download...'
		if asperaKey is not None:
			run_successfully, files = getFilesList(run_id)
			print files
			if run_successfully:
				aspera_run = downloadAspera(run_id, dir_with_gz, asperaKey, False, files)
			else:
				aspera_run = downloadAspera(run_id, dir_with_gz, asperaKey, True, None)
			if not aspera_run:
				print 'Trying download using FTP'
				ftp_down_insuc = download_ERR(run_id, dir_with_gz)
		else:
			ftp_down_insuc = download_ERR(run_id, dir_with_gz)

		download_step_performed = True

	else:
		print 'File ' + run_id + ' already exists...'

	downloadedFiles = searchDownloadedFiles(dir_with_gz)

	if ftp_down_insuc > 0 and aspera_run is False:
		shutil.rmtree(dir_with_gz)
		return False, False, downloadedFiles

	if len(downloadedFiles) > 2:
		filesToUse = []
		for i in downloadedFiles:
			removeFile = True
			for pairEnd_file_number in pairEnd_file:
				if pairEnd_file_number in i:
					filesToUse.append(i)
					removeFile = False
					break
			if removeFile and download_step_performed:
				shutil.rmtree(os.path.join(dir_with_gz, i))
		downloadedFiles = filesToUse
		print 'Files used: ' + str(downloadedFiles)

	resultsFolder = os.path.join(dir_sample, 'rematch_results')

	if not os.path.exists(resultsFolder):
		os.makedirs(resultsFolder)

	bowtie_output_file = os.path.join(resultsFolder, str(run_id + '.sam'))

	pairedOrSingle = "Single_end"

	print "Running bowtie..."

	bowtie_command = ['bowtie2', '-k', '2', '--quiet', '--no-unal', '-x', bowtieBuildFileName, str(''), str('--rg-id ENA --rg SM:' + run_id), '--sensitive-local', '--threads', str(threads), '--met-file', os.path.join(resultsFolder, str(run_id + '.bowtie_metrics.txt')), '-S', bowtie_output_file]

	if len(downloadedFiles) == 1:
		bowtie_command[7] = str('-U ' + os.path.join(dir_with_gz, downloadedFiles[0]))
		run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(bowtie_command, False, None)
		if not run_successfully:
			print 'Bowtie2 fails!'
			print stdout.decode("utf-8"), stderr.decode("utf-8")
			return False, False, downloadedFiles
	elif len(downloadedFiles) == 2:
		pairedOrSingle = "Paired_end"
		bowtie_command[7] = str('-1 ' + os.path.join(dir_with_gz, downloadedFiles[0]) + ' -2 ' + os.path.join(dir_with_gz, downloadedFiles[1]))
		run_successfully, stdout, stderr = rematch_utils.runCommandPopenCommunicate(bowtie_command, False, None)
		if not run_successfully:
			print 'Bowtie2 fails!'
			print stdout.decode("utf-8"), stderr.decode("utf-8")
			return False, False, downloadedFiles
	else:
		print "0 fastq files found. Aborting..."
		shutil.rmtree(dir_with_gz)
		return False, False, downloadedFiles

	toClear.append(os.path.join(resultsFolder, run_id + ".bowtie_metrics.txt"))
	return bowtie_output_file, pairedOrSingle, downloadedFiles

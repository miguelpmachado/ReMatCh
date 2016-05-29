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

import shlex, subprocess, ftplib
import os, os.path, glob
import shutil

from os import listdir
from os.path import isfile, join
from threading import Timer

def runTimeoutLimit(command_fragmented_in_list, timeout_sec):
	proc = subprocess.Popen(command_fragmented_in_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	timer = Timer(timeout_sec, proc.kill)
	timer.start()
	stdout,stderr = proc.communicate()
	timer.cancel()
	return proc.returncode, stdout.decode("utf-8"), stderr.decode("utf-8")

def download(dirs2,target_dir2,ref2,success2,f2,link2):
	insucess=0
	for item in dirs2:
					
		try:
			f2.cwd(link2)
			final_target_dir=target_dir2+"/"+ item
			file = open(final_target_dir, 'wb')
			print "Downloading: %s" % item

			f2.retrbinary('RETR %s' % item, file.write)
			file.close()
			print "Downloaded %s" % item
			success2+=1
		except Exception as e:
			print e
			insucess+=1

	return success2,insucess

def download_ERR(ERR_id,target_dir):
	
	f=ftplib.FTP('ftp.sra.ebi.ac.uk', timeout=3600)
	f.login()
	#for each reference id get the read files
	ref=ERR_id
	failed=0
	success=0
	insucess=0

	try:
				
		firstid=ref[0:6]
		#get the read files name from the reference id
		link='/vol1/fastq/'+firstid+"/"+ref
		f.cwd(link)
		dirs=f.nlst();
				
	except:
		try:
			firstid=ref[0:6]
			#get the read files name from the reference id
			link='/vol1/fastq/'+firstid+"/00"+ref[-1]+"/"+ref
			f.cwd(link)
			dirs=f.nlst();
					
					
						
		except Exception as e:
			failed +=1
			print "Bad ID: " + ref
		else:
			success,insucess=download(dirs,target_dir,ref,success,f,link)	
			
	else:	
		success,insucess=download(dirs,target_dir,ref,success,f,link)
	
	try:
		f.quit()
	except Exception as e:
		print e
		print 'Insucess: ' + str(insucess)
	print "Downloaded %s files successfully, %s fail and %s ID references were wrong" % (success,insucess,failed)
	return insucess

# Download using Aspera Connect
def downloadAspera(run_id, outdir, asperaKey):
	aspera_run = False
	aspera_command = ['ascp', '-QT', '-l', '300m', '-i', asperaKey, str('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/' + run_id[0:6] + '/' + run_id), outdir]
	aspera, std_out, std_err = runTimeoutLimit(aspera_command, 3600)
	if aspera == 0:
		aspera_run = True
		print ' '.join(aspera_command)
	else:
		print "Connecting using the following command didn't work:"
		print ' '.join(aspera_command)
		aspera_command[6] = str('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/' + run_id[0:6] + '/00' + run_id[-1] + '/' + run_id)
		print 'Trying: ' + ' '.join(aspera_command)
		aspera, std_out, std_err = runTimeoutLimit(aspera_command, 3600)
		if aspera == 0:
			aspera_run = True
		else:
			print "Download using Aspera didn't work"
			return aspera_run
	files = glob.glob1(os.path.join(outdir, run_id), '*')
	for file in files:
		shutil.move(os.path.join(outdir, run_id, file), outdir)
	shutil.rmtree(os.path.join(outdir, run_id))
	return aspera_run

# Search Fastq files (that were downloaded or already provided by the user)
def searchDownloadedFiles(directory):
	filesExtensions = ['fastq.gz', 'fq.gz']
	for extension in filesExtensions:
		downloadedFiles = glob.glob1(directory, str('*.' + extension))
		if len(downloadedFiles) > 0:
			break
	return downloadedFiles

def downloadAndBowtie(referencePath, run_id, target_dir, buildBowtie, picardJarPath, threads, toClear, asperaKey):

	bowtieBuildFileName, extension = os.path.splitext(referencePath)

	if buildBowtie == True:
		print "Indexing Reference file..."
		picardFileName, extension = os.path.splitext(referencePath)
		os.system('java -jar ' + picardJarPath +" CreateSequenceDictionary R= " + referencePath + " O= " + picardFileName + ".dict 2> "+picardFileName+"_picard_out.txt")

		bowtiBuildeLog=bowtieBuildFileName +"_bowtiBuildLog.txt"
		myoutput = open(bowtiBuildeLog, 'w')
		subprocess.call(["bowtie2-build", referencePath, bowtieBuildFileName],stdout=myoutput,stderr=myoutput)
	
	if not os.path.isdir(target_dir):
		os.makedirs(target_dir)
	dir_sample = os.path.join(target_dir,run_id)
	if not os.path.isdir(dir_sample):
		os.makedirs(dir_sample)

	dir_with_gz = os.path.join(target_dir,run_id, 'fastq')
	if not os.path.isdir(dir_with_gz):
		os.makedirs(dir_with_gz)
	
	#download ERR
	aspera_run = False
	ftp_down_insuc = 0
	
	downloadedFiles = searchDownloadedFiles(dir_with_gz)
	if len(downloadedFiles) < 1:
		print 'Trying download...'
		if asperaKey != None:
			aspera_run = downloadAspera(run_id, dir_with_gz, asperaKey[0])
			if aspera_run == False:
				print 'Trying download using FTP'
				ftp_down_insuc=download_ERR(run_id, dir_with_gz)
		else:
			ftp_down_insuc=download_ERR(run_id, dir_with_gz)
	else:
		print 'File '+ run_id + ' already exists...' 
	
	downloadedFiles = searchDownloadedFiles(dir_with_gz)
	
	if ftp_down_insuc > 0 and aspera_run == False:
		os.system('rm -r ' + dir_with_gz)
		return False, False, downloadedFiles, False
		
	
	if len(downloadedFiles) > 2:
		filesToUse = []
		for i in downloadedFiles:
			if '_1.f' in i or '_2.f' in i:
				filesToUse.append(i)
			else:
				os.system('rm ' + os.path.join(dir_with_gz, i))
		downloadedFiles = filesToUse
		print 'Files used: ' + str(downloadedFiles)
	
	

	resultsFolder = os.path.join(dir_sample, 'rematch_results')
	
	if not os.path.exists(resultsFolder):
		os.makedirs(resultsFolder)
	
	bowtie_output_file=os.path.join(resultsFolder, str(run_id + '.sam'))
	bowtieLog = os.path.join(resultsFolder, run_id + "_bowtie_error.txt")
	
	pairedOrSingle="Single_end"	

	print "Running bowtie..."
	
	bowtie_command = ['bowtie2', '-k', '2', '--quiet', '--no-unal', '-x', bowtieBuildFileName, str(''), str('--rg-id ENA --rg SM:' + run_id), '--sensitive-local', '--threads', str(threads), '--met-file', os.path.join(resultsFolder, str(run_id + '.bowtie_metrics.txt')), '-S', bowtie_output_file]
	
	print os.listdir(dir_with_gz)
	
	if len(downloadedFiles)==1:
		bowtie_command[7] = str('-U ' + os.path.join(dir_with_gz, downloadedFiles[0]))
		myoutput = open(bowtieLog, 'w')
		proc = subprocess.call(bowtie_command, stdout=myoutput, stderr=myoutput)
		myoutput.close()
		if proc != 0:
			print 'Bowtie2 fails! Find Bowtie2 log file at' + bowtieLog
			return False, False, downloadedFiles
	elif len(downloadedFiles)==2:
		bowtie_command[7] = str('-1 ' + os.path.join(dir_with_gz, downloadedFiles[0]) + ' -2 ' + os.path.join(dir_with_gz, downloadedFiles[1]))
		myoutput = open(bowtieLog, 'w')
		proc = subprocess.call(bowtie_command, stdout=myoutput, stderr=myoutput)
		myoutput.close()
		if proc != 0:
			print 'Bowtie2 fails! Find Bowtie2 log file at' + bowtieLog
			return False, False, downloadedFiles
	else:
		print "0 fastq files found. Aborting..."
		os.system('rm -r ' + dir_with_gz)
		return False, False, downloadedFiles

	toClear.append(os.path.join(resultsFolder, run_id+"_bowtie_error.txt"))
	toClear.append(os.path.join(resultsFolder, run_id+".bowtie_metrics.txt"))
	return bowtie_output_file, pairedOrSingle, downloadedFiles

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

from os import listdir
from os.path import isfile, join
from threading import Timer ## mpmachado ##

def runTimeoutLimit(command_fragmented_in_list, timeout_sec): ## mpmachado ##
	proc = subprocess.Popen(command_fragmented_in_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	timer = Timer(timeout_sec, proc.kill)
	timer.start()
	stdout,stderr = proc.communicate()
	timer.cancel()
	return proc.returncode, stdout.decode("utf-8"), stderr.decode("utf-8")

def download(dirs2,target_dir2,ref2,success2,f2,link2):
	#new folder for each reference with reference id name
	# subprocess.call(['mkdir', target_dir2+"/"+ref2]) # mpmachado #
	
	#get fasta file for each read file name
	numFilesInDir = len(dirs2)
	insucess2=0

	if numFilesInDir > 2:
		print "more than 2 files"
		#logFile.write("more than 2 files" + '\n')
		return success2

	for item in dirs2:
					
		try:
			f2.cwd(link2)
			final_target_dir=target_dir2+"/"+ item
			file = open(final_target_dir, 'wb')
			print "Downloading: %s" % item
			#logFile.write("Downloading: %s\n" % item)

			f2.retrbinary('RETR %s' % item, file.write)
			file.close()
			print "Downloaded %s" % item
			#logFile.write("Downloaded %s\n" % item)
			success2+=1		
		except Exception as e:
			print e
			insucess2+=1

	print 'download:' + str(success2,insucess2)
	return success2,insucess2

def download_ERR(ERR_id,target_dir):
	
	# if not os.path.isdir(target_dir): # mpmachado #
		# os.makedirs(target_dir) # mpmachado #
	
	f=ftplib.FTP('ftp.sra.ebi.ac.uk', timeout=3600) # mpmachado #
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
					
					
						
		except Exception, e:
			failed +=1
			print "Bad ID: " + ref
			#logFile.write("Bad ID: " + ref + '\n')
		else:
			print 'except2:' + str(success,insucess)
			success,insucess=download(dirs,target_dir,ref,success,f,link)	
			
	else:
		print 'done:' + str(success,insucess)	
		success,insucess=download(dirs,target_dir,ref,success,f,link)
	
	f.quit()	
	print "Successfully downloaded %s files and %s ID references were wrong" % (success,failed)	
	#logFile.write("Successfully downloaded %s files and %s ID references were wrong\n" % (success,failed))
	return insucess

def downloadAspera(run_id, outdir, asperaKey): ## mpmachado ##
	aspera_run = False
	aspera_command = ['ascp', '-QT', '-l', '300m', '-i', asperaKey, str('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/' + run_id[0:6] + '/' + run_id), outdir]
	aspera, std_out, std_err = runTimeoutLimit(aspera_command, 3600)
	if aspera == 0:
		aspera_run = True
	else:
		print "Connecting using the following command didn't work:"
		print ' '.join(aspera_command)
		aspera_command[6] = str('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/' + run_id[0:6] + '/00' + run_id[-1] + '/' + run_id)
		print 'Trying: ' + ' '.join(aspera_command)
		aspera, std_out, std_err = runTimeoutLimit(aspera_command, 3600)
		if aspera == 0:
			aspera_run = True
		else:
			print "Download using Aspera didn't work" + "\n"
	return aspera_run

# Search Fastq files (that were downloaded or already provided by the user)
def searchDownloadedFiles(directory): ## mpmachado ##
	filesExtensions = ['fastq.gz', 'fq.gz']
	downloadedFiles = []
	for extension in filesExtensions:
		downloadedFiles = glob.glob1(directory, str('*.' + extension))
		if len(downloadedFiles) == 0:
			downloadedFiles = []
		else:
			break
	return downloadedFiles


def downloadAndBowtie(referencePath, run_id, target_dir, buildBowtie, picardJarPath, threads, toClear, asperaKey): # mpmachado #

	bowtieBuildFileName, extension = os.path.splitext(referencePath) # mpmachado #

	if buildBowtie == True:
		print "Running picard"
		#logFile.write("Running picard" + '\n')
		picardFileName, extension = os.path.splitext(referencePath)
		os.system(picardJarPath +" CreateSequenceDictionary R= " + referencePath + " O= " + picardFileName + ".dict 2> "+picardFileName+"_picard_out.txt")
		toClear.append(picardFileName+"_picard_out.txt")
		toClear.append(picardFileName + ".dict")
		toClear.append(referencePath+'.fai')
		
	# if buildBowtie == True: # mpmachado #
		print "Running bowtie..."
		#logFile.write("Running bowtie..." + '\n')
		bowtiBuildeLog=bowtieBuildFileName +"_bowtiBuildLog.txt"
		myoutput = open(bowtiBuildeLog, 'w')
		subprocess.call(["bowtie2-build", referencePath, bowtieBuildFileName],stdout=myoutput,stderr=myoutput)
		toClear.append(bowtieBuildFileName + ".*.bt2")
		toClear.append(bowtiBuildeLog)
	
	if not os.path.isdir(target_dir): # mpmachado #
		os.makedirs(target_dir) # mpmachado #
	dir_sample = os.path.join(target_dir,run_id)
	if not os.path.isdir(dir_sample): # mpmachado #
		os.makedirs(dir_sample) # mpmachado #

	dir_with_gz = os.path.join(target_dir,run_id, 'fastq')
	if not os.path.isdir(dir_with_gz): # mpmachado #
		os.makedirs(dir_with_gz) # mpmachado #
	
	#download ERR
	aspera_run = False
	ftp_down_insuc = 0
	
	# numberFilesDowned= len(glob.glob1(dir_with_gz, "*.fastq.gz")) # mpmachado #
	downloadedFiles = searchDownloadedFiles(dir_with_gz) # mpmachado #
	if len(downloadedFiles) < 1: # mpmachado #
		if asperaKey != None: ## mpmachado ##
			aspera_run = downloadAspera(run_id, dir_with_gz, asperaKey[0]) ## mpmachado ##
			if aspera_run == False:
				print 'Trying download using FTP' + "\n" ## mpmachado ##
				ftp_down_insuc=download_ERR(run_id, dir_with_gz) ## mpmachado ##
				
				
		else: ## mpmachado ##
			ftp_down_insuc=download_ERR(run_id, dir_with_gz) # mpmachado #
	else:
		print 'File '+ run_id+' already exists...' 
		#logFile.write('File '+ run_id+' already exists...' + '\n')
	
	downloadedFiles = searchDownloadedFiles(dir_with_gz) # mpmachado #
	
	if ftp_down_insuc>0 and aspera_run == False:
		return False, False, downloadedFiles # mpmachado #
		
	
	if len(downloadedFiles) > 2:
		filesToUse = []
		for i in downloadedFiles:
			if '_1.f' in i or '_2.f' in i:
				filesToUse.append(i)

		downloadedFiles = filesToUse
	
	

	resultsFolder = os.path.join(dir_sample, 'rematch_results')
	
	if not os.path.exists(resultsFolder):
		os.makedirs(resultsFolder)
	
	bowtie_output_file=os.path.join(resultsFolder, run_id + ".sam")
	bowtieLog = os.path.join(resultsFolder, run_id + "_bowtie_error.txt")
	
	pairedOrSingle="Single_end"	

	
	if len(downloadedFiles)==1:


		command_line ="bowtie2 -k 2 --quiet --no-unal -x "+bowtieBuildFileName+" -U "+dir_with_gz+"/"+downloadedFiles[0] + " --rg-id ENA --rg SM:"+run_id+" --sensitive-local --threads "+ str(threads) +" --met-file "+ os.path.join(resultsFolder, run_id+".bowtie_metrics.txt") + " -S "+bowtie_output_file+" "


		myoutput = open(bowtieLog, 'w')
		args = shlex.split(command_line)
		p = subprocess.call(args,stdout=myoutput,stderr=myoutput)
		toClear.append(os.path.join(resultsFolder, run_id+".bowtie_metrics.txt"))
		toClear.append(os.path.join(resultsFolder, run_id+"_bowtie_error.txt"))


	elif len(downloadedFiles)==2:

		command_line ="bowtie2 -k 2 --quiet --no-unal -x "+bowtieBuildFileName+" -1 "+dir_with_gz+"/"+downloadedFiles[0]+ " -2 "+dir_with_gz+"/"+downloadedFiles[1] + " --rg-id ENA --rg SM:"+run_id+" --sensitive-local --threads "+ str(threads) +" --met-file "+ os.path.join(resultsFolder, run_id+".bowtie_metrics.txt") + " -S "+bowtie_output_file+" "
		
		myoutput = open(bowtieLog, 'w')
		args = shlex.split(command_line)
		p = subprocess.call(args,stdout=myoutput,stderr=myoutput)
		pairedOrSingle="Paired_end"	
		toClear.append(os.path.join(resultsFolder, run_id+".bowtie_metrics.txt"))
		toClear.append(os.path.join(resultsFolder, run_id+"_bowtie_error.txt"))

	
	else:
		
		print "0 fastq files found. Aborting...\n"
		#logFile.write("0 fastq files found. Aborting..." + '\n')
		os.system('rm -r ' + dir_with_gz)
		return False, False, downloadedFiles # mpmachado #

	
	return bowtie_output_file, pairedOrSingle, downloadedFiles # mpmachado #



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

import shlex, subprocess,ftplib
import os, os.path, glob

from os import listdir
from os.path import isfile, join

def download(dirs2,target_dir2,ref2,success2,f2,link2, logFile):
	#new folder for each reference with reference id name
	subprocess.call(['mkdir', target_dir2+"/"+ref2])
				
	#get fasta file for each read file name
	numFilesInDir = len(dirs2)

	if numFilesInDir > 2:
		print "more than 2 files"
		logFile.write("more than 2 files" + '\n')
		return success2

	for item in dirs2:
					
		f2.cwd(link2)
		final_target_dir=target_dir2+"/"+ref2 +"/"+ item
		file = open(final_target_dir, 'wb')
		print "Downloading: %s" % item
		logFile.write("Downloading: %s\n" % item)

		f2.retrbinary('RETR %s' % item, file.write)
		file.close()
		print "Downloaded %s" % item
		logFile.write("Downloaded %s\n" % item)
		success2+=1		

	return success2

def download_ERR(ERR_id,target_dir, logFile):
	
	if not os.path.isdir(target_dir):
		os.makedirs(target_dir)
	
	f=ftplib.FTP('ftp.sra.ebi.ac.uk')
	f.login()
	#for each reference id get the read files
	ref=ERR_id
	failed=0
	success=0

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
			logFile.write("Bad ID: " + ref + '\n')
		else:
			success=download(dirs,target_dir,ref,success,f,link, logFile)	
			
	else:
				
		success=download(dirs,target_dir,ref,success,f,link)
	
	f.quit()	
	print "Successfully downloaded %s files and %s ID references were wrong" % (success,failed)	
	logFile.write("Successfully downloaded %s files and %s ID references were wrong\n" % (success,failed))



def downloadAndBowtie(referencePath, run_id, target_dir, buildBowtie, picardJarPath, threads, logFile):

	if buildBowtie == True:
		print "run picard"
		logFile.write("run picard" + '\n')
		picardFileName, extension = os.path.splitext(referencePath)
		os.system("java -jar "+ picardJarPath +" CreateSequenceDictionary R= " + referencePath + " O= " + picardFileName + ".dict 2> "+picardFileName+"_picard_out.txt")



	bowtieBuildFileName, extension = os.path.splitext(referencePath)


	if buildBowtie == True:
		print "run bowtie"
		logFile.write("run bowtie" + '\n')
		bowtiBuildeLog=bowtieBuildFileName+"_bowtiBuildeLog.txt"
		myoutput = open(bowtiBuildeLog, 'w')
		subprocess.call(["bowtie2-build", referencePath, bowtieBuildFileName],stdout=myoutput,stderr=myoutput)
	
	#download ERR

	dir_with_gz = os.path.join(target_dir,run_id)

	numberFilesDowned= len(glob.glob1(dir_with_gz, "*.fastq.gz")) 

	if numberFilesDowned < 1:
		download_ERR(run_id, target_dir, logFile)
	else:
		print 'File '+ run_id+' already exists...' 
		logFile.write('File '+ run_id+' already exists...' + '\n')
	
	#download_ERR(run_id, target_dir)

	numberFilesDowned = len(glob.glob1(dir_with_gz, "*.fastq.gz")) 

	#print len(glob.glob1(dir_with_gz, "*.fastq.gz")) 
	
	bowtie_output_file=os.path.join(dir_with_gz, run_id + ".sam")
	
	bowtieLog = os.path.join(dir_with_gz, run_id + "_bowtie_error.txt")
	
	pairedOrSingle="Single_end"	

	
	if numberFilesDowned==1:

		command_line ="bowtie2 -k 2 --quiet --no-unal -x "+bowtieBuildFileName+" -U "+dir_with_gz+"/"+run_id+".fastq.gz --rg-id ENA --rg SM:"+run_id+" --sensitive-local --threads "+ str(threads) +" --met-file "+ os.path.join(dir_with_gz, run_id+".bowtie_metrics.txt") + " -S "+bowtie_output_file+" "


		myoutput = open(bowtieLog, 'w')
		args = shlex.split(command_line)
		p = subprocess.call(args,stdout=myoutput,stderr=myoutput)


	elif numberFilesDowned==2:
		command_line ="bowtie2 -k 2 --quiet --no-unal -x "+bowtieBuildFileName+" -1 "+dir_with_gz+"/"+run_id+"_1.fastq.gz -2 "+dir_with_gz+"/"+run_id+"_2.fastq.gz --rg-id ENA --rg SM:"+run_id+" --sensitive-local --threads "+ str(threads) +" --met-file "+ os.path.join(dir_with_gz, run_id+".bowtie_metrics.txt") + " -S "+bowtie_output_file+" "
		
		myoutput = open(bowtieLog, 'w')
		args = shlex.split(command_line)
		p = subprocess.call(args,stdout=myoutput,stderr=myoutput)
		pairedOrSingle="Paired_end"	


	
	else:
		
		print "0 or more than 2 fastq files"
		logFile.write("0 or more than 2 fastq files" + '\n')
		os.rmdir(dir_with_gz)
		return False, False, numberFilesDowned

	
	return bowtie_output_file, pairedOrSingle, numberFilesDowned



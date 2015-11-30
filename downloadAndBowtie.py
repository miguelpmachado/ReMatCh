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

def download(dirs2,target_dir2,ref2,success2,f2,link2):
	#new folder for each reference with reference id name
	subprocess.call(['mkdir', target_dir2+"/"+ref2])
				
	#get fasta file for each read file name
	for item in dirs2:
					
		f2.cwd(link2)
		final_target_dir=target_dir2+"/"+ref2 +"/"+ item
		file = open(final_target_dir, 'wb')
		print "Downloading: %s" % item

		f2.retrbinary('RETR %s' % item, file.write)
		file.close()
		print "Downloaded %s" % item
		success2+=1		

	return success2

def download_ERR(ERR_id,target_dir):
	
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
		else:
			success=download(dirs,target_dir,ref,success,f,link)	
			
	else:
				
		success=download(dirs,target_dir,ref,success,f,link)
	
	f.quit()	
	print "Successfully downloaded %s files and %s ID references were wrong" % (success,failed)	



def downloadAndBowtie(referencePath, run_id, target_dir, buildBowtie, picardJarPath, threads):

	print "run picard"
	if buildBowtie == True:
		picardFileName, extension = os.path.splitext(referencePath)
		os.system("java -jar "+ picardJarPath +" CreateSequenceDictionary R= " + referencePath + " O= " + picardFileName + ".dict")



	bowtieBuildFileName, extension = os.path.splitext(referencePath)
	
	print "run bowtie"

	if buildBowtie == True:
		subprocess.call(["bowtie2-build", referencePath, bowtieBuildFileName])
	
	#download ERR

	if not os.path.isdir(os.path.join(target_dir, run_id)):
		download_ERR(run_id, target_dir)
	else:
		print 'File already exists...'
	
	#download_ERR(run_id, target_dir)
	
	dir_with_gz = os.path.join(target_dir,run_id)

	#print len(glob.glob1(dir_with_gz, "*.fastq.gz")) 
	
	numberFilesDowned= len(glob.glob1(dir_with_gz, "*.fastq.gz")) 
	
	bowtie_output_file=os.path.join(dir_with_gz, run_id + ".sam")
	
	bowtieLog = os.path.join(dir_with_gz, run_id + "_bowtie_output.txt")
	
		
	
	
	if numberFilesDowned==1:

		command_line ="bowtie2 -k 2 --quiet --no-unal -x "+bowtieBuildFileName+" -U "+dir_with_gz+"/"+run_id+".fastq.gz --rg-id ENA --rg PL:illumina --rg SM:"+run_id+" --sensitive-local --threads "+ str(threads) +" --met-file "+ os.path.join(dir_with_gz, run_id+".bowtie_metrics.txt") + " -S "+bowtie_output_file+" "


		myoutput = open(bowtieLog, 'w')
		args = shlex.split(command_line)
		p = subprocess.call(args,stdout=myoutput,stderr=myoutput)



	elif numberFilesDowned==2:
		command_line ="bowtie2 -k 2 --quiet --no-unal -x "+bowtieBuildFileName+" -1 "+dir_with_gz+"/"+run_id+"_1.fastq.gz -2 "+dir_with_gz+"/"+run_id+"_2.fastq.gz --rg-id ENA --rg PL:illumina --rg SM:"+run_id+" --sensitive-local --threads "+ str(threads) +" --met-file "+ os.path.join(dir_with_gz, run_id+".bowtie_metrics.txt") + " -S "+bowtie_output_file+" "
		
		myoutput = open(bowtieLog, 'w')
		args = shlex.split(command_line)
		p = subprocess.call(args,stdout=myoutput,stderr=myoutput)

		
		
	
	else:
		
		print "0 or more than 2 fastq files"


	return bowtie_output_file



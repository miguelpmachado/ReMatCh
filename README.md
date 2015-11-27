# ReMaChe - Read Mapping and Check-Coverage

# Usage 

`remache.py [-h] [-r REFERENCE_PATH] [-t TARGET_DIR] [-gatk GATK_PATH] [-cov MIN_COVERAGE] [-qual MIN_MAP_QUALITY] [-mul MULTIPLE_ALLELES] [-ri RUN_ID]`

# Description 

This program performs a download and read mapping of a given *RUN_ID* against a reference sequence and returns a report based on coverage.

Arguments:
 
  -h show this help message and exit

  -r REFERENCE_PATH 
  			Path for the reference sequence

  -t TARGET_DIR 
  			Output directory path
  
  -gatk GATK_PATH
        Path for the Genome Analysis Toolkit jar file
  
  -cov MIN_COVERAGE
  			Minimum coverage
  
  -qual
        Minimum mapping quality
  
  -mul
        Multiple alleles
  
  -ri
        Id of the sequencing run

# Example of usage

Mapping of the run ERR012354 against the reference.fasta file:

`python remache.py -r /reference.fasta -t results/ -gatk gatk_path/GenomAnalysisTK.jar -cov 10 -qual 0.6 -mul 0.1 -ri ERR012354`

#Dependencies

* Bowtie http://bowtie-bio.sourceforge.net/index.shtml

* NumPy http://www.numpy.org/

* SAMtools http://samtools.sourceforge.net/

* BEDtools http://bedtools.readthedocs.org/en/latest/

#Modules used

* NumPy http://www.numpy.org/

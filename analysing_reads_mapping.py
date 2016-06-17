import os
import sys

# Running commands
def runCommandPopenCommunicate(command_list):
	run_successfully = False
	print 'Running: ' + ' '.join(command_list)
	proc = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	if proc.returncode == 0:
		run_successfully = True
	else:
		print 'STDOUT'
		print stdout.decode("utf-8")
		print 'STDERR'
		print stderr.decode("utf-8")
	return run_successfully

# Indexing reference file
def indexSequence(sequenceFile):
	if os.path.isfile(str(sequenceFile + '.fai')):
		run_successfully = True
	else:
		command = ['samtools', 'faidx', sequenceFile]
		run_successfully = runCommandPopenCommunicate(command)
	return run_successfully

# Create indexed bam and remove sam
def samToBam(samFile, sequenceFile, threads):
	bamFile = os.path.splitext(samFile)[0] + '.bam'

	# Index sequence file
	run_successfully = indexSequence(sequenceFile)

	if not os.path.isfile(bamFile) and run_successfully:
		# sam to bam
		command = ['samtools', 'view', '-b', '-u', '-h', '-o', str(samFile + 'temp.bam'), '-T', sequenceFile, '-@', threads, samFile]
		run_successfully = runCommandPopenCommunicate(command)

		# Sort bam
		if run_successfully:
			# Remove samFile
			os.remove(samFile)
			# Sort bam
			command = ['samtools', 'sort', '-l', '0', '-o', bamFile, '-@', threads, str(samFile + 'temp.bam')]
			run_successfully = runCommandPopenCommunicate(command)

			# Index bam
			if run_successfully:
				# Remove intermediate bam
				os.remove(str(samFile + 'temp.bam'))
				# Index bam
				command = ['samtools', 'index', bamFile]
				run_successfully = runCommandPopenCommunicate(command)

	if not run_successfully:
		bamFile = None

	return run_successfully, bamFile

# Create gene vcf
def createVCF(bamFile, sequenceFile, geneName, minCoverage):
	geneVCF = os.path.splitext(bamFile)[0] + '.' + geneName.replace(' ', '_') + '.vcf'
	command = ['samtools', 'mpileup', '--count-orphans', '--no-BAQ', '--min-BQ', '13', '--fasta-ref', sequenceFile, '--region', str('"' + geneName + '"'), '--output', geneVCF, '--VCF', '--uncompressed', '--output-tags', 'AD,DP', '--min-ireads', str(minCoverage), bamFile]
	run_successfully = runCommandPopenCommunicate(command)
	if not run_successfully:
		geneVCF = None
	return run_successfully, geneVCF

# Read vcf file
class Vcf:
	self.line = ''
	def __init__(self, vcfFile):
		self.vcf = open(vcfFile, 'rtU')
		self.line_read = self.vcf.readline()
		while self.line_read.startswith('#'):
			self.line_read = self.vcf.readline()
		self.line = self.line_read
	def readline(self):
		self.line_stored = self.line
		try:
			self.line = self.vcf.readline()
		return self.line_stored
	def close(self):
		pass

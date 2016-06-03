import urllib2
import sys
import urllib
import xml.etree.ElementTree as ET
import argparse
import time



def main():

	#run -i SeqFromWebTaxon.py -i "acinetobacter baumannii" -o "allfastq.txt" -g True

	parser = argparse.ArgumentParser(description="This program gets a list of sequencing runs and machine were the sequencing was performed, given a taxon name accepted by the European nucleotide Archive")
	parser.add_argument('-i', nargs='?', type=str, help='taxon name', required=True)
	parser.add_argument('-o', nargs='?', type=str, help='output file name', required=True)
	parser.add_argument('-g', nargs='?', type=bool, help='True to include sequencing machine in the output', required=False)
	parser.add_argument('--getOmicsDataType', help='Informs the programme to include OMICS data type (GENOMIC / TRANSCRIPTOMIC / SYNTHETIC) in the output', action='store_true', required=False)



	args = parser.parse_args()
	
	try:
		getmachine = args.g

	except:
		getmachine = False
	taxonname = args.i
	outputfile = args.o
	getOmicsDataType = args.getOmicsDataType

def GetSequencesFromTaxon(taxonname,outputfile,getmachine, getOmicsDataType):

	taxonname=urllib.quote(taxonname)
	url="http://www.ebi.ac.uk/ena/data/view/Taxon%3A"+taxonname+"&display=xml"
	try:
		content = urllib2.urlopen(url)
		xml = content.read()
		tree = ET.fromstring(xml)
		taxonid=''
	except:
		print "Ooops!There might be a problem with the ena service, try later or check if the xml is well formated at " + url
		raise
	for child in tree:
		taxonid=child.get('taxId')
	if (taxonid):
		print "\n" + "Taxon ID found: "+taxonid
		url="http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree%28"+taxonid+"%29%22&result=read_run&display=xml"

		content = urllib2.urlopen(url)
		xml = content.read()
		tree = ET.fromstring(xml)
		runid=''
		n=0
		with open(outputfile, "wb") as f:
			f.write(str(time.strftime("%d/%m/%Y"))+"\n")
			model=''
			prjid=''
			models={}
			length_line = 0
			omics = ''
			for child in tree:
				runid=child.get('accession')
				
				n+=1
				

				if getmachine is True or getOmicsDataType:
					for child2 in child:
						if child2.tag == 'EXPERIMENT_REF':
							expid=child2.get('accession')
							url2="http://www.ebi.ac.uk/ena/data/view/"+expid+"&display=xml"
							content = urllib2.urlopen(url2)
							xml = content.read()
							tree2 = ET.fromstring(xml)
							try:
								for child3 in tree2:
									#print child3
									for child4 in child3:
										#print child4
										if child4.tag == 'PLATFORM':
											for child5 in child4:
												#print child5
												for child6 in child5:
													if child6.tag == 'INSTRUMENT_MODEL':
														model=child6.text
										elif child4.tag == 'STUDY_REF':
											prjid=child4.get('accession')
										elif child4.tag == 'DESIGN':
											if getOmicsDataType:
												for child5 in child4:
													if child5.tag == 'LIBRARY_DESCRIPTOR':
														for child6 in child5:
															if child6.tag == 'LIBRARY_SOURCE':
																omics = child6.text
							except:
								model='not found'
								omics = 'not found'
					f.write(str(runid)+"\t"+model+"\t"+prjid+"\t"+omics+"\n")
					line = "run acession %s sequenced on %s from project %s for %s data" % (runid, model, prjid, omics)
					if length_line < len(line):
						length_line = len(line)
					sys.stderr.write("\r" + line + ' '*(length_line - len(line)))
					sys.stderr.flush()
				else:
					f.write(str(runid)+"\n")
					line = "run acession %s" % (runid,prjid)
					if length_line < len(line):
						length_line = len(line)
					sys.stderr.write("\r" + line + ' '*(length_line - len(line)))
					sys.stderr.flush()
		print "\n"
		print "\nfound %s run id's" % n
		
	else:
		print "taxon name does not exist"	
	
if __name__ == "__main__":
	main()

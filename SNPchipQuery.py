import argparse
import os
import sys
import math
import pysam
from pybedtools import BedTool
from string import maketrans




# This function translates single nucleotides to their compliments
def Translate(SEQUENCE):
	intab = "ATGC"
	outtab = "TACG"
	trantab = maketrans(intab,outtab)
	return SEQUENCE.translate(trantab)

# This function takes in the A, C, G, T pileups at a single position, and computes a genotype from them.
# This genotype will be in the form of AA, AT, GG, GC, etc., where the first allele (if possible) is the reference allele  
def GetGenotype(ACOUNT,TCOUNT,GCOUNT,CCOUNT,REFALLELE,STRAND):
	if STRAND == '-':
		REFALLELE = Translate(REFALLELE)
# Make a dictionary of the base pileups
	BASEPILEUPS = {'A':ACOUNT,'C':CCOUNT,'T':TCOUNT,'G':GCOUNT}
# get the total coverage at this position
	COVERAGE = sum([ACOUNT,CCOUNT,TCOUNT,GCOUNT])
# Find the allele with the most coverage
	MAXALLELE1 = max(BASEPILEUPS,key=BASEPILEUPS.get)
#	print 'Max allele:',MAXALLELE1	
# Remove it from the dicitonary
	del BASEPILEUPS[MAXALLELE1]
# Find the allele with the second highest coverage
	MAXALLELE2 = max(BASEPILEUPS,key=BASEPILEUPS.get)
# if that allele has less than 10% of the total coverage, nix it, and make the second max allele equal to the first.
# These are cases of homozygous variants so we don't care:
	if BASEPILEUPS[MAXALLELE2] < (COVERAGE / 10):
		MAXALLELE2 = MAXALLELE1
		GENOTYPE = MAXALLELE1 + MAXALLELE2
		return GENOTYPE
	
		
#	print 'Second Max allele:',MAXALLELE2	
	if MAXALLELE1 == REFALLELE:
		GENOTYPE = MAXALLELE1 + MAXALLELE2
	elif MAXALLELE2 == REFALLELE:
		GENOTYPE = MAXALLELE2 + MAXALLELE1	
	else:
		GENOTYPE = MAXALLELE1+MAXALLELE2
	return GENOTYPE

	# Unit tests
	#print GetGenotype(50,0,1,0)
	#print GetGenotype(50,0,15,0)
	#print GetGenotype(0,10,20,0)
	#sys.exit()



# I didn't actually write this function, Aaron Odell did, but I've made modifications and a couple comments
class Variant:
	def __init__(self,chromosome,position,samfile):
		self.chromosome = chromosome
		self.position = position
		self.samfile = samfile
	def getInfo(self):
		### Here we are utilizing the pysam pileupcolumn object which is filled with pileupread objects which are filled with aligment objects
		self.nucDict = {'A':0,'C':0,'T':0,'G':0,'N':0}
		self.coverage = 0
		self.indelReads = 0
		#print '#################newVar############################'
		for pileupcolumn in self.samfile.pileup(self.chromosome,self.position+1):
			if pileupcolumn.pos == (self.position-1):
				for pileUpRead in pileupcolumn.pileups:
					if pileUpRead.is_del == 1:
						self.indelReads = self.indelReads + 1
					else:
						self.nucDict[pileUpRead.alignment.seq[pileUpRead.query_position]] = self.nucDict[pileUpRead.alignment.seq[pileUpRead.query_position]] + 1
						self.coverage = self.coverage + 1
				self.ann = self.chromosome+'\t'+str(self.position)+'\t'+'Coverage:'+str(self.coverage)+'\t'
				for i in self.nucDict:
					self.ann = self.ann+i+":"+str(self.nucDict[i])+'\t'
				self.ann = self.ann+"INDELReads:"+str(self.indelReads)
				break
			elif pileupcolumn.pos > (self.position + 10):
				break
		if self.coverage == 0:
			self.ann = self.chromosome+'\t'+str(self.position)+'\t'+'Coverage:'+str(self.coverage)+'\t'+'A:0'+'\t'+'C:0'+'\t'+'T:0'+'\t'+'G:0'+'\t'+'N:0'+'\t'+'INDELReads:'+str(self.indelReads)


def SnpArray(ARGS):
	# first I'll get the necessary parameters from ARGS
	query = open(ARGS.bed,'r')
        outFile = open(ARGS.o,'w')
        sampleID = ARGS.ID
	

	# The table file we're reading in has this format:
	#bin	chrom	chromStart	chromEnd	name	score	strand	refNCBI	refUCSC	observed	molType	class	valid	avHet	avHetSE	func	locType	weight	exceptions	submitterCount	submitters	alleleFreqCount	alleles	alleleNs	alleleFreqs	bitfields
	samfile = pysam.Samfile(args.bam,"rb")
	#ignore header
	query.readline()
	outFile.write('#Chromosome\tPosition\trsID\tREF\tStrand\tALTs\t%s\tA\tT\tG\tC\n'%sampleID)
	#print '#Chromosome\tPosition\trsID\tREF\tALTs\t%s\tA\tT\tG\tC\n'%sampleID
	for line in query:
		line = line.strip('\r')
		line = line.strip('\n')
		if line[0] != '#':
			columns = line.split('\t')
			chromosome = columns[1] #chromosome
			varPos = int(columns[3]) #position of variant
			rsid = columns[4]	#rs id
			strand = columns[6]	# had to incorporate strand, which will be used downstream to flip the genotype inferred from the bam file
			ref = columns[7]	# reference base
			combo = columns[9]	# some combo of ref/alt, but not consistent, so will infer alt using ref from columns 6, see below
						#Looks like this also has the ability to be A/C/T for triallelic variants, and A/C/G/T for quadallelic variants.
						# Going to have to deal with this in a more sophisticated manner below, for now will just retain combo
			# This is where I actually get the variant information
			varMan = Variant(chromosome,varPos,samfile)

	# next we will get the information for this variant slice
	# This is the actual call to pysam object, for the variant class that we designed above
			try:
				varMan.getInfo()
			except(ValueError):
				outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t..\t.\t.\t.\t.\n"%(chromosome,varPos,rsid,ref,strand,combo)) # if I can't get the variant, just print the ref and the combo, i.e. the A/C or the A/C/T
				continue
			# extract counts at this position for A, C, G, T
			if strand == '+':
				Acount = varMan.nucDict['A']
				Tcount = varMan.nucDict['T']
				Ccount = varMan.nucDict['C']
				Gcount = varMan.nucDict['G']
			# here I'm dealing with probes from the - strand
			elif strand == '-':
				Acount = varMan.nucDict['T']
				Tcount = varMan.nucDict['A']
				Ccount = varMan.nucDict['G']
				Gcount = varMan.nucDict['C']
			# Get the genotype and write to the output file
			genotype = GetGenotype(Acount,Tcount,Gcount,Ccount,ref,strand)
			outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n"%(chromosome,varPos,rsid,ref,strand,combo,genotype,Acount,Tcount,Gcount,Ccount))


def CytoscanHD(ARGS):
        query = open(ARGS.bed,'r')
        outFile = open(ARGS.o,'w')
        sampleID = ARGS.ID
	
        fasta = ARGS.Genome
	# The table we're reading in has this format:
	# chrom	start	stop	probeID	value	strand	start	stop
	samfile = pysam.Samfile(ARGS.bam,"rb")		
	# for now, I'll use this script to write out the A,C,G,T counts for each variant on the chip
	outFile.write("#Probe pileups for: %s\n"%sampleID)
	outFile.write("#Chrom\tstart\tstop\tprobeID\tA\tT\tG\tC\n")
	for line in query:
		line = line.strip('\r')
                line = line.strip('\n')
                if line[0] != '#':
			columns = line.split('\t')
			chromosome = columns[0]
			varPos = int(columns[2])
			start = columns[1]
			stop = columns[2]
			probeID = columns[3]
			refBase = BedTool.seq((chrom,varPos-1,varPos),fasta)
                        
			#populate the variant class
			varMan = Variant(chromosome,varPos,samfile)

			try:
				varMan.getInfo()
			except(ValueError):
				outFile.write("%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\n"%(chromosome,start,stop,probeID))
			Acount = varMan.nucDict['A']
                        Tcount = varMan.nucDict['T']
                        Ccount = varMan.nucDict['C']
                        Gcount = varMan.nucDict['G']
                        
			outFile.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n"%(chromosome,start,stop,probeID,Acount,Tcount,Ccount,Gcount))

if __name__  == "__main__":	
	parser = argparse.ArgumentParser()
	parser.add_argument("-bam",help="input your bam file",type=str)
	parser.add_argument("-bed",help="input the bed file of snps",type=str)
	parser.add_argument("-o",help="name of the ouput file",type=str)
	parser.add_argument("--ID",help="name of the sample ID",type=str)
	parser.add_argument("--Filetype",help="Two options: dbSNP || CytoscanHD",type=str)
        parser.add_argument("--Genome",help="Path to your reference genome sequence",type=str)

	args = parser.parse_args()

	# Initialize from the given arguments
	query = open(args.bed,'r')
	outFile = open(args.o,'w')
	sampleID = args.ID
        

	if args.Filetype == 'dbSNP':
		SnpArray(args)
	elif args.Filetype == 'CytoscanHD':
		CytoscanHD(args)
	
	
	



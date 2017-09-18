# ArrayFromWGS

This toolkit is designed to extract read-depth and genotype information from specific genomic sites in the genome.  The primary initial use (for me at least) is to utilize the Affy Cytoscan HD probe set in order to 



## Recipe (AKA Set-up and Installation)

1. Get dataset you want to use.  I'll demonstrate here with the Affymetrix Cytoscan:

wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/affyCytoScan.txt.gz  
gunzip -c affyCytoScan.txt.gz > affyCytoScan.txt  
grep -v "^#" affyCytoScan.txt | cut -f2,3,4,5,6,7,8 > AffyCytoScan_hg19.bed  

*Currently only implemented for h19*

2. Install required python libraries:

argparse, os, sys, math, pysam, string

3. Obtain your BAM file of interest

4. Run SNPchipQuery.py

python SNPchipQuery.py -h
optional arguments:
  -h, --help          show this help message and exit  
  -bam BAM            input your bam file  
  -bed BED            input the bed file of snps  
  -o O                name of the ouput file  
  -ID ID              name of the sample ID  
  -Filetype FILETYPE  Two options: dbSNP || CytoscanHD  

Example: 
python SNPchipQuery.py -bam WatsonGenome.bam -bed AffyCytoScan_hg19.bed -o WatsonGenome.affyCov.txt -ID Watson -Filetype CytoscanHD



  




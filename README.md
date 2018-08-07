# ArrayFromWGS

This toolkit is designed to extract read-depth and genotype information from specific genomic sites in the genome.  The primary initial use (for me at least) is to utilize the Affy Cytoscan HD probe set in order to 



## Recipe (AKA Set-up and Installation)

1. Get dataset you want to use.  I'll demonstrate here with the Affymetrix Cytoscan:
```
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/affyCytoScan.txt.gz  
gunzip -c affyCytoScan.txt.gz > affyCytoScan.txt  
grep -v "^#" affyCytoScan.txt | cut -f2,3,4,5,6,7,8 > AffyCytoScan_hg19.bed  
# For GRCh37
sed -e 's/^chr//g' AffyCytoScan_hg19.bed > AffyCytoScan_GRCh37.bed
```

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

BAM=/mnt/causes-data01/data/RICHMOND/Platinum/NA12878/N878_BWAmem_dupremoved_realigned.sorted.bam  
SAMPLE=NA12878  
OUPTUT=$SAMPLE_AffyCytoscan.tsv  

python SNPchipQuery.py \  
-bam $BAM \  
-bed AffyCytoScan_hg19.bed \  
-Filetype CytoscanHD \  
-ID $SAMPLE \  
-o $OUTPUT  



  




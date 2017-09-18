wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/affyCytoScan.txt.gz
gunzip -c affyCytoScan.txt.gz > affyCytoScan.txt
grep -v "^#" affyCytoScan.txt | cut -f2,3,4,5,6,7,8 > AffyCytoScan_hg19.bed

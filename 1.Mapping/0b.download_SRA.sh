#! /bin/bash

#run this on head node for internet connection

#get the other Zaprionus reads off of the SRA

#SRA numbers downloaded

cd /scratch/perickso/private/raw_data/published

while read SRR; do
  echo $SRR
   /scratch/perickso/private/sw/sratoolkit.3.0.1-ubuntu64/bin/fastq-dump --split-files "$SRR"
 done < SRR_Acc_List.txt 

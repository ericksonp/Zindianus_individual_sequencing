#add to config file

echo '#Z. indianus version Zi.PO1791'  >> /usr/local/sw/snpEff/snpEff.config

echo 'Zi.PO1791 : Z_indianus' >> /usr/local/sw/snpEff/snpEff.config

#put genome in the genomes folder
cp /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz /usr/local/sw/snpEff/data
mv /usr/local/sw/snpEff/data/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz /usr/local/sw/snpEff/data/genomes/Zi.PO1791.fa.gz

#put gff in its own folder
cp /scratch/perickso/private/raw_data/annotation/PO1791_Zaprionus_indianus.annotation.gff.gz  /usr/local/sw/snpEff/data
mv /usr/local/sw/snpEff/data/PO1791_Zaprionus_indianus.annotation.gff.gz /usr/local/sw/snpEff/data/Zi.PO1791/genes.gff.gz

#put proteins and cds files in folders
cp  /scratch/perickso/private/raw_data/annotation/PO1791_Zaprionus_indianus.protein.fasta.gz /usr/local/sw/snpEff/data/Zi.PO1791/protein.fa.gz
cp  /scratch/perickso/private/raw_data/annotation/PO1791_Zaprionus_indianus.transcript.fasta.gz /usr/local/sw/snpEff/data/Zi.PO1791/cds.fa.gz

#build database
java17 -jar /usr/local/sw/snpEff/snpEff.jar build -gff3 -v Zi.PO1791


#now use this database to annotate master vcf file


java17 -Xmx48g -jar /usr/local/sw/snpEff/snpEff.jar -v -c /usr/local/sw/snpEff/snpEff.config Zi.PO1791 /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gz > /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.annotated.vcf

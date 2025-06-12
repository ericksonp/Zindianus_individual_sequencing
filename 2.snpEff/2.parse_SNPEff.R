library(vcfR)
library(data.table)




#open VCF
vcf <- read.vcfR("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.annotated.vcf")

#see vcf metadata
queryMETA(vcf)

#see first 10 rows of fixed data
head(getFIX(vcf))

#extract fixed data into a data.table
y<-vcfR2tidy(vcf, info_only=T)
#y.ann<-extract.info(vcf, element="ANN")


##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos >  / AA.length | Distance | ERRORS / WARNINGS / INFO' 

data<-as.data.table(y$fix)
data<-data[,.(CHROM, POS, REF, ALT, ANN, LOF)]

data[,gene:=tstrsplit(ANN, split= "\\|")[[4]]]
data[,type:=tstrsplit(ANN, split= "\\|")[[2]]]

write.csv(data, file="/scratch/perickso/private/ind_seq/SnpEff_annotations_data_table.csv")

library(data.table)
freqs<-fread("/scratch/perickso/private/ind_seq/popgen/allelefreqsbySNP.csv")

load("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_CMHPO.Rdat")
ihs.VA<-as.data.table(wgscan.ihs$ihs)
setnames(ihs.VA, c("CHR", "POSITION"), c("chr", "position"))
bp.va.af<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesVsAfricaFemales/updated_data_for_manhattan.txt")
bp.va.af[,chr:=paste0("Scaffold_", Scaffold)]
bp.va.af[,position:=as.integer(tstrsplit(locations, split="_")[[3]])]
eff<-fread("/scratch/perickso/private/ind_seq/SnpEff_annotations_data_table.csv")
setnames(eff, c("CHROM", "POS"), c("chr", "position"))

eff<-merge(eff, ihs.VA[,.(chr, position, IHS)], by=c("chr", 'position'), all=T)
eff<-merge(eff, bp.va.af[,.(chr, position, M_XtX)], by=c("chr", "position"), all=T)
eff<-merge(eff, freqs[,.(snp.id, FL.freq, VA.freq, Africa.freq, chr, position)], by=c("chr", "position"))
eff[,freq.dif:=abs(VA.freq-Africa.freq)]


ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

ann[,annotation:=tstrsplit(description, split=";")[[1]]]
ann[,annotation:=tstrsplit(annotation, split="=")[[2]]]
ann[,gene_name:=tstrsplit(annotation, split="-")[[1]]]

#look at things near scaffold 2 peak

eff[chr=="Scaffold_2"&(M_XtX>7|IHS>5)&type=="missense_variant"]
eff[chr=="Scaffold_2"&(M_XtX>7&IHS>5)] #nothing
eff[chr=="Scaffold_2"&(M_XtX>7)&position>26500000&position<27500000&type=="missense_variant"]
eff[chr=="Scaffold_2"&(M_XtX>7)&position>26500000&position<27500000&gene=="ANN06929"]

ann[gene_name=="ANN06929"]
eff[gene=="ANN06929"&freq.dif>0.5]

ann[gene_name=="ANN06930"]
ann[gene_name=="ANN06922"]

eff[gene=="ANN06929"&type=="missense_variant"]

#Scaffold_2 26607469 is  p.Pro363Leu mutation (moderate effect)
#ref frequency is 1 in africa, 0.43 in Virginia, no signal of IHS but XtX is 7.63


eff[chr=="Scaffold_5"&(M_XtX>7|IHS>5)&type=="missense_variant"]
eff[chr=="Scaffold_5"&(M_XtX>7&IHS>5)] #nothing
eff[chr=="Scaffold_5"&(M_XtX>7)&position>6500000&position<7500000&type=="missense_variant"]
eff[chr=="Scaffold_5"&(freq.dif>.5)&position>6500000&position<7500000&type=="missense_variant"]


ann[gene_name=="ANN03540"] #Cyp6a9
ann[gene_name=="ANN03541"] #Cyp317a
#Scaffold_5 7892185 is a missense variant in ANN03540. VA freq is 0.38, Africa freq is 1. Cys to Ser (moderate)
#Scaffold_5 7900184 is a missense variant in ANN03541. VA freq is 0.39, Africa freq is 1. Arg to Cys (moderate)

#the protein homologies hold up based on blast searches

#This work also implicates cytochrome P450s in resistance to permethrin. DGRP variants most strongly associated with permethrin map to a region on chromosome 2R containing nine P450 genes, with peaks over Cyp6a23 and Cyp317a1. 

#anythign at Ace?
ann[info=="gene"][grep("Ace",description)]
#ace is ANN01838
#ace is found at Scaffold 1, 21507820
ann[gene_name=="ANN02078"]


ann[grep(" Esterase", description)]


#what about the IHS peak on X

ann[chr=="Scaffold_3"&info=="gene"&((start>500000&start<1250000)|(end>500000&end<1250000))]

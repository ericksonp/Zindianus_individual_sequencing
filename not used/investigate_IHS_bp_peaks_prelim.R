bp.peaks<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/zap_ind_2023_outlier_SNPs_q999.xtx")
setnames(bp.peaks, "SNP.info[outlier_SNPs_999$MRK, ]", "snp")
bp.peaks[,CHR:=paste0(tstrsplit(snp, split="_")[[1]],"_", tstrsplit(snp, split="_")[[2]]) ]
bp.peaks[,POSITION:=as.integer(tstrsplit(snp, split="_")[[3]])]
bp.peaks[order(M_XtX, decreasing=T)][CHR!="Scaffold_3"]
#Scaffold_1 20306622
#Scaffold_2 20748414
#Scaffold_4  6573623
#Scaffold_1 16187187
#Scaffold_4 16209741
#Scaffold_2_26766536
ann<-fread("/scratch/perickso/private/annotation/PO1791_Zaprionus_indianus.annotation.gff")
setnames(ann, c("chr", "type", "info", "start", "end", "score", "strand", "other", "description"))

#what's within 10kb of each bayepass peak?

ann[chr=="Scaffold_1"&end>20306622-10000&start<20306622+10000&info=="gene"]
#adipocyte plasma protein
#unk (unkempt) -eye development

#this looks like it is close ot the scaffold 1 inversion breakpoint

ann[chr=="Scaffold_1"&end>16187187-10000&start<16187187+10000&info=="gene"]
#esterase-B1--insectiside resistance in Culex!!
# Scaffold_1	16177203 is a mis-sense variant in ANN01440 but is alanine to valine so not a big change
#maybe use vcfR to look for all NS variants in these genes

ann[chr=="Scaffold_2"&end>20748414-10000&start<20748414+10000&info=="gene"]

#amnionless - protein resporbtion and vitamin b12 in nephrocytes
#smoothened - hedgehog signaling

#this SNP is very close to the scaffold 2 inversion breakpoint

ann[chr=="Scaffold_2"&end>26766536-10000&start<26766536+10000&info=="gene"]
#no annotations

ann[chr=="Scaffold_4"&end>6573623-10000&start<6573623+10000&info=="gene"]
#rho guanine nucleotide exchange factor


ann[chr=="Scaffold_4"&end>16209741-10000&start<16209741+10000&info=="gene"]
#tbox transcription factor

#do the bayepass peaks come close to the inversion breakpoints?

inv<-fread("/scratch/perickso/private/ind_seq/popgen/LDdecay_inversions.bed", header=F)


#what about the scaffold 5 location that is a peak for all 3 selection tests?
wins<-fread("/scratch/perickso/private/ind_seq/popgen/ihs_bp_fst_1kbwindows.csv")
wins[CHR==5&all.3==T]
ann[chr=="Scaffold_5"&end>7844001-10000&start<7885001+10000&info=="gene"]

#this is a cluster of CYP450 genes

#read in snpeff data
eff<-fread("/scratch/perickso/private/ind_seq/SnpEff_annotations_data_table.csv")
eff[(gene=="ANN03537"|gene=="ANN03538"|gene=="ANN03539"|gene=="ANN03540")&type=="missense_variant"]
eff[(gene=="ANN03537")&type=="missense_variant"]

#IHS peak is scaffold 5 7875359
eff[POS==7875359] #synonymous variant but several missense variants nearby
load("/scratch/perickso/private/ind_seq/popgen/FST_VAvsFL.Rdata")
z[position>7844001&position<7886000&chromosome=="Scaffold_5"][order(fst.snp, decreasing=T)]
eff[POS==7844664]



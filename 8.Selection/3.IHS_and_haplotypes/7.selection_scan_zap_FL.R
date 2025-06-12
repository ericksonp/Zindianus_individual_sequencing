library(data.table)
library(rehh)
library(foreach)
library(vcfR)

#scan autosomes together
scaffolds<-fread("/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai", header=F)
setnames(scaffolds, "V1" ,"chr")
scaffolds[,index:=1:nrow(scaffolds)]

metadata<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)



wgscan<-foreach(i=c(1:5), .errorhandling="remove")%do% {
    print(i)
        hap_file=paste0("/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.", i, ".vcf") #change to ancestral after it's made
        hh <- data2haplohh(hap_file = hap_file,
                           chr.name = paste0("Scaffold_", i),
                           min_perc_geno.mrk = 90,
                           polarize_vcf = TRUE,
                           vcf_reader = "vcfR")
        samps<-data.table(id=hap.names(hh))
        samps[,haplotype.id:=1:(nrow(samps))]
        samps[,sample.id:=substr(id,1,nchar(id)-2 )]
        samps[,hapnum:=rep(c(1:2), times=.N/2)]

        samps<-merge(metadata, samps, by="sample.id")
        if(i==3){
        haps.to.use<-samps[(assigned_sex=="F"|(assigned_sex=="M"&hapnum==1))&(Location=="FL"|Location=="MIA"), haplotype.id] #can add other conditionals here for specific subsets 
        }
        else{
          haps.to.use<-samps[Location=="MIA"|Location=="FL", haplotype.id] #can add other conditionals here for specific subsets 
        }
    
        hh_subset = subset(hh, select.hap = haps.to.use, min_perc_geno.mrk = 90, min_maf=0)
        scan <- scan_hh(hh_subset, discard_integration_at_border = FALSE, polarized=TRUE)
        
        
        return(scan)
}

wgscan<-rbindlist(wgscan )
print("starting whole genome scan")
wgscan.ihs <- ihh2ihs(wgscan, freqbin = 0.05)

save(wgscan.ihs, file=paste0("/scratch/perickso/private/ind_seq/popgen/rehh_wgscan_bychr_FL.Rdat"))


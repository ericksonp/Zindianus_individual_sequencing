srun \
--pty \
-t 6:00:00 \
--mem=100G \
--partition erickson \
--ntasks-per-node=24 \
bash

# Code from https://github.com/JimWhiting91/guppy_convergence/blob/main/BayPass/scripts

require(ape) ; require(data.table)

setwd("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales")

#######################################################
# XtX calibration
#######################################################
#get the pod XtX
pod.xtx=read.table("zap_ind_2023_sim_summary_pi_xtx.out",h=T)$M_XtX
#anacore.snp.res = read.table("zaprionus_individual_2023_core_3a_summary_pi_xtx.out",h=T) #REAL DATAAAAAA ADDDDD
anacore.snp.res = read.table("zap_ind_2023_core_summary_pi_xtx.out",h=T)

#compute the 1% threshold
pod.thresh99=quantile(pod.xtx,probs=0.99) #USE THIS AS THE CUTOFF NUMBER
pod.thresh999=quantile(pod.xtx,probs=0.999)

thresh = data.table(pod.thresh99, pod.thresh999)
fwrite(thresh, "zap_thresholds.txt")
# #add the thresh to the actual XtX plot
# plot(zaprionus_individual_2023_2_core_summary_pi_xtx.out$M_XtX)
# abline(h=pod.thresh99,lty=2)

####################################################
# Highlight 'outlier' SNPs
####################################################
outlier_SNPs_99<-anacore.snp.res[anacore.snp.res$M_XtX > pod.thresh99,] #so that is in the file and is okayyyyyy
outlier_SNPs_999<-anacore.snp.res[anacore.snp.res$M_XtX > pod.thresh999,]

# Add in SNP info
outlier_SNPs_99<-cbind(outlier_SNPs_99,SNP.info[outlier_SNPs_99$MRK,]) #not quite sure what this SNP info is supposed to have
outlier_SNPs_999<-cbind(outlier_SNPs_999,SNP.info[outlier_SNPs_999$MRK,])

# Write to output
write.table(outlier_SNPs_99,"zap_ind_2023_outlier_SNPs_q99.xtx",row.names=F,sep="\t",quote=F)
write.table(outlier_SNPs_999,"zap_ind_2023_outlier_SNPs_q999.xtx",row.names=F,sep="\t",quote=F)

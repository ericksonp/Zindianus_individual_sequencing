# Introduction
This folder contains scripts for mapping raw reads, calling SNPs and assembling the vcf used for downstream analyses for individual sequencing of Z. indianus specimens for Erickson et al.

#Description of scripts

<0b.download_SRA.sh> downloads data published from Comeault et al 2020 and 2021

<1b.map_zap_published2.sh> maps raw reads to the genome, merges overlapping reads and produces final bam org_files

<1c.coverage_zap_published.sh> runs an analysis of sequencing coverage_zap_published

<1d.coverage_analysis.R> is a script to analyze and plot the results of 1c

<2.hap_zap_published.sh> runs GATK's HaplotypeCaller on mapped reads to produce gVCFs for each sample

<2b.vcfnames.R> produces the input file of sample names for haplotype calling

<3c.genomicsDBI_zap_updated.sh
> imports gVCFs into a genomicsDBI object

<4b.genotype_zap_addpublished.sh>

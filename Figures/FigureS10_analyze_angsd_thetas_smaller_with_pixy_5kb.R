library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggpubfigs)
library(ggsci)


pixy.pi<-foreach(i=c(1,2,4,5))%do%{

    data<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/pixy/pixy_subpops_5kb_Scaffold_", i, "_pi.txt"))

  return(data)
}
pixy.pi<-rbindlist(pixy.pi)
pixy.pi[,season:=tstrsplit(pop, split="_")[[3]]]
#pixy.pi[season=="mid", season:="early"]
pixy.pi[,year:=tstrsplit(pop, split="_")[[2]]]
pixy.pi[,loc:=tstrsplit(pop, split="_")[[1]]]
pixy.pi[loc=="MIA", season:="Florida"]
pixy.pi[loc=="MIA", loc:="FL"]
pixy.pi<-pixy.pi[no_sites>2000]
pixy.pi[,season:=factor(season, levels=c("Florida", "early", "mid", "late"))]
pixy.pi[,loc.year:=paste(loc, year, sep="-")]



pixy.sum<-pixy.pi[,.(mean.pi=mean(avg_pi, na.rm=T),
                     sd.pi=sd(avg_pi, na.rm=T),
                     median.pi=median(avg_pi,na.rm=T), 
                     mad.pi=mad(avg_pi, na.rm=T),
                     q05=quantile(avg_pi, .05, na.rm=T),
                     q95=quantile(avg_pi, .95, na.rm=T),
                     n=.N), .(year, season, loc, pop)]

pixy.sum[,season:=factor(season, levels=c("Florida", "early", "mid", "late"))]
pixy.sum[,loc.year:=paste(loc, year, sep="-")]

pixy.sum[,se.pi:=sd.pi/sqrt(n)]


a<-ggplot(pixy.sum)+
  geom_point(aes(x=as.factor(season), y=mean.pi, color=season))+
  geom_errorbar(aes(x=as.factor(season), ymin=mean.pi-se.pi, ymax=mean.pi+se.pi, color=season), width=0.5)+
  facet_grid(.~loc.year, scales="free_x", space='free')+
  labs(x="", y="π, 5 kb windows")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  guides(color="none")+                                                                # Change font size
  theme(strip.text.x = element_text(size = 10))+
  lims(y=c(0.0085, 0.0095))

a<-ggplot(pixy.sum)+
  geom_point(aes(x=as.factor(season), y=median.pi, color=season))+
  geom_errorbar(aes(x=as.factor(season), ymin=median.pi-mad.pi, ymax=median.pi+mad.pi, color=season), width=0.5)+
  facet_grid(.~loc.year, scales="free_x", space='free')+
  labs(x="", y="π, 5 kb windows")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  guides(color="none")+                                                                # Change font size
  theme(strip.text.x = element_text(size = 10))+
  lims(y=c(0.002,0.016))


#try some K-S tests comparing Florida to early and early to late in each year
ks.test(pixy.pi[pop=="MIA_2019_June",avg_pi], pixy.pi[pop=="VA-CM_2017_early", avg_pi])# 1 x 10-6
ks.test(pixy.pi[pop=="MIA_2019_June",avg_pi], pixy.pi[pop=="VA-CM_2018_early", avg_pi])#.0158
ks.test(pixy.pi[pop=="MIA_2019_June",avg_pi], pixy.pi[pop=="VA-CM_2019_mid", avg_pi])#4x10-7
ks.test(pixy.pi[pop=="MIA_2019_June",avg_pi], pixy.pi[pop=="VA-CM_2020_early", avg_pi])#2x10-16
ks.test(pixy.pi[pop=="MIA_2019_June",avg_pi], pixy.pi[pop=="VA-HPO_2019_mid", avg_pi])#.016
ks.test(pixy.pi[pop=="MIA_2019_June",avg_pi], pixy.pi[pop=="VA-HPO_2020_early", avg_pi])#2x10-16

ks.test(pixy.pi[pop=="VA-CM_2017_early",avg_pi], pixy.pi[pop=="VA-CM_2017_late", avg_pi])#.04
ks.test(pixy.pi[pop=="VA-CM_2018_early",avg_pi], pixy.pi[pop=="VA-CM_2018_late", avg_pi])#2x10-16
ks.test(pixy.pi[pop=="VA-CM_2019_mid",avg_pi], pixy.pi[pop=="VA-CM_2019_late", avg_pi])#7x10-9



groups<-fread("/scratch/perickso/private/ind_seq/popgen/org_files/populations.txt", header=F)$V1

thetas<-foreach(pop=groups, .errorhandling="remove")%do%{
  theta<-fread(paste0("/scratch/perickso/private/ind_seq/popgen/angsd/", pop, ".thetasWindow5000_5000.pestPG"))
  theta[,group:=pop]
  return(theta)
}

thetas<-rbindlist(thetas)
thetas<-thetas[Chr%in%c("Scaffold_1", "Scaffold_2", "Scaffold_4", "Scaffold_5")]

thetas[,thetaW:=tW/nSites]
thetas[,thetaPi:=tP/nSites]

thetas[,year:=tstrsplit(group, split="_")[[2]]]
thetas[,loc:=tstrsplit(group, split="_")[[1]]]
thetas[,season:=tstrsplit(group, split="_")[[3]]]

#remove windows wiht a lot of missing data (less than 9K sites in 10K window)
thetas<-thetas[nSites>=4500]

#thetas[season=="mid", season:="early"]
thetas[loc=="MIA", season:="Florida"]
thetas[loc=="MIA", loc:="FL"]
theta.sum<-thetas[,.(mean.theta=mean(thetaW), 
                     sd.theta=sd(thetaW), 
                     mean.pi=mean(thetaPi),
                     sd.pi=sd(thetaPi),
                     median.theta=median(thetaW, na.rm=T),
                     mad.theta=mad(thetaW, na.rm=T),
                     n=.N, 
                     mean.tj=mean(Tajima), 
                     sd.tj=sd(Tajima), 
                     median.tj=median(Tajima, na.rm=T),
                     mad.tj=mad(Tajima,na.rm=T)), 
                  .(year, season, loc, group)]
theta.sum[,se.theta:=sd.theta/sqrt(n)]
theta.sum[,se.pi:=sd.pi/sqrt(n)]
theta.sum[,se.tj:=sd.tj/sqrt(n)]

theta.sum[,season:=factor(season, levels=c("Florida", "early", "mid", "late"))]
theta.sum[,loc.year:=paste(loc, year, sep="-")]


b<-ggplot(theta.sum[!is.na(loc)])+
  geom_point(aes(x=as.factor(season), y=mean.theta, color=season))+
  geom_errorbar(aes(x=as.factor(season), ymin=mean.theta-se.theta, ymax=mean.theta+se.theta, color=season), width=0.5)+
  facet_grid(.~loc.year, scales="free_x", space='free')+
  labs(x="", y="θ, 5 kb windows")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  guides(color="none")+                                                                # Change font size
  theme(strip.text.x = element_text(size = 10))+
  lims(y=c(0.015, 0.01675))


b<-ggplot(theta.sum[!is.na(loc)])+
  geom_point(aes(x=as.factor(season), y=median.theta, color=season))+
  geom_errorbar(aes(x=as.factor(season), ymin=median.theta-mad.theta, ymax=median.theta+mad.theta, color=season), width=0.5)+
  facet_grid(.~loc.year, scales="free_x", space='free')+
  labs(x="", y="θ, 5 kb windows")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  guides(color="none")+                                                                # Change font size
  theme(strip.text.x = element_text(size = 10))+
  lims(y=c(0.006, 0.028))



ks.test(thetas[group=="MIA_2019_June",thetaW], thetas[group=="VA-CM_2017_early", thetaW])# 2 x 10-16
ks.test(thetas[group=="MIA_2019_June",thetaW], thetas[group=="VA-CM_2018_early", thetaW])#.002
ks.test(thetas[group=="MIA_2019_June",thetaW], thetas[group=="VA-CM_2019_mid", thetaW])#2.7x10-6
ks.test(thetas[group=="MIA_2019_June",thetaW], thetas[group=="VA-CM_2020_early", thetaW])#2x10-16
ks.test(thetas[group=="MIA_2019_June",thetaW], thetas[group=="VA-HPO_2019_mid", thetaW])#6x10-6
ks.test(thetas[group=="MIA_2019_June",thetaW], thetas[group=="VA-HPO_2020_early", thetaW])#.007

ks.test(thetas[group=="VA-CM_2017_early",thetaW], thetas[group=="VA-CM_2017_late", thetaW])#.104
ks.test(thetas[group=="VA-CM_2018_early",thetaW], thetas[group=="VA-CM_2018_late", thetaW])#2x10-16
ks.test(thetas[group=="VA-CM_2019_mid",thetaW], thetas[group=="VA-CM_2019_late", thetaW])#2x10-9=15

c<-ggplot(theta.sum[!is.na(loc)])+
  geom_point(aes(x=as.factor(season), y=mean.tj, color=season))+
  geom_errorbar(aes(x=as.factor(season), ymin=mean.tj-se.tj, ymax=mean.tj+se.tj, color=season), width=0.5)+
  facet_grid(.~loc.year, scales="free_x", space='free')+
  labs(x="", y="Tajima's D,\n 5 kb windows")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  guides(color="none")+                                                                # Change font size
  theme(strip.text.x = element_text(size = 10))+
  lims(y=c(0.35, 0.85))


c<-ggplot(theta.sum[!is.na(loc)])+
  geom_point(aes(x=as.factor(season), y=median.tj, color=season))+
  geom_errorbar(aes(x=as.factor(season), ymin=median.tj-mad.tj, ymax=median.tj+mad.tj, color=season), width=0.5)+
  facet_grid(.~loc.year, scales="free_x", space='free')+
  labs(x="", y="Tajima's D,\n 5 kb windows")+
  scale_color_manual(values = friendly_pal("ito_seven"))+
  guides(color="none")+                                                                # Change font size
  theme(strip.text.x = element_text(size = 10))+
  lims(y=c(-0.2, 1.5))




ks.test(thetas[group=="MIA_2019_June",Tajima], thetas[group=="VA-CM_2017_early", Tajima])# 2 x 10-16
ks.test(thetas[group=="MIA_2019_June",Tajima], thetas[group=="VA-CM_2018_early", Tajima])#.2x10-12
ks.test(thetas[group=="MIA_2019_June",Tajima], thetas[group=="VA-CM_2019_mid", Tajima])#2x10-16
ks.test(thetas[group=="MIA_2019_June",Tajima], thetas[group=="VA-CM_2020_early", Tajima])#2x10-16
ks.test(thetas[group=="MIA_2019_June",Tajima], thetas[group=="VA-HPO_2019_mid", Tajima])#2x10-16
ks.test(thetas[group=="MIA_2019_June",Tajima], thetas[group=="VA-HPO_2020_early", Tajima])#2x10-16

ks.test(thetas[group=="VA-CM_2017_early",Tajima], thetas[group=="VA-CM_2017_late", Tajima])#2x10-16
ks.test(thetas[group=="VA-CM_2018_early",Tajima], thetas[group=="VA-CM_2018_late", Tajima])#2x10-16
ks.test(thetas[group=="VA-CM_2019_mid",Tajima], thetas[group=="VA-CM_2019_late", Tajima])#2x10-16
jpeg("/scratch/perickso/private/ind_seq/popgen/plots/diversity_stats_5kbwindows_pixy_and_angsd_median_mad.jpg", height=9, width=10,units="in", res=300 )

plot_grid(a,b, c, nrow=3, align="v", axis="l", labels=c("A", "B", "C"))
dev.off()


jpeg("/scratch/perickso/private/ind_seq/popgen/plots/diversity_stats_5kbwindows_pixy_and_angsd_median_mad.jpg", height=9, width=10,units="in", res=300 )

plot_grid(a,b, c, nrow=3, align="v", axis="l", labels=c("A", "B", "C"))
dev.off()


jpeg("/scratch/perickso/private/ind_seq/popgen/plots/diversity_stats_5kbwindows_pixy_and_angsd_extraspace.jpg", height=9, width=10,units="in", res=300 )

plot_grid(a,b, c, nrow=3, align="v", axis="l", labels=c("A", "B", "C"))
dev.off()



library(data.table)
library(foreach)
library(ggplot2)
library(ggpubfigs)
library(cowplot)
theme_set(theme_cowplot())

psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
  X<-scan(file=file,what="",sep="\n",quiet=TRUE)
  
  START<-grep("^RD",X)
  END<-grep("^//",X)
  
  X<-X[START[i.iteration+1]:END[i.iteration+1]]
  
  TR<-grep("^TR",X,value=TRUE)
  RS<-grep("^RS",X,value=TRUE)
  
  write(TR,"temp.psmc.result")
  theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
  N0<-theta0/4/mu/s
  
  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-as.numeric(2*N0*a[,3])
  Ne<-as.numeric(N0*a[,4])
  
  file.remove("temp.psmc.result")
  
  n.points<-length(Ne)
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
  data.frame(YearsAgo,Ne)
}

psmc.final<-foreach(p=pops, .combine="rbind", .errorhandling="remove")%do%{
  files<-list.files(path=paste0("/scratch/perickso/private/ind_seq/popgen/phlash/", p, "/"), pattern="1111.psmc$")
#making a loop to read files
  psmc<-foreach(i=files, .combine="rbind")%do%{
    x<-as.data.table(psmc.result(file=paste0("/scratch/perickso/private/ind_seq/popgen/phlash/", p, "/", i), i.iteration = 25, mu=2.8e-9, s=100, g=0.08))
    x[,sample:=i]
    return(x)
  }
#combine tables
return(psmc)
  }

#adding names
info<-fread("/scratch/perickso/private/ind_seq/zap_full_info_updated_v2.csv", drop=1)

#create new column with the sample name coded above
psmc.final[,sample.id:=tstrsplit(sample,"[.]")[[1]]]
#merge tables based on shared column name
psmc.merge<-merge(psmc.final,info, by="sample.id")
ggplot(psmc.merge)+geom_line(aes(x=log10(YearsAgo), y=log10(Ne),group=sample,color=Location))+facet_wrap(~Location)

psmc.merge[loc.spec=="Africa", type2:="native"]
psmc.merge[loc.spec=="VA-CM"|loc.spec=="VA-HPO", type2:="introduced-Virginia"]
psmc.merge[loc.spec=="FL", type2:="introduced-Florida"]
psmc.merge[,type2:=factor(type2, levels=c("introduced-Florida", "introduced-Virginia",  "native" ))]
psmc.merge<-psmc.merge[order(type2)]

ggplot(psmc.merge)+geom_line(aes(x=log10(YearsAgo), y=log10(Ne),group=sample,color=type2))+ 
  scale_color_manual(values = friendly_pal("contrast_three"))





#specifying locations to be plotted
ggplot(psmc.merge[Location=="MIA"|Location=="HPO"|Location=="CM"])+geom_line(aes(x=log10(YearsAgo), y=log10(Ne),group=sample,color=Location))

#africa plot
ggplot(psmc.merge[Location%in%c("Kenya","SaoTome","SenegalDesert","SenegalForest","Zambia")])+geom_line(aes(x=log10(YearsAgo), y=log10(Ne),group=sample,color=Location))
#TN, HI, NC
ggplot(psmc.merge[Location%in%c("NC", "TN", "HI")])+geom_line(aes(x=log10(YearsAgo), y=log10(Ne),group=sample,color=Location))
#iterations of africa
psmc.iterations<-foreach(i=c(1:25))%do%{
  x<-as.data.table(psmc.result(file="SRR11077360.psmc", i.iteration = i, mu=8.4e-9, s=100, g=0.08))
  x[,iteration:=i]
  return(x)
}
psmc.final.it<-rbindlist(psmc.iterations)

ggplot(psmc.final.it)+geom_line(aes(x=log10(YearsAgo), y=log10(Ne),group=iteration,color=iteration))


#plot based on african or north american

psmc.merge[,continent:=ifelse(Location%in%c("Kenya","SaoTome","SenegalDesert","SenegalForest","Zambia"), "Africa", "NorthAmerica")]
ggplot(psmc.merge)+geom_line(aes(x=log10(YearsAgo), y=log10(Ne),group=sample,color=continent))

                             
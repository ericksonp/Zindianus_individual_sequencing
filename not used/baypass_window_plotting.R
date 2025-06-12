
library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)



bp.va.fl<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/updated_data_for_manhattan.txt")
bp.va.af<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesVsAfricaFemales/updated_data_for_manhattan.txt")
bp.fl.af<-fread("/scratch/perickso_shared/alexandra/BayPassAlexEdited/allFLFemalesVsAfricaFemales/updated_data_for_manhattan.txt")

bp.va.fl[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]
bp.va.af[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]
bp.fl.af[,pos:=as.integer(tstrsplit(locations, split="_")[[3]])]

win.bp <- 5e3
step.bp <- 5e3

#FL VS VA
wins<-foreach(chr.i=c(1:5),
              .combine="rbind", 
              .errorhandling="remove")%dopar%{
                
                tmp <- bp.va.fl[Scaffold==chr.i]
                data.table(CHR=chr.i,
                           start=seq(from=1, to=max(tmp$pos)-win.bp, by=step.bp),
                           end=seq(from=1, to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
              }

wins[,index:=c(1:nrow(wins))]
setkey(wins,index)

bp.sum<-foreach(window=wins$index, .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  if(window%%1000==0){
                    print(window)
                  }
                  bp.tmp<-bp.va.fl[Scaffold==wins[J(window),CHR] & pos>=wins[J(window), start] & pos<wins[J(window),end]]
                  stat<-mean(bp.tmp$M_XtX, na.rm=T)
                  data.table(index=window,
                             window.xtxst=stat,
                             n.snp=nrow(bp.tmp))
                }



bp.sum<-merge(bp.sum, wins, by="index")

write.csv(bp.sum, file="/scratch/perickso/private/ind_seq/popgen/baypass/allVAFemalesvsFLfemales_5kb_window_mean.csv")

#VA vs Africa

wins<-foreach(chr.i=c(1:5),
              .combine="rbind", 
              .errorhandling="remove")%dopar%{
                
                tmp <- bp.va.af[Scaffold==chr.i]
                data.table(CHR=chr.i,
                           start=seq(from=1, to=max(tmp$pos)-win.bp, by=step.bp),
                           end=seq(from=1, to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
              }

wins[,index:=c(1:nrow(wins))]
setkey(wins,index)

bp.sum<-foreach(window=wins$index, .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  if(window%%1000==0){
                    print(window)
                  }
                  bp.tmp<-bp.va.af[Scaffold==wins[J(window),CHR] & pos>=wins[J(window), start] & pos<wins[J(window),end]]
                  stat<-mean(bp.tmp$M_XtX, na.rm=T)
                  data.table(index=window,
                             window.xtxst=stat,
                             n.snp=nrow(bp.tmp))
                }


bp.sum<-merge(bp.sum, wins, by="index")

write.csv(bp.sum, file="/scratch/perickso/private/ind_seq/popgen/baypass/allVAFemalesvsAfricafemales_5kb_window_mean.csv")


#FL vs Africa

wins<-foreach(chr.i=c(1:5),
              .combine="rbind", 
              .errorhandling="remove")%dopar%{
                
                tmp <- bp.fl.af[Scaffold==chr.i]
                data.table(CHR=chr.i,
                           start=seq(from=1, to=max(tmp$pos)-win.bp, by=step.bp),
                           end=seq(from=1, to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
              }

wins[,index:=c(1:nrow(wins))]
setkey(wins,index)

bp.sum<-foreach(window=wins$index, .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  if(window%%1000==0){
                    print(window)
                  }
                  bp.tmp<-bp.fl.af[Scaffold==wins[J(window),CHR] & pos>=wins[J(window), start] & pos<wins[J(window),end]]
                  stat<-mean(bp.tmp$M_XtX, na.rm=T)
                  data.table(index=window,
                             window.xtxst=stat,
                             n.snp=nrow(bp.tmp))
                }


bp.sum<-merge(bp.sum, wins, by="index")

write.csv(bp.sum, file="/scratch/perickso/private/ind_seq/popgen/baypass/allFLFemalesvsAfricafemales_5kb_window_mean.csv")

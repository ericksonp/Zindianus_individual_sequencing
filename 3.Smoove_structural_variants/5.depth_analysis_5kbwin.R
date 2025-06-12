library(data.table)
library(runner)
library(foreach)
library(doMC)
registerDoMC(50)

samps<-fread("/scratch/perickso/private/ind_seq/sv/all_samples.txt", header=F)$V1


#cycle through all samples

all.depth<-foreach(x=samps, .errorhandling = "remove") %dopar% {
  print(x)
  y<-foreach(chr=c(1:5), .errorhandling="remove") %do% {
    chr.depth<-fread(paste0("/scratch/perickso/private/ind_seq/mapping_stats/", x, "_scaffold_", chr, ".depth.txt"), header=F)
    mean.depth=mean(chr.depth$V3)
    rel.win.depth=runner(chr.depth$V3, f = function(x) mean(x), at=seq(from=1, to=nrow(chr.depth), by=5000), k=5000)/mean.depth
    data<-data.table(sample.id=x,
                    chr=chr,
                     win=seq(from=1, to=nrow(chr.depth), by=5000),
                     rel_depth=rel.win.depth)
    return(data)
  }
  y<-rbindlist(y)
  return(y)
  }

all.depth<-rbindlist(all.depth)

fwrite(all.depth, file= "/scratch/perickso/private/ind_seq/sv/depth_analysis_5kb.csv")
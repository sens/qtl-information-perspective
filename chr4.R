#######################################################
#                Figure 2
#######################################################

library(qtl)
data(hyper)
idx <- order(hyper$pheno[,1])
plot.geno(hyper,chr=4,ind=idx[seq(0,250,by=5)],hor=T)

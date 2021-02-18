#####################################################
#             Figure 6
#####################################################

library(qtlDesign)

d <- 100/seq(2,100,by=0.2)

opt.alpha1 <- d
opt.alpha01 <- d
opt.alpha001 <- d
opt.alpha0001 <- d

for( i in 1:length(d) )
  {
    opt.alpha1[i] <- optalpha.bc(cost=1,d=d[i]/100,G=14.5)
    opt.alpha01[i] <- optalpha.bc(cost=0.1,d=d[i]/100,G=14.5)
    opt.alpha001[i] <- optalpha.bc(cost=0.01,d=d[i]/100,G=14.5)
    opt.alpha0001[i] <- optalpha.bc(cost=0.001,d=d[i]/100,G=14.5)
  }


###########
plot(d,opt.alpha1,ylim=c(0,1),pch=".",ylab="",xlab="",axes=F)
axis(1, at = seq(0,50,by=10), labels = formatC(seq(0,50,by=10), format="fg"))
axis(2, at = seq(0,100,by=10)/100
     , labels = formatC(seq(0,100,by=10), format="fg"))
mtext(expression(paste(d," (cM)")),side=1,line=2)
mtext(expression(paste(alpha," (%)")),side=2,line=2)                    
lines(d,opt.alpha1)
lines(d,opt.alpha01)
lines(d,opt.alpha001)
lines(d,opt.alpha0001)
text(40,0.08,"c=1")
text(35,0.25,"c=0.1")
text(30,0.58,"c=0.01")
text(25,0.84,"c=0.001")
#dev.copy2eps(file="opt-alpha.eps")

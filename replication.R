#######################################################
#              Figure 7
#######################################################

i <- 0:10
plot(i*(i+1),i+1,type="s",axes=F,ylim=c(1,10),
     xlab=expression(tau^2/c),ylab=expression(t),cex.lab=1.4)
axis(2,at=i+1,cex.axis=1.4)
axis(1,at=i*(i+1),cex.axis=1.4)


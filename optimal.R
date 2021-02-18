recomb <- function(d)
  {
    0.5*(1-exp(-2*d))
  }

genetic.dist <- function(theta)
  {
    -0.5*log(1-2*theta)
  }

info <- function(alpha,theta)
  {
    theta1 <- recomb(0.5*genetic.dist(theta))
    q <- ((1-theta1)^2)/(theta1^2 + (1-theta1)^2)
    A <- (1-4*q*(1-q))*(1-theta)
    z <- -qnorm(alpha/2)
    i <- 2*A*( z*dnorm(z) + alpha/2 )
    i
  }


info0 <- function(alpha)
  {
    z <- -qnorm(alpha/2)
    2*( z*dnorm(z) + alpha/2 )
  }


info2cost <- function(alpha,theta,cost)
  {
    info(alpha,theta)/(1+cost*alpha)
  }

info2cost2 <- function(alpha,d,G,cost)
  {
    # alpha = selection fraction
    # d = distance between two markers in Morgans
    # G = genome size in Morgans
    # cost = cost of genotyping a single marker
    theta <- recomb(d)
    info(alpha,theta)/(1+cost*alpha*G/d)
  }

optalpha <- function(cost,theta)
  {
    optimize(f=info2cost,interval=c(0.0001,0.9999),max=T,
             theta=theta,cost=cost)$maximum
  }

optalpha2 <- function(d,G,cost)
  {
    optimize(f=info2cost2,interval=c(0.0001,0.9999),max=T,
             G=G,d=d,cost=cost)$maximum
  }

########################################################
#                      Figure 4
########################################################
# plot the optimal selection fraction for given cost when delta=0
# cost is defined as cost of complete genotyping to rearing cost
c <- 10^(-seq(-3,3,by=0.025))
a <- c
for( i in 1:length(c) )
  a[i] <- optalpha(c[i],0)
plot(log10(c),a,pch=".",axes=F,ylim=c(0,1),
     ylab="",xlab="Cost, c",yaxs="i",type="l")
mtext(side=2,line=2,expression(paste("Optimal selection fraction, ",alpha,
    ", (%)")))
axis(1,at=-3:3,labels=as.character(10^(-3:3)))
axis(2,at=seq(0,1,by=0.1),labels=as.character(100*seq(0,1,by=0.1)))
# dev.copy2eps(file="opt.sel.frac.eps")




########################################################
#                      Figure 5
########################################################

a <- seq(0.01,1,len=1000)
d <- seq(0.01,1,len=1000)
z001 <- outer(a,d,info2cost2,G=14.5,cost=0.001)
z01 <- outer(a,d,info2cost2,G=14.5,cost=0.01)
z1 <- outer(a,d,info2cost2,G=14.5,cost=0.1)
z <- outer(a,d,info2cost2,G=14.5,cost=1)
z10 <- outer(a,d,info2cost2,G=14.5,cost=10)

at.y <- 100*seq(0,1,by=0.1)
at.x <- 100*seq(0,2,by=0.1)

layout(matrix(c(1,2,3,4),nr=2,byrow=T))
par(mgp=c(2,1,0))
# contlabels <- c("1/2","1/4","1/8","1/16","1/32","1/64","1/128")
# contlevels <- -(1:7)

par(mar=c(0,3,4,0)+0.1)
contour(z/max(z),x=100*a,y=100*d,axes=F,frame.plot=T)
axis(2, at = at.y, labels = formatC(at.y, format="fg"))
axis(3, at = at.x, labels = formatC(at.x, format="fg"))
mtext(expression(paste(d," (cM)")),side=2,line=2)
mtext(expression(paste(alpha," (%)")),side=3,line=2)

par(mar=c(0,0,4,3)+0.1)
contour(z1/max(z1),x=100*a,y=100*d,axes=F,frame.plot=T)
axis(4, at = at.y, labels = formatC(at.y, format="fg"))
axis(3, at = at.x, labels = formatC(at.x, format="fg"))
mtext(expression(paste(d," (cM)")),side=4,line=2)
mtext(expression(paste(alpha," (%)")),side=3,line=2)

par(mar=c(4,3,0,0)+0.1)
contour(z01/max(z01),x=100*a,y=100*d,axes=F,frame.plot=T)
axis(2, at = at.y, labels = formatC(at.y, format="fg"))
axis(1, at = at.x, labels = formatC(at.x, format="fg"))
mtext(expression(paste(d," (cM)")),side=2,line=2)
mtext(expression(paste(alpha," (%)")),side=1,line=2)

par(mar=c(4,0,0,3)+0.1)
contour(z001/max(z001),x=100*a,y=100*d,axes=F,frame.plot=T)
axis(4, at = at.y, labels = formatC(at.y, format="fg"))
axis(1, at = at.x, labels = formatC(at.x, format="fg"))
mtext(expression(paste(d," (cM)")),side=4,line=2)
mtext(expression(paste(alpha," (%)")),side=1,line=2)

# dev.copy2eps(file="info2cost.ratio.eps")

#########################################################
#                Figure 3 
#########################################################
pmixnorm <- function(x,mm=c(0,0),ss=c(1,1),alpha=0.5,level=0)
{
  alpha * pnorm(x,mm[1],ss[1]) + (1-alpha) * pnorm(x,mm[2],ss[2]) - level
}

missinfo0 <- function(delta,alpha)
{
  limit <- uniroot(pmixnorm,interval=c(0,delta+5),mm=c(-delta,delta),
               level=1-alpha/2)$root
  mi <- 2*integrate(fracmiss0,rel.tol=1e-7,lower=0,upper=limit,
                    delta=delta)$value
  list( mi=mi, lim=limit )
}

missinfo1 <- function(delta,alpha,theta)
{
  limit <- uniroot(pmixnorm,interval=c(0,delta+5),mm=c(-delta,delta),
               level=1-alpha/2)$root
  mi.mid <- 2*integrate(fracmiss0,sub=10000,rel.tol=1e-7,lower=0,upper=limit,
                         delta=delta)$value
  mi0.upp <- 2*integrate(fracmiss1,sub=10000,rel.tol=1e-7,lower=limit,
                         upper=delta+10,
                         m=0,delta=delta,theta=theta)$value
  mi1.upp <- 2*integrate(fracmiss1,sub=10000,rel.tol=1e-7,lower=limit,
                         upper=delta+10,
                         m=1,delta=delta,theta=theta)$value

  list( mi=mi.mid+mi0.upp+mi1.upp, lim=limit )
}

missinfo2 <- function(delta,alpha,theta1,theta2)
{
  limit <- uniroot(pmixnorm,interval=c(0,delta+5),mm=c(-delta,delta),
               level=1-alpha/2)$root
  mi.mid <- 2*integrate(fracmiss0,sub=10000,rel.tol=1e-7,lower=0,upper=limit,
                         delta=delta)$value

  mi00.upp <- 2*integrate(fracmiss2,sub=10000,rel.tol=1e-7,lower=limit,
                          upper=delta+10,m1=0,m2=0,delta=delta,
                          theta1=theta1,theta2=theta2)$value

  mi01.upp <- 2*integrate(fracmiss2,sub=10000,rel.tol=1e-7,lower=limit,
                          upper=delta+10,m1=0,m2=1,delta=delta,
                          theta1=theta1,theta2=theta2)$value

  mi10.upp <- 2*integrate(fracmiss2,sub=10000,rel.tol=1e-7,lower=limit,
                          upper=delta+10,m1=1,m2=0,delta=delta,
                          theta1=theta1,theta2=theta2)$value

  mi11.upp <- 2*integrate(fracmiss2,sub=10000,rel.tol=1e-7,lower=limit,
                          upper=delta+10,m1=1,m2=1,delta=delta,
                          theta1=theta1,theta2=theta2)$value


  list( mi=mi.mid+mi00.upp+mi01.upp+mi10.upp+mi11.upp, lim=limit )
}


# ---------------------------------------------------------------------
# fraction of missing info when there is no genotype information
# ---------------------------------------------------------------------

fracmiss0 <- function(y,delta)
{
  density <- 0.5 * ( dnorm(y,mean=delta) + dnorm(y,mean=-delta) )
  A <- exp(delta*y)
  B <- exp(-delta*y)
  genovar <- (A*B)/((A+B)^2)
  # print(c(density,genovar))
  4*y*y*density*genovar
}

# ---------------------------------------------------------------------
# fraction of missing info when there is complete genotype at nearby marker
# ---------------------------------------------------------------------

fracmiss1 <- function(y,m,delta,theta)
{
  if(m==0)
    {
      density <- 0.5 * ( theta * dnorm(y,mean=delta)
                        + (1-theta) * dnorm(y,mean=-delta) )
      A <- exp(delta*y) * theta
      B <- exp(-delta*y) * (1-theta)
    }
  if(m==1)
    {
      density <- 0.5 * ( (1-theta) * dnorm(y,mean=delta)
                        + theta * dnorm(y,mean=-delta) )
      A <- exp(delta*y) * (1-theta)
      B <- exp(-delta*y) * theta
    }

  genovar <- (A*B)/((A+B)^2)
  # print(c(density,genovar))
  4*y*y*density*genovar
}

# ---------------------------------------------------------------------
# fraction of missing info when there is complete genotype at flanking markers
# ---------------------------------------------------------------------

fracmiss2 <- function(y,m1,m2,delta,theta1,theta2)
{
  # y = phenotype
  # m1 = left marker genotype
  # m2 = right marker genotype
  # delta = QTL effect; means are -delta and +delta
  # theta1 = recombination fraction to left marker
  # theta2 = recombination fraction to right marker

  lik0.pheno <- dnorm(y,mean=-delta)
  lik0.m1 <- theta1*m1 + (1-theta1)*(1-m1)
  lik0.m2 <- theta2*m2 + (1-theta2)*(1-m2)      
  lik0 <- lik0.pheno * lik0.m1 * lik0.m2 * 0.5

  lik1.pheno <- dnorm(y,mean=delta)
  lik1.m1 <- (1-theta1)*m1 + theta1*(1-m1)
  lik1.m2 <- (1-theta2)*m2 + theta2*(1-m2)      
  lik1 <- lik1.pheno * lik1.m1 * lik1.m2 * 0.5

  density <- lik0 + lik1
  genovar <- (lik0*lik1)/((lik0+lik1)^2)
  # print(c(density,genovar))
  4*y*y*density*genovar
}



# ---------------------------------------------------------------------

simpson <- function(f,lower,upper,nint=1000,...)
  {
    x <- seq(lower,upper,length=nint+1)
    ff <- match.fun(f)
    y <- ff(x,...)
    mult <- ( ( (0:nint) %% 2 ) + 1 ) * 2
    mult[1] <- 1
    mult[nint+1] <- 1
    h <- (upper-lower)/nint
    result <- sum(mult*y)*h/3
    result
  }


plotmissinfo2 <- function(delta,alpha)
{
n1 <- length(delta)
n2 <- length(alpha)
aaa <- matrix(nr=n1,nc=n2)

for(i1 in 1:n1 )
  {
   for(i2 in 1:n2 )
    {
    d <- delta[i1]
    a <- alpha[i2]
    # print(c(d,a))
    x <- missinfo2(d,a)	
    aaa[i1,i2] <- x$mi
    } 	
   }
list(mi=aaa,delta=delta,alpha=alpha)
}


plotmissinfo3 <- function(delta,alpha,theta)
{
n1 <- length(delta)
n2 <- length(alpha)
aaa <- matrix(nr=n1,nc=n2)

for(i1 in 1:n1 )
  {
   for(i2 in 1:n2 )
    {
    d <- delta[i1]
    a <- alpha[i2]
    # print(c(d,a))
    x <- missinfo3(d,a,theta)	
    aaa[i1,i2] <- x$mi
    } 	
   }
list(mi=aaa,delta=delta,alpha=alpha,theta)
}


plotmissinfo4 <- function(delta,alpha,theta1,theta2)
{
n1 <- length(delta)
n2 <- length(alpha)
aaa <- matrix(nr=n1,nc=n2)

for(i1 in 1:n1 )
  {
   for(i2 in 1:n2 )
    {
    d <- delta[i1]
    a <- alpha[i2]
    # print(c(d,a))
    x <- missinfo2(d,a,theta1,theta2)	
    aaa[i1,i2] <- x$mi
    } 	
   }
list(mi=aaa,delta=delta,alpha=alpha,theta1,theta2)
}




a2 <- plotmissinfo4(seq(0,2,len=321),seq(0.001,0.999,len=321),
                    theta1=0.2,theta2=0.2)
a1 <- plotmissinfo4(seq(0,2,len=321),seq(0.001,0.999,len=321),
                    theta1=0.1,theta2=0.1)
a05 <- plotmissinfo4(seq(0,2,len=321),seq(0.001,0.999,len=321),
                     theta1=0.05,theta2=0.05)
a01 <- plotmissinfo4(seq(0,2,len=321),seq(0.001,0.999,len=321),
                     theta1=0.01,theta2=0.01)

####################################
at.y <- seq(0,1,by=0.2)
at.x <- seq(0,2,by=0.2)

layout(matrix(c(1,2,3,4),nr=2,byrow=T))
par(mgp=c(2,1,0))
contlabels <- c("1/2","1/4","1/8","1/16","1/32","1/64","1/128")
contlevels <- -(1:7)

par(mar=c(0,3,4,0)+0.1)
contour(log(a2$mi)/log(2),x=a1$delta,y=a1$alpha,levels=contlevels,
        labels=contlabels,axes=F,frame.plot=T)
axis(2, at = at.y, labels = formatC(at.y, format="fg"))
axis(3, at = at.x, labels = formatC(at.x, format="fg"))
mtext(expression(alpha),side=2,line=2)
mtext(expression(delta),side=3,line=2)
mtext(expression(theta[1]==theta[2]==0.2),side=3,line=3)

par(mar=c(0,0,4,3)+0.1)
contour(log(a1$mi)/log(2),x=a1$delta,y=a1$alpha,levels=contlevels,
        labels=contlabels,axes=F,frame.plot=T)
axis(4, at = at.y, labels = formatC(at.y, format="fg"))
axis(3, at = at.x, labels = formatC(at.x, format="fg"))
mtext(expression(alpha),side=4,line=2)
mtext(expression(delta),side=3,line=2)
mtext(expression(theta[1]==theta[2]==0.1),side=3,line=3)

par(mar=c(4,3,0,0)+0.1)
contour(log(a05$mi)/log(2),x=a1$delta,y=a1$alpha,levels=contlevels,
        labels=contlabels,axes=F,frame.plot=T)
axis(2, at = at.y, labels = formatC(at.y, format="fg"))
axis(1, at = at.x, labels = formatC(at.x, format="fg"))
mtext(expression(alpha),side=2,line=2)
mtext(expression(delta),side=1,line=2)
mtext(expression(theta[1]==theta[2]==0.05),side=1,line=3)

par(mar=c(4,0,0,3)+0.1)
contour(log(a01$mi)/log(2),x=a1$delta,y=a1$alpha,levels=contlevels,
        labels=contlabels,axes=F,frame.plot=T)
axis(4, at = at.y, labels = formatC(at.y, format="fg"))
axis(1, at = at.x, labels = formatC(at.x, format="fg"))
mtext(expression(alpha),side=4,line=2)
mtext(expression(delta),side=1,line=2)
mtext(expression(theta[1]==theta[2]==0.01),side=1,line=3)

# dev.print(file="missinfo.eps",height=8,width=8)

########################################################
#          Figure 8
########################################################

# function for missing information for first QTL as a function of the
# phenotype y, and the effect of the second QTL, b

H <- function(y,b){
Im <- (y^2+b^2) - 2*b*y*tanh(b*y)
Im
}

# cdf of a mixture of normal distributions with means, mm, standard
# deviations, ss, mixture proportion, alpha, at the point, x; we can
# subtract a number, level, from it (this is for equation solving)

pmixnorm <- function(x,mm=c(0,0),ss=c(1,1),alpha=0.5,level=0)
{
  alpha * pnorm(x,mm[1],ss[1]) + (1-alpha) * pnorm(x,mm[2],ss[2]) - level
}

# density of mixture of normal distributions

dmixnorm <- function(x,mm=c(0,0),ss=c(1,1),alpha=0.5)
{
  alpha * dnorm(x,mm[1],ss[1]) + (1-alpha) * dnorm(x,mm[2],ss[2]) 
}

# quantiles of mixture of normal distributions

qmixnorm <- function(q,mm=c(0,0),ss=c(1,1),alpha=0.5)
  {
    limit <- uniroot(pmixnorm,interval=c(mm[1]-5,mm[2]+5),mm=mm,
                     level=q)$root
    limit
  }

# function to be integrated which is the product of the missing
# information function, H, and the density of the mixture of normal
# distributions, dmixnorm

f <- function(y,b)
  {
    Im <- H(y,b)
    Im * dmixnorm(y,mm=c(-b,b))
  }

# information as a function of the selection fraction, alpha, and the
# effect of the second QTL, b

info <- function(alpha,b)
  {
    limit <- qmixnorm(alpha/2,mm=c(-b,b))
    1-integrate(f,-abs(limit),abs(limit),b=b)$value
  }


# make figure

# make the selection fractions
alpha <- seq(1:999)/1000

# initialize the information vectors
i0 <- alpha
i0.5 <- alpha
i1 <- alpha
i1.5 <- alpha
i2 <- alpha
i4 <- alpha
i8 <- alpha
i10 <- alpha

# calculate the information functions
for( i in 1:999 )
  {
    i0[i] <- info(alpha[i],b=0)
    i0.5[i] <- info(alpha[i],b=0.5)
    i1[i] <- info(alpha[i],b=1)
    i1.5[i] <- info(alpha[i],b=1.5)    
    i2[i] <- info(alpha[i],b=2)
    i4[i] <- info(alpha[i],b=4)
    i8[i] <- info(alpha[i],b=8)
    i10[i] <- info(alpha[i],b=10)    
  }

plot(100*alpha,100*(1-i0),type="n",xlab="",ylab="",ylim=c(0,100),
     xlim=c(0,100),axes=F,yaxs="i")
lines(100*alpha,100*(1-i0))
lines(100*alpha,100*(1-i0.5),lty=2)
lines(100*alpha,100*(1-i1),lty=2)
lines(100*alpha,100*(1-i1.5),lty=2)
lines(100*alpha,100*(1-i8))
byten <- seq(0,100,by=10)
axis(1,at=byten,labels=as.character(byten))
axis(2,at=byten,labels=as.character(byten))
axis(4,at=byten,labels=as.character(byten))

mtext("Fraction of missing information (%)",side=2,line=2)
mtext(expression(paste("Selection fraction,",alpha," (%)")),side=1,line=2)
text(80,49,expression(paste(beta,"=",infinity)))
#text(70,47,expression(paste(beta,"=2")))
text(68,40,expression(paste(beta,"=3/2")))
text(50,34,expression(paste(beta,"=1")))
text(55,20,expression(paste(beta,"=1/2")))
text(36,11,expression(paste(beta,"=0")))
#dev.copy2eps(file="2qtl.eps")


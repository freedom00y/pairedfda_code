# Calcute the propotion of outliers in 
# normal,t(2),t(5),t(10),slash(1),slash(2),slash(5)
# with variance 0.04 and 0.25
rt_def <- function(n,df,sd)
{
  u = rgamma(n,df/2,df/2)
  sapply(u,function(x) rnorm(1,mean=0,sd = sd/sqrt(x)))
}

rslash_def <- function(n,df,sd)
{
  u = rbeta(n,df,1)
  sapply(u,function(x) rnorm(1,mean=0,sd = sd/sqrt(x)))
}

set.seed(1234)
n = 10000
a = c() 
for(sd in c(0.2,0.5))
{
  nsamp = rnorm(n,mean=0,sd=sd)
  nbox = boxplot(nsamp,plot=F)
  print(paste0("Outlier% for N(0,",sd,"^2) is ",length(nbox$out)/n))
  a = c(a,length(nbox$out)/n)
  for(t in c(2,5,10))
  {
    tsamp = rt_def(n,df=t,sd=0.2)
    tbox = boxplot(tsamp,plot=F)
    print(paste0("Outlier% for T(0,",t,",",sd,"^2) is ",length(tbox$out)/n))
    a = c(a,length(tbox$out)/n)
  }
  
  for(s in c(1,2,5))
  {
    ssamp = rt_def(n,df=s,sd=0.2)
    sbox = boxplot(ssamp,plot=F)
    print(paste0 ("Outlier% for Slash(0,",s,",",sd,"^2) is ",length(sbox$out)/n))
    a = c(a,length(sbox$out)/n)
  }
}

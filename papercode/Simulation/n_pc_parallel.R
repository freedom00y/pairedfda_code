library(pairedfda)
load("astronomy.RData")
outpath='~/Desktop/papercode4.0/real_analysis_ri3/'

bknots <<- c(c(-0.1,-0.05),
             seq(0, 1, length.out = 15),
             c(1.05,1.1))
bknots = 50*bknots-5

data=predata(nobs_y = nobs_r,
             nobs_z = nobs_i,
             time_y = R$Phase,
             time_z = II$Phase,
             y = R$Mag,
             z = II$Mag,
             knots = bknots, order = 3)

kakb<-function(data,ka,kb,type)
{
  source("real_analysis_ri3/Kfold_CV.R")
  lambda = simplex(data,10,ka,kb,type,tol = 1e-4,maxiter = 200)
  lambda = round(lambda,6)
  value = Kf_CV(data,lambda,K=10,ka,kb,type)
  return(list(value=value,lambda=lambda))
}

library(parallel)
n  = data$n
npc = 4
#for(type in c('t','n'))
for(type in c('n'))
{
  cvtable=matrix(nrow = npc,ncol = npc)
  record_lambda=matrix(nrow=npc*npc,ncol=6)
  flag=1
  print("start")
  no_cores <-detectCores()-1
  cl <- makeCluster(no_cores)
  clusterExport(cl,c("data","type","kakb","npc"))
  res_collect=parSapply(cl, 
                1:npc**2, 
                function(i){
                  library(pairedfda)
                  ka=ceiling(i/npc)
                  kb=(i-1)%%npc+1
                  res=try(kakb(data,ka,kb,type))
                  if(class(res)=='try-error')
                    rep('NA',8)
                  else
                    c(ka,kb,res$lambda,res$value)
                })
  filename=paste0(outpath,"cvtable_",type,".RData")
  save(res_collect,file = filename)
}

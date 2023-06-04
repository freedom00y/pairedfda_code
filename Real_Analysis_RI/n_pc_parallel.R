library(pairedfda)
load("astronomy.RData")

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
  source("Kfold_CV.R")
  tt = system.time({
  lambda = simplex(data,10,ka,kb,type,tol = 1e-5,maxiter = 400)
  lambda = round(lambda,6)
  value = Kf_CV(data,lambda,K=10,ka,kb,type)
  })
  return(list(value=value,lambda=lambda,time = tt[3]))
}

set.seed(1234)
library(parallel)
n  = data$n
npc = 4

no_cores <- 16
cl <- makeCluster(no_cores)

for(type in c("t","n"))
{
  cvtable=matrix(nrow = npc,ncol = npc)
  record_lambda=matrix(nrow=npc*npc,ncol=6)

  clusterExport(cl,c("data","type","kakb","npc"))
  res_collect=parSapply(cl, 
                1:npc**2, 
                function(i){
                  library(pairedfda)
                  ka=ceiling(i/npc)
                  kb=(i-1)%%npc+1
                  res=try(kakb(data,ka,kb,type))
                  if(class(res)=='try-error')
                    rep('NA',9)
                  else
                    c(ka,kb,res$lambda,res$value,res$time)
                })
  
  filename=paste0("cvtable_",type,".RData")
  save(res_collect,file = filename)
}
stopCluster(cl)


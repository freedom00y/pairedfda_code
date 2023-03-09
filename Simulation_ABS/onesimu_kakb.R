# select number of pc for simulation data

simu_kakb <- function(seed,type,freedom,varres)
{
  library(pairedfda)
  set.seed(seed)
  
  ## set hyperparameters
  n = 100
  
  ## generate data
  data = gen_data(n,varres,freedom,type,ka=2,kb=2)
  colnames(data$dataset)  = c("subject","obs_time","y_t","z_t")
  tmp.data = predata(nobs_y = data$obs_times,
                     nobs_z = data$obs_times,
                     time_y = data$dataset[,2],
                     time_z = data$dataset[,2],
                     y      = data$dataset[,3],
                     z      = data$dataset[,4],
                     knots  = 10,
                     order  = 3)
  tmp.data$score = data$score
  tmp.data$latent = data$latent
  data=tmp.data
  
  print("ka,kb,lambda1to4,CVV")
  ## different ka,kb 
  kakb<-function(data,ka,kb,type)
  {
    lambda = simplex(data,5,ka,kb,type,tol = 1e-4,maxiter = 100)
    #lambda=c(0,0,0,0)
    value = Kf_CV(data,lambda,K=5,ka,kb,type) #5-fold CV
    print(c(ka,kb,lambda,value))
    return(list(value=value,lambda=lambda))
  }
  
  best_ka=c()
  best_kb=c()
  for(type in c('t'))
  {
    res = rep(0,9)
    flag =1
    for(ka in 1:3)
      for(kb in 1:3)
      {
        res1 = try(kakb(data,ka,kb,type))
        if(class(res1)=='try-error')
          res[flag]=NA
        else
          res[flag] = kakb(data,ka,kb,type)$value
        flag = flag+1
      }
    # ind = which.min(res)
    # best_ka = c(best_ka, floor((ind-1)/3)+1)
    # best_kb = c(best_kb, (ind-1)%%3+1)
  }
  return(res)
}

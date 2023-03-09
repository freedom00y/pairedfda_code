# "seed" and "type" are used for generating data
simu <- function(seed,type,freedom,varres,percent,filename=NA)
{
  #################### Set Parameters ######################################
  
  ka = 2
  kb = 2
  n  = 100
  Kfold=5 ## The fold number for choosing tuning parameters
  
  #################### Source packages ######################################
  library(pairedfda)
  #source('IMSE_each.R') 
  source('IMSE_each_abs.R')          ## Compute the integral between estimation curve and the true one
  source('plot_curve.R')
  Rcpp::sourceCpp('dataGen_CurveOutlier.cpp')
  #################### Data Analysis #########################################
  ## Generate Data
  set.seed(seed)
  data = gen_data_curveoutlier(n=100,varres=varres,percent,ka=2,kb=2)
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
  ## Choose tuning parameters
  print("start")
  tlambda = simplex(data,Kfold,ka,kb,'t',tol = 1e-4,maxiter = 100)
  print("finish t")
  print(tlambda)
  slambda = simplex(data,Kfold,ka,kb,'s',tol = 1e-4,maxiter = 100)
  print("finish s")
  print(slambda)
  nlambda = simplex(data,Kfold,ka,kb,'n',tol = 1e-4,maxiter = 100)
  print("finish n")
  print(nlambda)
  #tlambda = slambda=nlambda=c(0,0,0,0)

  ## Estimate Parameters
  pt      = minEM(data, tlambda, type='t',ka,kb, tol = 1e-4)
  ps      = minEM(data, slambda, type='s',ka,kb, tol = 1e-4)
  pn      = minEM(data, nlambda, type='n',ka,kb, tol = 1e-4)
  if(is.na(filename)==FALSE)
    plot_curve(filename,pt,ps,pn)
  
  ## Compute IMSE
  int_t   = IMSE_each(data,pt,n)
  int_s   = IMSE_each(data,ps,n)
  int_n   = IMSE_each(data,pn,n)
  res     = c(int_t,int_s,int_n)
  return(res)
}

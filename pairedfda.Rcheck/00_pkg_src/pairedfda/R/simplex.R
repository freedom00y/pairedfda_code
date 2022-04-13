#' Choosing tuning parameter
#' 
#' Use simplex method to choose tuning parameter. 
#' Here we use K-fold cross-validation to find the tuning parameter such that mean square error (MSE) reachs the minimum value.
#'
#' @param data Processed data. Use "predata" to preprocess the raw data first.
#' @param Kfold Number of folds
#' @param ka Number of pcs of the first reponse variable (Y)
#' @param kb Number of pcs of the second reponse variable (Z)
#' @param type Model type. 'n' means normal; 't' means student-t; 's' means slash. The default is normal model.
#' @param tol Tolerance of the simplex algorithm. The default one is 1e-3.
#' @param maxiter Maximum iteration time. The default is 50 times.
#' 
#' @return The function returns the tuning parameter which is a vector with 4 components. The tuning parameters shows as the following order, the first mean curve \eqn{(\lambda_\mu)}, the second mean curve \eqn{(\lambda_\nu)}, the first pcs \eqn{(\lambda_f)}, the second pcs \eqn{(\lambda_g)}.
#' 
#' @details We suppose the model is 
#' @details \deqn{Y_i = B_i \theta_\mu + B_i f \alpha_i + \epsilon_i, Z_i = B_i \theta_\nu + B_i g \beta_i + \xi_i,}
#' @details where \eqn{(\alpha_i, \beta_i)} and the residuals follow normal, t or slash distribution. It is meaning of parameter `type` in this function.
#' @details The penalty term is 
#' @details \deqn{\lambda_\mu \theta_\mu^T \Omega \theta_\mu + \lambda_\nu \theta_\nu^T \Omega \theta_\nu + \lambda_f \sum_{j=1}^ka \theta_{fj}^T \Omega \theta_{fj} + \lambda_g \sum_{j=1}^kb \theta_{gj}^T \Omega \theta_{gj} }
#' 
#' @export
#'
#' @examples
#' rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)
#' data = predata(nobs_y = rawdata$obs_times, 
#'                nobs_z = rawdata$obs_times, 
#'                time_y = rawdata$dataset[,2], 
#'                time_z = rawdata$dataset[,2], 
#'                y = rawdata$dataset[,3], 
#'                z = rawdata$dataset[,4], 
#'                knots = 10, 
#'                order=3)
#' lambda_n = simplex(data,Kfold=5,ka=2,kb=2,type='n')
#' lambda_t = simplex(data,Kfold=5,ka=2,kb=2,type='t')
#' lambda_s = simplex(data,Kfold=5,ka=2,kb=2,type='s')
simplex<-function(data,Kfold,ka,kb,type='n',tol=1e-3,maxiter=50)
{
  ## source('Kfolder_tuning.R')     ## Computer K fold cross validation value
  
  # Using simplex method to choose tuning parameters
  ## Set hyperparameters
  a     = 1
  toler = tol
  alpha = 1
  beta  = 0.5
  gamma = 2
  
  ## Initialize Design Matrix
  design = rbind(diag(rep(a,4)),0.25*rep(a,4))
  y=rep(0,5)
  for(i in 1:5)
  {
    #y[i]=Kf_CV(data,exp(design[i,]),Kfold,ka,kb,type)
    temp = try(Kf_CV(data,exp(design[i,]),Kfold,ka,kb,type))
    if(class(temp)=='try-error' | is.nan(temp))
    {
      y[i] = 100
    }else{
      y[i] = temp
    }
  }
  if(min(y)==100)
  {
    return(c(0,0,0,0))
  }
  ind = which.min(y)
  optimal = c(design[ind,], y[ind]) 
   
  ## Simplex Algorithm
  iter=1
  while(sd(y/mean(y)-1)>toler & iter<maxiter)
  {
    ind_h=(1:5)[y==max(y)]
    ind_l=(1:5)[y==min(y)]
    ord = order(y)
    ind_l=ord[1]
    ind_h=ord[5]
    barP = apply(design[-ind_h,],2,mean)
    barP[barP<0]=0
    starP = (1+alpha)*barP-alpha*design[ind_h,]
    starP[starP<0]=0
    #stary=Kf_CV(data,exp(starP),Kfold,ka,kb,type)
    temp=try(Kf_CV(data,exp(starP),Kfold,ka,kb,type))
    if(class(temp)=='try-error' | is.nan(temp))
    {
      stary = 100
    }else{
      stary = temp
    }
    
    if(stary<y[ind_l]){
      ssP=(1+gamma)*starP
      ssP[ssP<0]=0
      temp=try(Kf_CV(data,exp(ssP),Kfold,ka,kb,type))
      if(class(temp)=='try-error' | is.nan(temp))
      {
        # ind = which.min(y)
        # optimal=c(design[ind,],y[ind])
        # break
        ssy = 100
      }else{
        ssy=temp
      }
      #print(paste0("ssy=",ssy))
      if(ssy<y[ind_l]){
        design[ind_h,]=ssP
        y[ind_h]=ssy
      }else{
        design[ind_h,]=starP
        y[ind_h]=stary
      }## end if ssy
      
    }else{
      
      if(stary<=y[ord[4]]){
        design[ind_h,]=starP
        y[ind_h]=stary
      }else{
        if(stary<=y[ind_h]){
          design[ind_h,]=starP
          y[ind_h]=stary
        }## end if stary
        
        ssP=beta*design[ind_h,]+(1-beta)*barP
        ssP[ssP<0]=0
        temp=try(Kf_CV(data,exp(ssP),Kfold,ka,kb,type))
        if(class(temp)=='try-error'| is.nan(temp))
        {
          # ind = which.min(y)
          # optimal=c(design[ind,],y[ind])
          # break
          ssy = 100
        }else{
          ssy=temp
        }######
 
        if(ssy>y[ind_h]){
          # flag=1
          # record_design=design[ind_l,]
          # record_y = y[ind_l]
          design=(design+design[ind_l,])/2
          for(i in 1:5)
          {
            temp=try(Kf_CV(data,exp(design[i,]),Kfold,ka,kb,type))
            if(class(temp)=='try-error'| is.nan(temp))
            {
              # flag=0
              # break
              y[i]=100
            }else{
              y[i]=temp
            }
          }
          # if(flag==0)
          # {
          #   optimal=c(record_design,record_y)
          #   break
          # }
            
        }else{
          design[ind_h,]=ssP
          y[ind_h]=ssy
        }## end if ssy
        
      }## end if stary
    }
    # print(cbind(design,y))
    if(min(y)<optimal[5])
    {
      ind = which.min(y)
      optimal=c(design[ind,],y[ind])
    }
    iter=iter+1
  }
  # print(iter)
  return(exp(optimal[1:4]))
}

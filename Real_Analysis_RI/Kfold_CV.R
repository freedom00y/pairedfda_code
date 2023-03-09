#' K-fold Cross-Validation
#' @description Divide the dataset into K fold. Use K-1 folds to estimate the parameters and use the left one fold to evaluate the mean square error each time. 
#' 
#' @param data Processed data. Use predata to preprocess first.
#' @param lambda The penalty patameters, a vector with 4 components showing in the following order, the first mean curve \eqn{(\lambda_\mu)}, the second mean curve \eqn{(\lambda_\nu)}, the first pcs \eqn{(\lambda_f)}, the second pcs \eqn{(\lambda_g)}.
#' @param K Number of folds
#' @param ka Number of pcs for the first reponse variable (Y)
#' @param kb Number of pcs for the second reponse variable (Z)
#' @param type Model type. 'n' means normal distribution, 't' means student-t distribution, 's' means slash distribution
#'
#' @return The sum of total square error (TSE) of all folds
#' @export
#'
#' @examples
#' rawdata = gen_data(n=100,
#'                    varres=0.01, 
#'                    gama=2, 
#'                    type='t',
#'                    ka=2,
#'                    kb=2)
#' data = predata(nobs_y = rawdata$obs_times, 
#'                nobs_z = rawdata$obs_times, 
#'                time_y = rawdata$dataset[,2], 
#'                time_z = rawdata$dataset[,2], 
#'                y = rawdata$dataset[,3], 
#'                z = rawdata$dataset[,4], 
#'                knots = 10, 
#'                order=3)
#' lambda_nopen = c(0,0,0,0)
#' MSE_nopen = Kf_CV(data,lambda_nopen,K=5,ka=2,kb=2,'t')
#' lambda_pen = rep(0.618,4)
#' MSE_pen = Kf_CV(data,lambda_pen,K=5,ka=2,kb=2,'t')
Kf_CV<-function(data,lambda,K,ka,kb,type)
{
  ## Split the data set into K folders
  ## Use likelihood as K-fold cross-validation value
  n    = data$n
  ind_y  = c(0, cumsum(data$nobs_y))
  ind_z  = c(0, cumsum(data$nobs_z))
  y = data$y
  z = data$z
  ny = data$nobs_y
  nz = data$nobs_z
  sumn_y = sum(data$nobs_y)
  sumn_z = sum(data$nobs_z)
  
  ## label all observations 
  y_group = c(0,sumn_y)
  z_group = c(0,sumn_z)
  num_mat_y = matrix(nrow=n,ncol=K)
  num_mat_z = matrix(nrow=n,ncol=K)
  for(i in 1:n)
  {
    ind1 = ind_y[i]+1
    ind2 = ind_y[i+1]
    divide_res_y = divide_group(ny[i],K)
    y_group[ind1:ind2]=divide_res_y$group_label
    num_mat_y[i,] = divide_res_y$num_in_group
    ind1 = ind_z[i]+1
    ind2 = ind_z[i+1]
    divide_res_z = divide_group(nz[i],K)
    z_group[ind1:ind2]=divide_res_z$group_label
    temp1=as.vector(table(divide_res_z$group_label))
    num_mat_z[i,] = divide_res_z$num_in_group
    temp2 = divide_res_z$num_in_group
  }

  value=c(0,0)
  for(i in 1:K)
  {
    ## create training data
    ### find the index 
    y_train_ind = (y_group!=i)
    z_train_ind = (z_group!=i)
    
    y_train = y[y_train_ind]
    z_train = z[z_train_ind]
    time_y_train = data$time_y[y_train_ind]
    time_z_train = data$time_z[z_train_ind]
    ny_train = ny - num_mat_y[,i]
    nz_train = nz - num_mat_z[,i]
    By_train = data$B_y[y_train_ind,]
    Bz_train = data$B_z[z_train_ind,]
    
    train_data     = list( n         = n,
                           nobs_y    = ny_train,
                           nobs_z    = nz_train,
                           time_y    = time_y_train,
                           time_z    = time_z_train,
                           y         = y_train,
                           z         = z_train,
                           B_y       = By_train,
                           B_z       = Bz_train,
                           knots     = data$knots,
                           Omega     = data$Omega)
    
    para = minEM(train_data, lambda, type,ka,kb, tol = 1e-4,maxiter = 50)
    
    ## create testing data
    y_test_ind = (y_group==i)
    z_test_ind = (z_group==i)
    
    y_test = y[y_test_ind]
    z_test = z[z_test_ind]
    time_y_test = data$time_y[y_test_ind]
    time_z_test = data$time_z[z_test_ind]
    ny_test = num_mat_y[,i]
    nz_test = num_mat_z[,i]
    By_test = data$B_y[y_test_ind,]
    Bz_test = data$B_z[z_test_ind,]
    
    test_data = list( n         = n,
                      nobs_y    = ny_test,
                      nobs_z    = nz_test,
                      time_y    = time_y_test,
                      time_z    = time_z_test,
                      y         = y_test,
                      z         = z_test,
                      B_y       = By_test,
                      B_z       = Bz_test,
                      knots     = data$knots,
                      Omega     = data$Omega)
    
    ## Compute the Criteria Value
    #temp = TSE(test_data,para,ny,nz)
    temp = TAE(test_data,para,ny,nz)
    value = value + c(temp$sy2,temp$sz2)
  }
  
  return(value/n)
}

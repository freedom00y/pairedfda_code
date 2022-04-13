TSE <- function(data,para,ny,nz)
{
  ## Load dataset
  n = data$n
  y = data$y
  z = data$z
  By = data$B_y
  Bz = data$B_z
  nobs_y = data$nobs_y
  nobs_z = data$nobs_z
  
  ## Load parameters
  mu = para$theta_mu
  nu = para$theta_nu
  f  = para$theta_f
  g  = para$theta_g
  ka = ncol(f)
  kb = ncol(g)
  alpha = para$alpha
  beta  = para$beta 
  # Da = para$Da
  # Db = para$Db
  # C  = para$C
  # eps = para$sig_eps
  # xi  = para$sig_xi
  n_total_y = nrow(By)
  n_total_z = nrow(Bz)
  
  # kab=ka+kb
  # ind_ka = as.vector(1:ka)
  # ind_kb = as.vector((ka+1):kab)
  # sig_ab = matrix(nrow = kab, ncol = kab)
  # sig_ab[ind_ka, ind_ka] = Da
  # sig_ab[ind_ka, ind_kb] = C
  # sig_ab[ind_kb, ind_ka] = t(C)
  # sig_ab[ind_kb, ind_kb] = Db
  # eig = eigen(sig_ab)
  # inv_sig_ab = eig$vectors%*%diag(1/eig$values)%*%t(eig$vectors)
  # 
  
  ## Estimate alpha & beta for test curves
  Bf    = By%*%f
  Bg    = Bz%*%g
  res_y = y - By%*%mu
  res_z = z - Bz%*%nu
  # alpha = matrix(nrow = n, ncol = ka)
  # beta  = matrix(nrow = n, ncol = kb)
  # 
  # ind_y = c(0, cumsum(nobs_y))      
  # ind_z = c(0, cumsum(nobs_z)) 
  # for (i in 1 : n )
  # {
  #   ni_y   = nobs_y[i]
  #   ni_z   = nobs_z[i]
  #   ind1_y = ind_y[i] + 1
  #   ind2_y = ind_y[i+1]
  #   ind1_z = ind_z[i] + 1
  #   ind2_z = ind_z[i+1]
  #   Bfi  = as.matrix(Bf[ind1_y:ind2_y,])
  #   Bgi  = as.matrix(Bg[ind1_z:ind2_z,])
  #   inv_sig = diag( c( rep(1/eps,ni_y), rep(1/xi,ni_z) ) )
  #   ind3    = 1:ni_y;
  #   ind4    = (ni_y+1):(ni_y+ni_z);
  #   Bth              = matrix(0,ncol = (ka+kb), nrow = ni_y+ni_z)
  #   Bth[ind3,ind_ka] = Bfi
  #   Bth[ind4,ind_kb] = Bgi
  #   resid     = c(res_y[ind1_y:ind2_y],res_z[ind1_z:ind2_z])
  #   
  #   inner     = inv_sig_ab + t(Bth)%*%inv_sig%*%Bth
  #   eig=eigen(inner)
  #   inv_inner = eig$vectors%*%diag(1/eig$values)%*%t(eig$vectors)
  #   inv_sigi  = inv_sig-inv_sig%*%Bth%*%inv_inner%*%t(Bth)%*%inv_sig
  #   bar_ab    = sig_ab%*% t(Bth) %*%inv_sigi%*%resid
  #   alpha[i,] = bar_ab[ind_ka]
  #   beta[i,]  = bar_ab[ind_kb]
  # }
  # 
  
  ## Create Augmented matrix Alpha and Beta
  Alpha = matrix(ncol = ka, nrow = n_total_y)
  for(i in 1:ka)
  {
    Alpha[,i]=rep(alpha[,i],nobs_y)
  }
  
  Beta = matrix(ncol = kb, nrow = n_total_z)
  for(i in 1:kb)
  {
    Beta[,i]=rep(beta[,i],nobs_z)
  }
  
  Bmu  = By%*%mu
  aBf  = apply(Alpha*(By%*%f),1,sum)
  yres = Bmu+aBf-y 
  w_y  = rep(1/ny, nobs_y)
  sy2  = sum(yres^2 * w_y)
  
  Bnu  = Bz%*%nu
  bBg  = apply(Beta*(Bz%*%g),1,sum)
  zres = Bnu+bBg-z 
  w_z  = rep(1/nz, nobs_z)
  sz2  = sum(zres^2 * w_z)
  
  return(list(sy2 = sy2,
              sz2=sz2))
}
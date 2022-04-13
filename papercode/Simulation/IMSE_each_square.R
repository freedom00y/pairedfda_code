## IMSE compute the integrated difference 
## between the estimated function and the real functions
IMSE_each <- function(data,pt1,n)
{
  x  = seq(0,1,0.01)
  y  = orthbasis(x,knots=10)$B
  mu = mu_t(x)
  nu = nu_t(x) 
  f1 = f_y1(x)
  f2 = f_y2(x)
  g1 = f_z1(x)
  g2 = f_z2(x)
  int= rep(0,8)
  res_mu_m = y%*%pt1$theta_mu-mu
  res_nu_m = y%*%pt1$theta_nu-nu
  res_f1_m = y%*%pt1$theta_f[,1]-f1
  res_f1_p = y%*%pt1$theta_f[,1]+f1
  res_f2_m = y%*%pt1$theta_f[,2]-f2
  res_f2_p = y%*%pt1$theta_f[,2]+f2
  res_g1_m = y%*%pt1$theta_g[,1]-g1
  res_g1_p = y%*%pt1$theta_g[,1]+g1
  res_g2_m = y%*%pt1$theta_g[,2]-g2
  res_g2_p = y%*%pt1$theta_g[,2]+g2
  res_af1_m= y%*%pt1$theta_f[,1]%*%t(pt1$alpha[,1])-f1%*%t(data$score[,1])
  res_af1_p= y%*%pt1$theta_f[,1]%*%t(pt1$alpha[,1])+f1%*%t(data$score[,1])
  res_af2_m= y%*%pt1$theta_f[,2]%*%t(pt1$alpha[,2])-f2%*%t(data$score[,2])
  res_af2_p= y%*%pt1$theta_f[,2]%*%t(pt1$alpha[,2])+f2%*%t(data$score[,2])
  res_bg1_m= y%*%pt1$theta_g[,1]%*%t(pt1$beta[,1])-g1%*%t(data$score[,3])
  res_bg1_p= y%*%pt1$theta_g[,1]%*%t(pt1$beta[,1])+g1%*%t(data$score[,3])
  res_bg2_m= y%*%pt1$theta_g[,2]%*%t(pt1$beta[,2])-g2%*%t(data$score[,4])
  res_bg2_p= y%*%pt1$theta_g[,2]%*%t(pt1$beta[,2])+g2%*%t(data$score[,4])
  

  int[1] = mean(res_mu_m^2)
  int[2] = mean(res_nu_m^2)
  int[3] = min( mean(res_f1_p^2), mean(res_f1_m^2) )
  int[4] = min( mean(res_f2_p^2), mean(res_f2_m^2) )
  int[5] = min( mean(res_g1_p^2), mean(res_g1_m^2) )
  int[6] = min( mean(res_g2_p^2), mean(res_g2_m^2) )
  int[7] = mean((res_mu_m + res_af1_m + res_af2_m)^2)
  int[8] = mean((res_nu_m + res_bg1_m + res_bg2_m)^2)
  #plot all curves and the real curve
  # pdf("debug/all.pdf")
  # par(mfrow=c(4,4))
  # for(i in 1:100)
  # {
  #   plot(x,mu+cbind(f1,f2)%*%data$score[i,1:2],'l')
  #   lines(x,y%*%pt1$theta_mu+y%*%pt1$theta_f%*%pt1$alpha[i,],col='red')
  # }
  # dev.off()
  return(int)
}
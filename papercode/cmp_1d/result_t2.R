library(R.matlab)
Rcpp::sourceCpp('dataGen.cpp')
# True value
x  = seq(0,1,0.01)
mu = mu_t(x)
nu = nu_t(x) 
f1 = f_y1(x)
f2 = f_y2(x)
g1 = f_z1(x)
g2 = f_z2(x)

# Estimated value
result = matrix(nrow=100,ncol=8)
for(i in 1:100)
{
  res = readMat(paste0("data_t2/res",i,".mat"))
  data = readMat(paste0("data_t2/data",i,".mat"))
  score = data$score
  res_mu_m = mu - res$mu
  res_nu_m = nu - res$nu
  res_f1_p = f1 + res$f[,1]
  res_f1_m = f1 - res$f[,1]
  res_f2_p = f2 + res$f[,2]
  res_f2_m = f2 - res$f[,2]
  res_g1_p = g1 + res$g[,1]
  res_g1_m = g1 - res$g[,1]
  res_g2_p = g2 + res$g[,2]
  res_g2_m = g2 - res$g[,2]
  res_y_m= mu%*%rep(1,100) + f1%*%t(score[,1])  + f2%*%t(score[,2]) - res$yhat
  res_z_m= mu%*%rep(1,100) + g1%*%t(score[,3])  + g2%*%t(score[,4]) - res$zhat
  
  int = rep(0,8)
  int[1] = mean(res_mu_m^2)
  int[2] = mean(res_nu_m^2)
  int[3] = min( mean(res_f1_p^2), mean(res_f1_m^2) )
  int[4] = min( mean(res_f2_p^2), mean(res_f2_m^2) )
  int[5] = min( mean(res_g1_p^2), mean(res_g1_m^2) )
  int[6] = min( mean(res_g2_p^2), mean(res_g2_m^2) )
  int[7] = mean((res_y_m)^2)
  int[8] = mean((res_z_m)^2)
  
  result[i,] = int
}

seppca = 1000*apply(result,2,mean)


filename = "/Users/huiya/Documents/Research/pairedfda/code/papercode/Simulation/simu_output/t_0.04_550.txt"
A=read.table(filename,skip=4,nrow=550)
A=as.matrix(A)
B = A[order(A[,1]),]
jointpca = apply(B[1:100,1:8],2,mean)



 
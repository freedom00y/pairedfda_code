library(pairedfda)
library(R.matlab)
for(i in 1:100)
{
  set.seed(i+7)
  data = gen_data(n=100,varre=0.01,gama=2,type="t",ka=2,kb=2)
  writeMat(paste0("data_t2/data",i,".mat"),data=data$dataset,score = data$score)
}

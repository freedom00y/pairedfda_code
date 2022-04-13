varres=0.25
simu_times=530
type='s'
freedom_range=c(1,2)
filename=paste('simu_kakb_output/',type,'_',varres,'_',simu_times,'.txt', sep = "")
title=paste('Generate from',type, 'distribution.','Repeat',simu_times,'times.',
            'The variance of residual is',varres)
write.table(title,filename)
for(free in freedom_range)
{
  write.table(paste("== degree of freedom: ",free,"=="),filename,append=TRUE)
  library(parallel)
  no_cores <-detectCores()-1
  cl <- makeCluster(no_cores)
  clusterExport(cl,c("free","filename","varres","type"))
  res1=parSapply(cl, 
                 1:simu_times+7, 
                 function(i){
                   source('onesimu_kakb.R')
                   temp=try(simu_kakb(seed = i,type,free,varres))
                   if(class(temp)=='try-error')
                     a=rep(NA,9)
                   else
                     a=temp
                   
                   write.table( t( c(i,round(a,4)) ), filename,append = TRUE,row.names = F,col.names = F)
                   a
                 })
  stopCluster(cl)
  # res = res1[,!is.na(res1[1,])]
  # res=apply(res[,1:100],1,mean)
  # res=matrix(res,nrow = 2,byrow = TRUE)
  # colnames(res)=c('t','slash','norm')
  # rownames(res)=c('ka',"kb")
  # write.table(round(1000*res,2),filename,append=TRUE)
}
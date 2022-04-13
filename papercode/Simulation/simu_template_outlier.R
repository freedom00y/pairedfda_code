###############################################################################
###############################################################################
varres=0.04
simu_times=550
type='n'
freedom_range=c(1)
filename=paste('simu_output_outlier/',type,'_',varres,'_',simu_times,'.txt', sep = "")
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
                   source('onesimu_outlier.R')
                   temp=try(simu(seed = i,type,free,varres))
                   if(class(temp)=='try-error')
                     a=rep(NA,24)
                   else
                     a=temp
                   
                   write.table( t( c(i,round(a,4)) ), filename,append = TRUE,row.names = F,col.names = F)
                   a
                 })
  stopCluster(cl)
  res = res1[,!is.na(res1[1,])]
  res=apply(res[,1:500],1,mean)
  res=matrix(res,nrow = 3,byrow = TRUE)
  colnames(res)=c("Mean1","Mean2","YPC1","YPC2","ZPC1","ZPC2","IND1","IND2")
  rownames(res)=c('t',"slash","norm")
  write.table(round(1000*res,2),filename,append=TRUE)
}  



###############################################################################
###############################################################################
varres=0.25
simu_times=550
type='n'
freedom_range=c(1)
filename=paste('simu_output_outlier/',type,'_',varres,'_',simu_times,'.txt', sep = "")
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
                  source('onesimu_outlier.R')
                  temp=try(simu(seed = i,type,free,varres))
                  if(class(temp)=='try-error')
                    a=rep(NA,24)
                  else
                    a=temp
                  
                  write.table( t( c(i,round(a,4)) ), filename,append = TRUE,row.names = F,col.names = F)
                  a
                })
  stopCluster(cl)
  res = res1[,!is.na(res1[1,])]
  res=apply(res[,1:500],1,mean)
  res=matrix(res,nrow = 3,byrow = TRUE)
  colnames(res)=c("Mean1","Mean2","YPC1","YPC2","ZPC1","ZPC2","IND1","IND2")
  rownames(res)=c('t',"slash","norm")
  write.table(round(1000*res,2),filename,append=TRUE)
}  
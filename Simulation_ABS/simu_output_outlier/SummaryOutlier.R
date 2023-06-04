SumTable <- function(filename,n=550)
{
  res = read.table(paste0(filename,".txt"),skip=n+5,nrow=3)
  colnames(res) = c("Fit_Method","Mean1","Mean2",
                    "YPC1","YPC2","ZPC1","ZPC2","Ind1","Ind2")
  Gen_Method = filename
  res = cbind(Gen_Method,res)
  print(apply(res[,-c(1,2)],2,which.max))
  return(res)
}

files=c(#"n_0.04_550_0.01","n_0.25_550_0.01",
        #"n_0.04_550_0.03","n_0.25_550_0.03",
        "n_0.04_550_0.05","n_0.25_550_0.05")

AllRes = data.frame()
for(k in 1:2)
{
  filename = files[k]
  OneRes = SumTable(filename)
  AllRes = rbind(AllRes,OneRes)
}
write.csv(AllRes,file="PointOutlier.csv")

files = paste0("curveoutlier_",c(#"0.04_550_0.01","0.25_550_0.01",
                                #"0.04_550_0.03","0.25_550_0.03",
                                 "0.04_550_0.05","0.25_550_0.05"))
AllRes = data.frame()
for(k in 1:2)
{
  filename = files[k]
  OneRes = SumTable(filename,n=550)
  AllRes = rbind(AllRes,OneRes)
}
write.csv(AllRes,file="CurveOutlier.csv")

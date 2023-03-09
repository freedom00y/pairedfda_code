files=c("t_0.04_550","t_0.25_550",
        "s_0.04_550","s_0.25_550",
        "n_0.04_550","n_0.25_550")

SumTable <- function(filename)
{
  nsamp = 3
  method = substring(filename, 1, 1)
  if(method != "n")
  {
    A = read.table(paste0(filename,".txt"),skip=555,nrow=nsamp)
    B = read.table(paste0(filename,".txt"),skip=1111,nrow=nsamp)
    C = read.table(paste0(filename,".txt"),skip=1667,nrow=nsamp)
    res = rbind(A,B,C)
    print(apply(A[,-1],2,which.max))
    print(apply(B[,-1],2,which.max))
    print(apply(C[,-1],2,which.max))
  }else{
    res = read.table(paste0(filename,".txt"),skip=555,nrow=nsamp)
    print(apply(res[,-1],2,which.max))
  }
  
  colnames(res) = c("Fit_Method","Mean1","Mean2",
                    "YPC1","YPC2","ZPC1","ZPC2","Ind1","Ind2")
  Gen_Method = filename
  res = cbind(Gen_Method,res)
  return(res)
}

AllRes = data.frame()
for(k in 1:6)
{
  print(k)
  filename = files[k]
  OneRes = SumTable(filename)
  AllRes = rbind(AllRes,OneRes)
}
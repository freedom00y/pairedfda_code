mu_t <- function(t){
  return( 2.5 + 2.5*t + 2.5 * exp(-20* (t-0.6)**2) ) 
}

nu_t <- function(t){
  len = length(t)
  return( 7.5 - 1.5*t - 2.5 * exp(-20* (t-0.3)**2) ) 
}

fy1_t <- function(t){
  return( sqrt(15)/(1+sqrt(5)) * (t**2 + 1/sqrt(5)) )
}

fy2_t <- function(t){
  return( sqrt(15)/(sqrt(5)-1) * (t**2 - 1/sqrt(5)) )
}


fz1_t <- function(t){
  return(sqrt(2) * cos(2*pi*t))
}

fz2_t <- function(t){
  return(sqrt(2) * sin(2*pi*t))
}


ttime = 0.01*0:100
mu = mu_t(ttime)
nu = nu_t(ttime)
f1 = fy1_t(ttime)
f2 = fy2_t(ttime)
g1 = fz1_t(ttime)
g2 = fz2_t(ttime)


data2frame <- function(n=100,varres,gama,type,ka=2,kb=2,outlier=F)
{
  rawdata = gen_data(n=n,varres=varres, gama=gama, type=type,ka=ka,kb=kb)
  data = rawdata$dataset
  colnames(data) = c("sub","time","y","z")
  if(outlier == T)
  {
    Nobs=dim(data)[1]
    nout = round(0.03*Nobs)
    ind_cg_y = sample(1:Nobs,nout,replace = FALSE)
    data[ind_cg_y,3] = data[ind_cg_y,3]+
      sample(c(-1,1),nout,replace=TRUE)*runif(nout,8,10)
    ind_cg_z = sample(1:Nobs,nout,replace = FALSE)
    data[ind_cg_z,4] = data[ind_cg_z,4]+
      sample(c(-1,1),nout,replace=TRUE)*runif(nout,8,10)
  }
  data = as.data.frame(data)
  data["u"] = rep(rawdata$latent,rawdata$obs_times)
  data_y = data[,c(1,2,3,5)]
  data_y["variable"] = "Y"
  colnames(data_y)[3] = "value"
  data_z = data[,c(1,2,4,5)]
  data_z["variable"] = "Z"
  colnames(data_z)[3] = "value"
  data_yz = rbind(data_y,data_z)
  
  
  nt  = length(ttime)
  ty  = c(apply(rawdata$score,1,function(u) mu+u[1]*f1+u[2]*f2))
  tz  = c(apply(rawdata$score,1,function(u) nu+u[3]*g1+u[4]*g2))
  tdata_yz = data.frame(sub   = rep(1:n-1,each=length(ttime)),
                        time  = ttime,
                        value = c(ty,tz),
                        u     = rep(rawdata$latent,each=nt),
                        variable = rep(c("Y","Z"), each=n*nt))
  
  return(list(obs=data_yz,true=tdata_yz))
}



library(ggplot2)

nsample = 3
set.seed(1)
data_t  = data2frame(n=nsample,varres=0.04,gama=2,type="t",ka=2,kb=2)

set.seed(15)
data_s  = data2frame(n=nsample,varres=0.04,gama=1,type="s",ka=2,kb=2)

set.seed(20)
data_n  = data2frame(n=nsample,varres=0.04,gama=1,type="n",ka=2,kb=2)

set.seed(111)
data_nc = data2frame(n=nsample,varres=0.04,gama=1,type="n",ka=2,kb=2,outlier=T)

data = rbind(data_t$obs,data_s$obs,data_nc$obs,data_n$obs)
method = rep(c("RRME-t(2)","RRME-slash(1)","Contaminated RRME-normal","RRME-normal"),
             c(dim(data_t$obs)[1],dim(data_s$obs)[1],dim(data_nc$obs)[1],dim(data_n$obs)[1]))
method = factor(method,levels=c("RRME-t(2)","RRME-slash(1)","Contaminated RRME-normal","RRME-normal"))
data = cbind(data,method)
data['sub'] = factor(data$sub,levels=c("0","1","2"))

data0 = rbind(data_t$true,data_s$true,data_nc$true,data_n$true)
method = rep(c("RRME-t(2)","RRME-slash(1)","Contaminated RRME-normal","RRME-normal"),
             c(dim(data_t$true)[1],dim(data_s$true)[1],dim(data_nc$true)[1],dim(data_n$true)[1]))
method = factor(method,levels=c("RRME-t(2)","RRME-slash(1)","Contaminated RRME-normal","RRME-normal"))
data0 = cbind(data0,method)
data0['sub'] = factor(data0$sub,levels=c("0","1","2"))

cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


o1 = dim(data_t$obs)[1]+dim(data_s$obs)[1]+which.min(data_nc$obs$value)
o2 = dim(data_t$obs)[1]+dim(data_s$obs)[1]+which.max(data_nc$obs$value)
datasub = data[c(o1,o2),]

p = ggplot(data,aes(x=time,y=value,group=sub,color=sub)) + 
  geom_point(size=1.5) +
  geom_line(data=data0,aes(x=time,y=value,group=sub,color=sub),linetype = "dashed") + 
  geom_point(data = datasub,aes(x=time,y=value,group=sub,color=sub),size=3,shape=8) +
  xlab("t") + ylab("") + 
  facet_grid(cols=vars(variable),rows=vars(method)) + 
  theme_bw() +
  #theme_set(theme_bw()) +
  theme(panel.grid=element_blank(),legend.position = "none") +
  scale_colour_manual(values = cbPalette) 
p
ggsave(filename="simulation_data_withtrue.pdf",width = 5.7,height = 7.8)

# p = ggplot(data,aes(x=time,y=value,group=sub,color=sub)) + geom_line(linetype = "dashed") +
#   geom_point(size=1) +
#   xlab("t") + ylab("") + 
#   facet_grid(cols=vars(variable),rows=vars(method)) + 
#   theme_set(theme_bw()) +
#   theme(panel.grid=element_blank(),legend.position = "none") +
#   scale_colour_manual(values = cbPalette) 
# p
# ggsave(filename="simulation_data.pdf",width = 5.7,height = 7.5)
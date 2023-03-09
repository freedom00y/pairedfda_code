ka = 2
kb = 2
n  = 5
Kfold=5 ## The fold number for choosing tuning parameters

varres=0.25
freedom=2
type='t'
#################### Source packages ######################################
library(pairedfda)

#################### Data Analysis #########################################
## Generate Data
set.seed(202011)
data = gen_data(n,varres,freedom,type,ka,kb)
colnames(data$dataset)  = c("subject","obs_time","y_t","z_t")
Nobs=dim(data$dataset)[1]
nout = round(0.01*Nobs)
ind_cg_y = sample(1:Nobs,nout,replace = FALSE)
data$dataset[ind_cg_y,3] = data$dataset[ind_cg_y,3]+
  sample(c(-1,1),nout,replace=TRUE)*runif(nout,8,10)
ind_cg_z = sample(1:Nobs,nout,replace = FALSE)
data$dataset[ind_cg_z,4] = data$dataset[ind_cg_z,4]+
  sample(c(-1,1),nout,replace=TRUE)*runif(nout,8,10)

tmp.data = predata(nobs_y = data$obs_times,
                   nobs_z = data$obs_times,
                   time_y = data$dataset[,2],
                   time_z = data$dataset[,2],
                   y      = data$dataset[,3],
                   z      = data$dataset[,4],
                   knots  = 10,
                   order  = 3)
tmp.data$score = data$score
tmp.data$latent = data$latent
data=tmp.data

pdf("plot_outlier.pdf", width=9, height=5)
nobs_y = data$nobs_y
nobs_z = data$nobs_z
ind_y = c(0,cumsum(nobs_y))
ind_z = c(0,cumsum(nobs_z))
color = c("blue","blue","red","orange")
opar = par(no.readonly=TRUE)
par(mfrow=c(1,2))
ty = data$time_y
tz = data$time_z
y = data$y
z = data$z


i=2
indrange = (ind_y[i]+1):ind_y[i+1]
plot(ty[indrange],y[indrange],'o',col=color[1],ylim=range(y), xlab='Time', ylab='1st response variable',lty=2)
for(i in 3:4)
{
  indrange = (ind_y[i]+1):ind_y[i+1]
  lines(ty[indrange],y[indrange],'o',col=color[i],lty=2)
}
points(ty[ind_cg_y],y[ind_cg_y],pch=8)
legend("bottomleft",c("normal observations","artifical outliers"),pch=c(1,8),cex=0.8,bty = "n")

i=1
indrange = (ind_z[i]+1):ind_z[i+1]
plot(tz[indrange],z[indrange],'o',col=color[1],ylim=range(z), xlab='Time', ylab='2nd response variable',lty=2)
for(i in 3:4)
{
  indrange = (ind_z[i]+1):ind_z[i+1]
  lines(tz[indrange],z[indrange],'o',col=color[i],lty=2)
}
points(tz[ind_cg_z],z[ind_cg_z],pch=8)
on.exit(par(opar))
dev.off()

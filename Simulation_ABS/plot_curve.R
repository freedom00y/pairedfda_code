plot_curve<-function(string,pt,ps,pn)
{
  #Rcpp::sourceCpp('simu_funcs.cpp')
x = seq(0,1,0.01)
y = orthbasis(x)$B
pdf(paste(string,'.pdf',sep = ''))
split.screen(c(3,2))
#par(mfrow=c(3,2))

screen(1)
par(mar=c(5.1,4.1,0.1,0.1))
plot(x,y%*%pt$theta_mu,'l',col='red',xlab='',ylab='Mean',lty=3,ylim=c(2.5,6.5))
points(x,y%*%pn$theta_mu,'l',col='blue',lty=2)
points(x,y%*%ps$theta_mu,'l',col='green',lty=4)
curve(mu_t,xlim=c(0,1),add = TRUE,lty=1)
screen(2)
par(mar=c(5.1,4.1,0.1,0.5))
plot(x,y%*%pt$theta_nu,'l',col='red',xlab='',ylab='',lty=3,ylim=c(4,7.5))
points(x,y%*%pn$theta_nu,'l',col='blue',lty=2)
points(x,y%*%ps$theta_nu,'l',col='green',lty=4)
curve(nu_t,xlim=c(0,1),add = TRUE,lty=1)

legend("topright",c("true","norm","t","slash"),
       col=c("black","blue","red","green"),bty="n",
       text.col=c("black","blue","red","green"),
       lty=1:4,
       cex=0.5)


screen(3)
par(mar=c(5.1,4.1,0.1,0.1))
# plot(x,y%*%pt1$theta_f[,1],'l',col='red',xlab='',
#      ylab='PC function 1',lty=3,ylim=c(0,2))
plot(x,y%*%pt$theta_f[,1],'l',col='red',xlab='',ylab='PC function 1',lty=3,ylim=c(0,2))
points(x,y%*%pn$theta_f[,1],'l',col='blue',lty=2)
points(x,y%*%ps$theta_f[,1],'l',col='green',lty=4)
curve(f_y1,xlim=c(0,1),add = TRUE)

screen(4)
par(mar=c(5.1,4.1,0.1,0.1))
plot(x,y%*%pt$theta_g[,1],'l',col='red',xlab='',ylab='',lty=3)
points(x,y%*%pn$theta_g[,1],'l',col='blue',lty=2)
points(x,y%*%ps$theta_g[,1],'l',col='green',lty=4)
curve(f_z1,xlim=c(0,1),add = TRUE,lty=1)
screen(5)
par(mar=c(5.1,4.1,0.1,0.5))
# plot(x,-y%*%pt1$theta_f[,2],'l',col='red',xlab='Y',ylab='PC function 2',lty=3,ylim=c(-1.5,1.5))
plot(x,-y%*%pt$theta_f[,2],'l',col='red',xlab='Y',ylab='PC function 2',lty=3,ylim=c(-1.5,1.5))
points(x,-y%*%pn$theta_f[,2],'l',col='blue',lty=2,ylim=c(-2,2))
points(x,-y%*%ps$theta_f[,2],'l',col='green',lty=4)
curve(f_y2,xlim=c(0,1),add = TRUE,lty=1)

screen(6)
par(mar=c(5.1,4.1,0.1,0.5))
plot(x,y%*%pt$theta_g[,2],'l',col='red',xlab='Z',ylab='',lty=3,ylim=c(-1.5,1.5))
points(x,-y%*%pn$theta_g[,2],'l',col='blue',lty=2)
points(x,y%*%ps$theta_g[,2],'l',col='green',lty=4)
curve(f_z2,xlim=c(0,1),add = TRUE,lty=1)

close.screen(all=TRUE)
erase.screen() 
dev.off()
}
lw = 0
par(mfrow=c(1,2))
load("cvtable_t.RData")
t_table = res_collect
Snpc = as.double(res_collect[1,])+as.double(res_collect[2,])
Smae = as.double(res_collect[7,])+as.double(res_collect[8,])
plot(Snpc,Smae,pch=20,
     xlab="ka+kb",ylab = "MAE",xlim=c(2,8.7),ylim=c(0.10,0.18),
     main = "RRME-t model")                                     
#text(Snpc+0.4,Smae+0.001,paste0("(",res_collect[1,],",",res_collect[2,],")"))
text(Snpc-lw,Smae,paste0("(",res_collect[1,],",",res_collect[2,],")"),pos=4,cex=0.5)
lowind = c(1,2,3,7,14,15,16)
lines(Snpc[lowind],Smae[lowind])


load("cvtable_n.RData")
n_table = res_collect

Snpc = as.double(res_collect[1,])+as.double(res_collect[2,])
Smae = as.double(res_collect[7,])+as.double(res_collect[8,])
plot(Snpc,Smae,pch=20,
     xlab="ka+kb",ylab = "MAE",xlim=c(2,8.7),ylim=c(0.11,0.18),
     main = "RRME-normal model")                                     
#text(Snpc+0.4,Smae+0.001,paste0("(",res_collect[1,],",",res_collect[2,],")"))
text(Snpc-lw,Smae,paste0("(",res_collect[1,],",",res_collect[2,],")"),pos=4,cex=0.6)
lowind = c(1,2,3,7,11,12,16)
lines(Snpc[lowind],Smae[lowind])
legend("bottomright", c("(ka,kb)"),
       bg="grey90",
       pch = 20, ncol = 1,box.lwd=1,cex=0.6)

############################################################
library(ggplot2)
load("cvtable_t.RData")
t_table = res_collect
Snpc = as.double(res_collect[1,])+as.double(res_collect[2,])
Smae = as.double(res_collect[7,])+as.double(res_collect[8,])
name = paste0("(",res_collect[1,],",",res_collect[2,],")")
df_t = data.frame(Snpc = Snpc, Smae=Smae,name=name,method="RRME-t")
lowind_t = c(1,2,3,7,14,15,16)
load("cvtable_n.RData")
n_table = res_collect
Snpc = as.double(res_collect[1,])+as.double(res_collect[2,])
Smae = as.double(res_collect[7,])+as.double(res_collect[8,])
df_n = data.frame(Snpc = Snpc, Smae=Smae,name=name,method="RRME-normal")
lowind_n = c(1,2,3,7,11,12,16)
df = rbind(df_t,df_n)
df$method = factor(df$method,levels=c("RRME-t","RRME-normal"))
df_sub = rbind(df_t[lowind_t,],df_n[lowind_n,])
df_sub$method = factor(df_sub$method,levels=c("RRME-t","RRME-normal"))
p_kakb = ggplot(df,aes(x=Snpc,y=Smae)) + 
  geom_text(aes(label = name),label.size = NA,nudge_x = 0.4) +
  geom_point() + 
  geom_line(data=df_sub)+
  geom_label(aes(x=8,0.11),label=expression("(k"*alpha*",k"*beta*")"))+
  facet_grid(cols =vars(method))+
  xlab(expression("k"*alpha*"+k"*beta)) + ylab("MAE")
ggsave("select_kakb.pdf",p_kakb, width=9, height=4, units="in")

############################################################
# par(mfrow=c(1,2))
# load(paste0(outpath,"cvtable_t.RData"))
# t_table = res_collect
# Snpc = as.double(res_collect[1,])+as.double(res_collect[2,])
# Smae = as.double(res_collect[7,])+as.double(res_collect[8,])
# pos = paste0("(",res_collect[1,],",",res_collect[2,],")")
# ind = order(Smae,decreasing = T)
# plot(1:16,Smae[ind],"l",xlab="(ka,kb)",ylab = "MAE",main = "RRME-t model",xaxt="n")
# axis(1,at=1:16,pos[ind],font=0.5,las=2)                                     
# 
# 
# load(paste0(outpath,"cvtable_n.RData"))
# n_table = res_collect
# Snpc = as.double(res_collect[1,])+as.double(res_collect[2,])
# Smae = as.double(res_collect[7,])+as.double(res_collect[8,])
# pos = paste0("(",res_collect[1,],",",res_collect[2,],")")
# ind = order(Smae,decreasing = T)[-16]
# plot(1:15,Smae[ind],"l",xlab="(ka,kb)",ylab = "MAE",main = "RRME-n model",xaxt="n")
# axis(1,at=1:15,pos[ind],font=0.5,las=2)  
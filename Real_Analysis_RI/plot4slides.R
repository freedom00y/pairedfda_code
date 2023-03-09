library(pairedfda)
load("../astronomy.RData")

bknots <<- c(c(-0.1,-0.05),
             seq(0, 1, length.out = 15),
             c(1.05,1.1))
bknots = 50*bknots-5


data=predata(nobs_y = nobs_r,
             nobs_z = nobs_i,
             time_y = R$Phase,
             time_z = II$Phase,
             y = R$Mag,
             z = II$Mag,
             knots = bknots, order = 3)

filenames <- list.files(path="../ExtinctionCorrectedAlignedData/",pattern = ".*.txt")


library(ggplot2)
x = seq(-10,50,0.1)
Bval = orthbasis(x,bknots)$B
titles = c("SN2004eo","SN2005iq","SN2002fb","SN2005ke")
title_ind=c()
for(i in 1:length(titles))
{
  title_ind = c(title_ind,which(filenames==paste0(titles[i],".txt")))
}


Band = factor( rep( rep(c("R","I"),each=2*length(x)), 4) ,levels=c("R","I") )
method = rep( rep(c("T","Normal","T","Normal"), each = length(x)), 4) 

n = data$n
ind_y = c(0, cumsum(data$nobs_y))
ind_z = c(0, cumsum(data$nobs_z))
time_y=data$time_y
time_z=data$time_z

pt_x=c()
pt_y=c()
pt_band=c()
pt_sub=c()

flag=1
for(j in title_ind)
{
  temp_x = c(time_y[(ind_y[j]+1):ind_y[j+1]], time_z[(ind_z[j]+1):ind_z[j+1]])
  temp_y = c(data$y[(ind_y[j]+1):ind_y[j+1]], data$z[(ind_z[j]+1):ind_z[j+1]])
  temp_band = rep(c("R","I"), c(nobs_r[j],nobs_i[j]))
  temp_sub = rep(titles[flag],nobs_r[j]+nobs_i[j])
  flag = flag+1
  pt_x = c(pt_x, temp_x)
  pt_y = c(pt_y, temp_y)
  pt_band = c(pt_band, temp_band)
  pt_sub = c(pt_sub, temp_sub)
}
pt_sub=factor(pt_sub,levels=c("SN2004eo","SN2005iq","SN2002fb","SN2005ke"))
pt_band=factor(pt_band, levels = c("R","I") )
df_sub_pt = data.frame(pt_x,pt_y,Band = pt_band, sub=pt_sub)


p_sub = ggplot(df_sub_pt) + 
  geom_point(aes(x=pt_x, y=pt_y)) +
  facet_grid(vars(Band), vars(sub)) +
  ylim(3.4,-0.2) + labs(x = "Phase", y="Magnitude")


ggsave("test_sub_curve.pdf",p_sub, width=10, height=5, units="in")
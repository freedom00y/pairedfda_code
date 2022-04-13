library(R.matlab)
library(ggplot2)
fittedy = readMat("fitted_y.mat")
yhat = as.numeric(fittedy$xhatC2)
truey = readMat("data_t2/data1.mat")
data = apply(truey$data, 2, as.numeric)
data_plot = rbind(data[,1:3],
                  cbind(data[,1:2],yhat))
data_plot = as.data.frame(data_plot)
colnames(data_plot) = c("sub","time","y")
data_plot["label"] = rep(c("Fitted","True"),each=length(yhat))
ggplot(data_plot,aes(x=time,y=y,col=sub,group=sub)) + geom_line() +
  facet_grid(cols=vars(label))

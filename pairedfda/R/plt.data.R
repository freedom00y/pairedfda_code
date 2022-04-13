#' Plot the raw data
#' @description display two plots of the paired raw data
#' @param data Processed data. Use "predata" to preprocess the raw data first.
#'
#' @importFrom grDevices graphics
#' @export
#' @examples
#' rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)
#' data = predata(nobs_y = rawdata$obs_times, 
#'                nobs_z = rawdata$obs_times, 
#'                time_y = rawdata$dataset[,2], 
#'                time_z = rawdata$dataset[,2], 
#'                y = rawdata$dataset[,3], 
#'                z = rawdata$dataset[,4], 
#'                knots = 10, 
#'                order=3)
#' plt.data(data)
plt.data <- function(data)
{
  n = data$n
  nobs_y = data$nobs_y
  nobs_z = data$nobs_z
  ind_y = c(0,cumsum(nobs_y))
  ind_z = c(0,cumsum(nobs_z))
  #color = rainbow(n)
  color= heat.colors(n, alpha = 1)
  opar = par(no.readonly=TRUE)
  par(mfrow=c(1,2))
  ty = data$time_y
  tz = data$time_z
  y = data$y
  z = data$z
  
  indrange = (ind_y[1]+1):ind_y[2]
  plot(ty[indrange],y[indrange],'l',col=color[1],ylim=range(y), xlab='Time', ylab='1st Response Variable')
  for(i in 2:n)
  {
    indrange = (ind_y[i]+1):ind_y[i+1]
    lines(ty[indrange],y[indrange],'l',col=color[i])
  }
  
  indrange = (ind_z[1]+1):ind_z[2]
  plot(tz[indrange],z[indrange],'l',col=color[1],ylim=range(z), xlab='Time', ylab='2nd Response Variable')
  for(i in 2:n)
  {
    indrange = (ind_z[i]+1):ind_z[i+1]
    lines(tz[indrange],z[indrange],'l',col=color[i])
  }
  on.exit(par(opar))
}
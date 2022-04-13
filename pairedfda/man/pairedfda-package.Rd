\name{pairedfda-package}
\alias{pairedfda-package}
\alias{pairedfda}
\docType{package}
\title{
\packageTitle{pairedfda}
}
\description{
\packageDescription{pairedfda}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{pairedfda}
\packageIndices{pairedfda}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{pairedfda}

Maintainer: \packageMaintainer{pairedfda}
}
\references{
Zhou, L., Huang, J. Z., & Carroll, R. J. (2008). Joint modelling of paired sparse functional data using principal components. Biometrika, 95(3), 601-619.
}

\keyword{ package }
\seealso{
\code{\link{predata}, \link{plt.data}, \link{Kf_CV}, \link{minEM}}
}
\examples{
rawdata = gen_data(n=100,
                   varres=0.01, 
                   gama=2, 
                   type='t',
                   ka=2,
                   kb=2)
data = predata(nobs_y = rawdata$obs_times,
               nobs_z = rawdata$obs_times,
               time_y = rawdata$dataset[,2],
               time_z = rawdata$dataset[,2],
               y = rawdata$dataset[,3],
               z = rawdata$dataset[,4],
               knots = 10,
               order=3)
## without penalty
lambda = c(0,0,0,0)
pt_nopen = minEM(data, 
                 lambda, 
                 type='t', 
                 ka=2, 
                 kb=2, 
                 tol = 1e-4, 
                 maxiter = 100)
## with penalty
lambda_t = simplex(data,Kfold=5,ka=2,kb=2,type='t')
pt_pen = minEM(data, 
               lambda_t, 
               type='t', 
               ka=2, 
               kb=2, 
               tol = 1e-4, 
               maxiter = 100)
}

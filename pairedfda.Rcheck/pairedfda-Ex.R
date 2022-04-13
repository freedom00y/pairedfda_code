pkgname <- "pairedfda"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "pairedfda-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('pairedfda')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Kf_CV")
### * Kf_CV

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Kf_CV
### Title: K-fold Cross-Validation
### Aliases: Kf_CV

### ** Examples

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
lambda_nopen = c(0,0,0,0)
MSE_nopen = Kf_CV(data,lambda_nopen,K=5,ka=2,kb=2,'t')
lambda_pen = rep(0.618,4)
MSE_pen = Kf_CV(data,lambda_pen,K=5,ka=2,kb=2,'t')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Kf_CV", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gen_data")
### * gen_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gen_data
### Title: Generate simulation data
### Aliases: gen_data

### ** Examples

rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gen_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("minEM")
### * minEM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: minEM
### Title: Estimations
### Aliases: minEM

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("minEM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("orthbasis")
### * orthbasis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: orthbasis
### Title: Create bsplines
### Aliases: orthbasis

### ** Examples

rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)
time = rawdata$dataset[,2]
spl_info = orthbasis(time,knots = 8, order =3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("orthbasis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pairedfda-package")
### * pairedfda-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pairedfda-package
### Title: Paired Functional Data Analysis
### Aliases: pairedfda-package pairedfda
### Keywords: package

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pairedfda-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plt.data")
### * plt.data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plt.data
### Title: Plot the raw data
### Aliases: plt.data

### ** Examples

rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)
data = predata(nobs_y = rawdata$obs_times, 
               nobs_z = rawdata$obs_times, 
               time_y = rawdata$dataset[,2], 
               time_z = rawdata$dataset[,2], 
               y = rawdata$dataset[,3], 
               z = rawdata$dataset[,4], 
               knots = 10, 
               order=3)
plt.data(data)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plt.data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predata")
### * predata

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predata
### Title: Preprocessing data
### Aliases: predata

### ** Examples

n = 5
nobs = 1+rbinom(n,5,0.9)
time=c()
for(i in 1:n)
{
  time=c(time,0,runif(nobs[i]-1))
}
sumobs = sum(nobs)
y = rnorm(sumobs,0,1.5)
z = rnorm(sumobs,1,1)
data = predata(nobs,time,y,z,knots = 8,order=3)
rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)
data = predata(nobs_y = rawdata$obs_times, 
               nobs_z = rawdata$obs_times, 
               time_y = rawdata$dataset[,2], 
               time_z = rawdata$dataset[,2], 
               y = rawdata$dataset[,3], 
               z = rawdata$dataset[,4], 
               knots = 10, 
               order=3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predata", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simplex")
### * simplex

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simplex
### Title: Choosing tuning parameter
### Aliases: simplex

### ** Examples

rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)
data = predata(nobs_y = rawdata$obs_times, 
               nobs_z = rawdata$obs_times, 
               time_y = rawdata$dataset[,2], 
               time_z = rawdata$dataset[,2], 
               y = rawdata$dataset[,3], 
               z = rawdata$dataset[,4], 
               knots = 10, 
               order=3)
lambda_n = simplex(data,Kfold=5,ka=2,kb=2,type='n')
lambda_t = simplex(data,Kfold=5,ka=2,kb=2,type='t')
lambda_s = simplex(data,Kfold=5,ka=2,kb=2,type='s')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simplex", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

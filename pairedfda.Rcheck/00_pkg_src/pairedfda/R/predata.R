#' Preprocessing data
#' 
#' @description Convert the raw data to a data structure that can be used in this package. The new data structure includes two important entries, the design matrix and the penalty matrix.
#'
#' @param nobs_y A vector whose component is the number of Y observations of each subject
#' @param nobs_z A vector whose component is the number of Z observations of each subject
#' @param time_y A vector records Y observation time of all subjects
#' @param time_z A vector records Z observation time of all subjects
#' @param y A vector records first response variable of all subjects
#' @param z A vector records second response variable of all subjects
#' @param knots The full set of knots used to define the basis functions. It can be a number or a vector. If it is a number, the full set of knots is  seq(0,1,1/knots); if it is a vector, the full set of knots is this vector.
#' @param order Order of bsplines
#'
#' @return 
#' A list including number of subjects (n), obs times (nobs_y and nobs_z), the time we observe data (time_y and time_z),paired data (y and z), the full set of knots (knots), design matrix (B_y and B_z) and penalty matrix (Omega) 
#' 
#' @details Suppose the data include the information of n subjects. 
#' @details For each subject, there are m_i^{(y)} (i=1,...n) observations of Y which are observed at time \eqn{t_i =  (t_{i1}^{(y)},...t_{i m_i}^{(y)})}.
#' @details For each subject, there are m_i^{(z)} (i=1,...n) observations of Z which are observed at time \eqn{t_i =  (t_{i1}^{(z)},...t_{i m_i}^{(z)})}.
#' @details Suppose the data are generated from potential curves Y_i and Z_i. Then 
#' @details nobs_y = (m_1^{(y)},...,m_n^{(y)}), time_y = (t_1^{(y)},...t_n^{(y)}), y = (Y_1(t_1^{(y)}),...,Y_n(t_n^{(y)})),
#' @details nobs_z = (m_1^{(z)},...,m_n^{(z)}), time_z = (t_1^{(z)},...t_n^{(z)}), z = (Z_1(t_1^{(z)}),...,Z_n(t_n^{(z)})).
#' 
#' @seealso \code{\link{orthbasis}}
#' @export
#'
#' @examples
#' n = 5
#' nobs = 1+rbinom(n,5,0.9)
#' time=c()
#' for(i in 1:n)
#' {
#'   time=c(time,0,runif(nobs[i]-1))
#' }
#' sumobs = sum(nobs)
#' y = rnorm(sumobs,0,1.5)
#' z = rnorm(sumobs,1,1)
#' data = predata(nobs,time,y,z,knots = 8,order=3)
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

predata <- function(nobs_y,nobs_z,time_y,time_z,y,z,knots = 10, order =3)
{
  nsum_y = sum(nobs_y)
  if(length(y)!=nsum_y | length(time_y)!=nsum_y )
  {
    print("The obs times of Y and the length of data don't match!")
    break
  }
  nsum_z = sum(nobs_z)
  if(length(z)!=nsum_z | length(time_z)!=nsum_z )
  {
    print("The obs times of Z and the length of data don't match!")
    break
  }
  
  if(length(nobs_y)!=length(nobs_z) )
  {
    print("The number of subjects don't match!")
    break
  }
  
  n = length(nobs_y)
  bas_y = orthbasis(time_y,knots,order)
  bas_z = orthbasis(time_z,knots,order)
  
  data = list( n = n,
               nobs_y = nobs_y,
               nobs_z = nobs_z,
               time_y = time_y,
               time_z = time_z,
               y      = y,
               z      = z,
               B_y    = bas_y$B,
               B_z    = bas_z$B, 
               knots  = bas_y$knots,
               Omega  = bas_y$Omega)

  return(data)
}
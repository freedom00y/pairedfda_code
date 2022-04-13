#' Create bsplines
#'
#' @param x The points you are interested in.
#' @param knots The full set of knots used to define the basis functions. It can be a number or a vector. If it is a number, the full set of knots is seq(0,1,1/knots); if it is a vector, knots is this vector.
#' @param order Order of bsplines. Default is 3 which means cubic splines.
#'
#' @return A list including three elements: The full set of knots (knots), design matrix (B) and penalty matrix (Omega).  
#' 
#' @details Suppose we create n othogonal bspline basis b=(b_1,...,b_n), where b_i is a function. 
#' Then the design matrix B = b(x) = (b_1(x),...,b_n(x)). 
#' The number of rows is equal to the length of x. The number of columns is equal to the number of basis.
#' The penalty matrix Omega = integral b''(t)b''(t)^T dt.
#' 
#' @import orthogonalsplinebasis
#' @export
#'
#' @examples
#' rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)
#' time = rawdata$dataset[,2]
#' spl_info = orthbasis(time,knots = 8, order =3)
orthbasis <- function(x, knots = 10, order = 3) {
  len = length(knots)
  if(len == 1){
    innknots = seq(0,1,1/knots)
  }else{
    innknots = knots
  }
  
  knots = expand.knots(innknots,order+1) # set knots
  obase = OBasis(knots,order + 1)        # get orthogonal basis

  B     = evaluate(obase, x)
  Omega = OuterProdSecondDerivative(obase)

  result = list( knots = innknots,
                 B = B,
                 Omega = Omega)
  return(result)
}

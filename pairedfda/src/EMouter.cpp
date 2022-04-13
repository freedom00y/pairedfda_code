#include <RcppArmadillo.h>
#include "EMinner.h"
#include "loglike.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//' Estimations
//' 
//' Use this function to estimate parameters
//' 
//' @param data Processed data. Use predata to preprocess first.
//' @param lambda The penalty patameters, a vector with 4 components showing in the following order, the first mean curve \eqn{(\lambda_\mu)}, the second mean curve \eqn{(\lambda_\nu)}, the first pcs \eqn{(\lambda_f)}, the second pcs \eqn{(\lambda_g)}.
//' @param type Model type. 'n' means normal distribution, 't' means student-t distribution, 's' means slash distribution
//' @param ka Number of pcs for the first reponse variable (Y)
//' @param kb Number of pcs for the second reponse variable (Z)
//' @param tol Tolerance of the EM algorithm. The default one is 1e-4.
//' @param maxiter Maximum iteration time. The default is 100 times.
//' 
//' @return The estimated parameters
//' 
//' @details  We suppose the model is 
//' \deqn{Y_i = B_i \theta_\mu + B_i f \alpha_i + \epsilon_i, Z_i = B_i \theta_\nu + B_i g \beta_i + \xi_i,}
//' where \eqn{(\alpha_i, \beta_i)} and residuals follow normal, t or slash distribution. We denote that 
//' \deqn{(\alpha_i, \beta_i)|u_i \sim N(0,\Sigma_{\alpha\beta}/u_i), \epsilon_i|u_i\sim N(0,\sigma_\epsilon^2/u_i), \xi_i|u_i\sim N(0,\sigma_\xi^2/u_i), u_i\sim H(\gamma).}
//' If \eqn{(\alpha_i, \beta_i)} and residuals follow normal, it means \eqn{H(\gamma)=1} all the time; 
//' if \eqn{(\alpha_i, \beta_i)} and residuals follow t distribution, it means \eqn{H(\gamma)} is the distribution \eqn{\Gamma(\gamma/2,\gamma/2)};
//' if \eqn{(\alpha_i, \beta_i)} and residuals follow slash distribution, it means \eqn{H(\gamma)} is the distribution \eqn{\beta(\gamma,1)},
//' where \eqn{\Sigma_{\alpha\beta} = (D_a & C\\C^T & Db)}.
//' This function returns the estimation of \eqn{\sigma^2_\epsilon,\sigma^2_\xi,\theta_\mu,\theta_\nu,\theta_f,\theta_g,Da,Db,C}, the degree of freedom \eqn{\gamma}.
//' 
//' @examples 
//' rawdata = gen_data(n=100,
//'                    varres=0.01, 
//'                    gama=2, 
//'                    type='t',
//'                    ka=2,
//'                    kb=2)
//' data = predata(nobs_y = rawdata$obs_times,
//'                nobs_z = rawdata$obs_times,
//'                time_y = rawdata$dataset[,2],
//'                time_z = rawdata$dataset[,2],
//'                y = rawdata$dataset[,3],
//'                z = rawdata$dataset[,4],
//'                knots = 10,
//'                order=3)
//' ## without penalty
//' lambda = c(0,0,0,0)
//' pt_nopen = minEM(data, 
//'                  lambda, 
//'                  type='t', 
//'                  ka=2, 
//'                  kb=2, 
//'                  tol = 1e-4, 
//'                  maxiter = 100)
//' ## with penalty
//' lambda_t = simplex(data,Kfold=5,ka=2,kb=2,type='t')
//' pt_pen = minEM(data, 
//'                lambda_t, 
//'                type='t', 
//'                ka=2, 
//'                kb=2, 
//'                tol = 1e-4, 
//'                maxiter = 100)
//[[Rcpp::export]]
const List minEM(const List data, const arma::vec lambda, const char type, const int ka, const int kb, const double tol, int maxiter = 100)
{
  // Data
  mat By      = data["B_y"];
  int ncolB   = By.n_cols; // # of basis
  
  // Initialize Parameters
  vec theta_mu0 = ones<vec>(ncolB);
  vec theta_nu0 = ones<vec>(ncolB);
  mat X(ncolB,ncolB,fill::eye);
  mat theta_f0  = X.head_cols(ka);
  mat theta_g0  = X.head_cols(kb);
  double sig_eps0 = 1;
  double sig_xi0 = 1;
  mat Da0(ka,ka,fill::eye);
  mat Db0(kb,kb,fill::eye);
  mat C0(ka,kb,fill::zeros);
  double gama0 = 5;
  vec llv = zeros<vec>(maxiter);
  
  List para0 = List::create(Named("sig_eps")  = sig_eps0,
                            Named("sig_xi")   = sig_xi0,
                            Named("theta_mu") = theta_mu0,
                            Named("theta_nu") = theta_nu0,
                            Named("theta_f")  = theta_f0,
                            Named("theta_g")  = theta_g0,
                            Named("Da")       = Da0,
                            Named("Db")       = Db0,
                            Named("C")        = C0,
                            Named("gama")     = gama0);
  
  // EM steps
  List para1    = EM(para0,data,lambda,type);
  double value0 = loglike(data,para0,type);
  double value1 = loglike(data,para1,type);
  double differ = 1.0-value1/value0;
  llv(0) = value0;
  
  int iter = 1;
  llv(iter) = value1;
  while( iter<(maxiter-1) && fabs(differ)>tol)
  {
    para0  = para1;
    value0 = value1;
    para1  = EM(para0,data,lambda,type);
    value1 = loglike(data,para1,type);
    differ = (value0-value1)/fabs(value0);
    if(differ<-0.05)
    {
      cout<<"value0="<<value0<<endl;
      cout<<"value1="<<value1<<endl;
      cout<<"differ="<<differ<<endl;
      cout<<"Increase more than 5%. Force to break!"<<endl;
      para1 = para0;
      value1 = value0;
      break;
    }
    iter   = iter+1;
    llv(iter) = value1;
    //cout<<"value1="<<value1<<endl;
  }
  para1["-2ll"] = llv.head(iter+1);
  para1["iter"] = iter+1;
  
  if(type=='n') para1["gama"]="NA";
  
  return para1;
}

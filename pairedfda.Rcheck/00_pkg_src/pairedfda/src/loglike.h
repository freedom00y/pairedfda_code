#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//' Evaluate the loglikelihood value
//' @param data postprocessed data
//' @param para parameter set
//' @return -2*log-likelihood value
//' @keywords internal
//' @importFrom heavy pgamma.deriv
//[[Rcpp::export]]
double loglike(const List data, const List para, const char type)
{
  // Data
  int n       = data["n"];
  vec y       = data["y"];
  vec z       = data["z"];
  vec nobs_y  = data["nobs_y"];
  vec nobs_z  = data["nobs_z"];
  mat By      = data["B_y"];
  mat Bz      = data["B_z"];
  
  // Parameter
  double gama = para["gama"];
  vec    mu   = para["theta_mu"];
  vec    nu   = para["theta_nu"];
  mat    f    = para["theta_f"];
  mat    g    = para["theta_g"];
  mat    Da   = para["Da"];
  mat    Db   = para["Db"];
  mat    C    = para["C"];
  double eps  = para["sig_eps"];
  double xi   = para["sig_xi"];
  
  int ka = Da.n_rows;
  int kb = Db.n_rows;
  int kab= ka+kb;
  double total = 0;
  vec temp0 = "0";
  vec ind_y = join_vert(temp0,cumsum(nobs_y));
  vec ind_z = join_vert(temp0,cumsum(nobs_z));
  mat Bf = By*f;
  mat Bg = Bz*g;
  vec res_y = y - By*mu;
  vec res_z = z - Bz*nu;
  
  mat sig_ab = zeros(kab,kab);
  sig_ab(span(0,ka-1),span(0,ka-1))     = Da;
  sig_ab(span(0,ka-1),span(ka,kab-1))   = C;
  sig_ab(span(ka,kab-1),span(0,ka-1))   = C.t();
  sig_ab(span(ka,kab-1),span(ka,kab-1)) = Db;
  mat inv_sig_ab = inv(sig_ab); 
  
  
  for(int i=0;i<n;i++)
  {
    int ni_y  = nobs_y(i);
    int ni_z  = nobs_z(i);
    int ind1_y = ind_y(i);
    int ind2_y = ind_y(i+1)-1;
    int ind1_z = ind_z(i);
    int ind2_z = ind_z(i+1)-1;
    mat Bfi  = Bf.rows(ind1_y,ind2_y);
    mat Bgi  = Bg.rows(ind1_z,ind2_z);
    mat Bth  = zeros<mat>(ni_y+ni_z,ka+kb);
    Bth(span(0,ni_y-1),span(0,ka-1))        = Bfi;
    Bth(span(ni_y,ni_y+ni_z-1),span(ka,ka+kb-1)) = Bgi;
    vec resid  = join_cols(res_y.subvec(ind1_y,ind2_y),res_z.subvec(ind1_z,ind2_z));
    vec sigvec = join_cols(eps*ones<vec>(ni_y),xi*ones<vec>(ni_z) );
    
    mat sig_i  = zeros<mat>(ni_y+ni_z,ni_y+ni_z);
    sig_i(span(0,ni_y-1),span(0,ni_y-1))      = Bfi*Da*Bfi.t() + eps*eye(ni_y,ni_y);
    sig_i(span(0,ni_y-1),span(ni_y,ni_y+ni_z-1))   = Bfi*C*Bgi.t();
    sig_i(span(ni_y,ni_y+ni_z-1),span(ni_y,ni_y+ni_z-1))= Bgi*Db*Bgi.t() + xi*eye(ni_z,ni_z);
    sig_i(span(ni_y,ni_y+ni_z-1),span(0,ni_y-1))   = Bgi*C.t()*Bfi.t();
    
    mat inv_Ssig     = diagmat(1/sigvec);
    mat BthSs        = Bth.t()*inv_Ssig;
    mat inner        = BthSs*Bth+inv_sig_ab;
    mat inv_inner    = inv(inner);
    inv_inner.elem( find(abs(inv_inner)<1e-10) ).zeros();
    mat inv_Sigi     = inv_Ssig-BthSs.t()*inv_inner*BthSs;
    mat sig_inv      = inv_Sigi*resid;
    double delta2    = sum(resid%sig_inv);
    
    switch(type)
    {
    case 'n':
      total += log(det(sig_i))+delta2+(ni_y+ni_z)*log(2*datum::pi);
      break;
    case 't':
      total += log(det(sig_i))+(ni_y+ni_z)*log(gama) + (gama+ni_y+ni_z)*log(1+delta2/gama) -
        2*lgamma( 0.5*(gama+ni_y+ni_z) ) + 2*lgamma(0.5*gama)+(ni_y+ni_z)*log(datum::pi);
      break;
    case 's':
      double m = gama+(ni_y+ni_z)/2;
      double nn = delta2/2;
      Environment pkg = Environment::namespace_env("heavy");
      Function pg = pkg["pgamma.deriv"];
      
      Rcpp::NumericVector p = pg(1, Rcpp::Named("shape")=m, Rcpp::_["scale"]=nn, Rcpp::_["deriv"]=0);
      total += log(det(sig_i))+(ni_y+ni_z)*log(2*datum::pi)-2*log(gama)-2*( -m*log(nn) + lgamma(m) + log(p(0)) );
    }
  }
  return total;
}

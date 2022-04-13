#include <RcppArmadillo.h>
#include "orth_algo.h"
using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]

//////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Sub Functions //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Estimate gamma for t dist
//' @keywords internal 
//[[Rcpp::export]]
double f(double x, int n, double slogu, double su)
{
  return -slogu-n*log(x/2)-n+su+n*R::digamma(x/2);
}

//' @keywords internal 
//[[Rcpp::export]]
double df(double x, int n, double slogu, double su)
{
  return -n/x+n/2*R::trigamma(x/2);
}

//' @keywords internal 
//[[Rcpp::export]]
double newton_md(int n, double slogu, double su)
{
  double x=1;
  double err=1.0e-5;
  while(abs(f(x,n,slogu,su))>err)
  {
    double step=1;
    double x_temp=x-step*f(x,n,slogu,su)/df(x,n,slogu,su);
    while(x_temp<0)
    {
      step=step/2;
      x_temp=x-step*f(x,n,slogu,su)/df(x,n,slogu,su);
    }
    x=x_temp;
  }
  return x;
}

// Esitimate u and logu for slash dist
//' @importFrom heavy pgamma.deriv
//' @keywords internal 
//[[Rcpp::export]]
Rcpp::List slash_u(double gama,int ni,double delta){
  Environment pkg = Environment::namespace_env("heavy");
  Function pg = pkg["pgamma.deriv"];
  //std::cout<<"shape="<<gama+ni<<std::endl;
  Rcpp::NumericVector ne = pg(1, Rcpp::Named("shape")=gama+ni+1, Rcpp::_["scale"]=delta/2, Rcpp::_["deriv"]=0);
  Rcpp::NumericVector de = pg(1, Rcpp::Named("shape")=gama+ni, Rcpp::_["scale"]=delta/2, Rcpp::_["deriv"]=0);
  vec temp = "0,1";
  Rcpp::NumericVector para_logu = pg(1, Rcpp::Named("shape")=gama+ni, Rcpp::_["scale"]=delta/2, Rcpp::_["deriv"]=temp );
  double nne = ne(0);
  double dde = de(0);
  double hat_u = nne / dde * 2 *(gama+ni) /delta;
  double hat_logu = R::digamma(gama+ni) - log(delta/2) + para_logu(1)/para_logu(0);
  return List::create(Named("u") = hat_u,
                      Named("logu")  = hat_logu);
}

//' One EM iterate
//' @param oldpara parameter set from the last step
//' @param data postprocessed data
//' @param lambda tuning parameter
//' @return new parameter set
//' @keywords internal
//' [[Rcpp::export]]
Rcpp::List EM(const Rcpp::List oldpara, const Rcpp::List data, const arma::vec lambda, const char type){
  // Parameters
  double eps0 = oldpara["sig_eps"];
  double xi0  = oldpara["sig_xi"];
  vec thmu0   = oldpara["theta_mu"];
  vec thnu0   = oldpara["theta_nu"];
  mat thf0    = oldpara["theta_f"];
  mat thg0    = oldpara["theta_g"];
  mat Da0     = oldpara["Da"];
  mat Db0     = oldpara["Db"];
  mat C0      = oldpara["C"];
  double gama0= oldpara["gama"];
  
  // Data
  int n = data["n"];
  vec y=data["y"];
  vec z=data["z"];
  vec nobs_y = data["nobs_y"];
  vec nobs_z = data["nobs_z"];
  mat By=data["B_y"];
  mat Bz=data["B_z"];
  mat Omega=data["Omega"];
  int ncolB=By.n_cols;
  
  // E Step
  mat Bf = By*thf0;
  mat Bg = Bz*thg0;
  vec res_y = y - By*thmu0;
  vec res_z = z - Bz*thnu0;
  int ka = Da0.n_rows;
  int kb = Db0.n_rows;
  int kab= ka+kb;
  
  //Initialize sig_ab, declare (u, logu, delta2), (alpha, beta), bar_sig
  mat sig_ab = zeros(kab,kab);
  sig_ab(span(0,ka-1),span(0,ka-1))     = Da0;
  sig_ab(span(0,ka-1),span(ka,kab-1))   = C0;
  sig_ab(span(ka,kab-1),span(0,ka-1))   = C0.t();
  sig_ab(span(ka,kab-1),span(ka,kab-1)) = Db0;
  mat inv_sig_ab = inv(sig_ab); 
  
  vec  u      = zeros<vec>(n);
  vec  logu   = ones<vec>(n);
  vec  delta2 = zeros<vec>(n);
  mat  alpha  = zeros<mat>(n,ka);
  mat  beta   = zeros<mat>(n,kb);
  cube bar_sig(kab,kab,n);
  vec  temp0  = "0";
  vec  ind_y    = join_vert(temp0,cumsum(nobs_y));
  vec  ind_z    = join_vert(temp0,cumsum(nobs_z));
  
  // For each subject, compute their latent variables u,alpha,beta
  for(int i=0;i<n;i++)
  {
    int ni_y = nobs_y(i);
    int ind1_y = ind_y(i);
    int ind2_y = ind_y(i+1)-1;
    int ni_z = nobs_z(i);
    int ind1_z = ind_z(i);
    int ind2_z = ind_z(i+1)-1;
    mat Bfi  = Bf.rows(ind1_y,ind2_y);
    mat Bgi  = Bg.rows(ind1_z,ind2_z);
    mat Bth  = zeros<mat>(ni_y+ni_z,ka+kb);
    Bth(span(0,ni_y-1),span(0,ka-1))        = Bfi;
    Bth(span(ni_y,ni_y+ni_z-1),span(ka,ka+kb-1)) = Bgi;
    vec resid  = join_cols(res_y.subvec(ind1_y,ind2_y),res_z.subvec(ind1_z,ind2_z));
    vec sigvec = join_cols(eps0*ones<vec>(ni_y),xi0*ones<vec>(ni_z) );
    mat sig_i  = zeros<mat>(ni_y+ni_z,ni_y+ni_z);
    sig_i(span(0,ni_y-1),span(0,ni_y-1))           = Bfi*Da0*Bfi.t() + eps0*eye(ni_y,ni_y);
    sig_i(span(0,ni_y-1),span(ni_y,ni_y+ni_z-1))   = Bfi*C0*Bgi.t();
    sig_i(span(ni_y,ni_y+ni_z-1),span(ni_y,ni_y+ni_z-1))= Bgi*Db0*Bgi.t() + xi0*eye(ni_z,ni_z);
    sig_i(span(ni_y,ni_y+ni_z-1),span(0,ni_y-1))   = Bgi*C0.t()*Bfi.t();
    mat sig_abi = sig_ab*Bth.t();
    
    // Compute inverse of Sigma_i
    // Sigma_i = Bth*Sigma_ab*Bth.t()+Sigma_sigma, where Sigma_sigma is a diagonal matrix
    // inverse(Sigma_i) = inv(Sigma_sigma) - inv(Sigma_sigma)*Bth* 
    // ( inv(Sigma_ab + Bth.t()*inv(Sigma_sigma)*Bth) *Bth.t()*inv(Sigma_sigma))
    mat inv_Ssig     = diagmat(1/sigvec);
    mat BthSs        = Bth.t()*inv_Ssig;
    mat inner        = BthSs*Bth+inv_sig_ab;
    mat inv_inner    = inv(inner);//good
    inv_inner.elem( find(abs(inv_inner)<1e-10) ).zeros();
    bar_sig.slice(i) = inv_inner;
    mat inv_Sigi     = inv_Ssig-BthSs.t()*bar_sig.slice(i)*BthSs;
    mat sig_inv      = inv_Sigi*resid;
    mat bar_ab       = sig_abi*sig_inv;
    alpha.row(i)     = bar_ab.rows(0,ka-1).t();
    beta.row(i)      = bar_ab.rows(ka,kab-1).t();
    delta2(i)        = sum(resid%sig_inv);
   
    if(type=='s')
    {
      List temp   = slash_u(gama0,ni_y/2+ni_z/2,delta2(i));
      u(i) = temp["u"];
      logu(i)= temp["logu"];
    }else if(type=='t'){
      logu(i) = R::digamma(0.5*(gama0+ni_y+ni_z))-log(0.5*(gama0+delta2(i)));
    }
  }
  
  if(type=='t'){
    u = (gama0 + nobs_y+nobs_z)/(gama0 + delta2);
  }if(type=='n'){
    u = ones<vec>(n);
  }
  u.replace(datum::nan, 1); 
  logu.replace(datum::nan, 0); 
  mat ab = join_rows(alpha,beta);
  
  // Center a and b
  rowvec a_mean = mean(alpha,0);
  rowvec b_mean = mean(beta,0);
  alpha = alpha.each_row()-a_mean;
  beta  = beta.each_row()-b_mean;
  
  // M step
  // Compute residual variance
  double sum_eps  = 0;
  double sum_xi   = 0;
  mat    sum_ubtb_y = zeros<mat>(ncolB,ncolB);
  mat    sum_ubtb_z = zeros<mat>(ncolB,ncolB);
  vec    sum_uby  = zeros<vec>(ncolB);
  vec    sum_ubz  = zeros<vec>(ncolB);
  vec sum_btbfua  = zeros<vec>(ncolB);
  vec sum_btbgub  = zeros<vec>(ncolB);
  for(int i=0;i<n;i++)
  {
    int ind1_y = ind_y(i);
    int ind2_y = ind_y(i+1)-1;
    int ind1_z = ind_z(i);
    int ind2_z = ind_z(i+1)-1;
    mat Bi_y   = By.rows(ind1_y,ind2_y);
    mat Bi_z   = Bz.rows(ind1_z,ind2_z);
    mat btb_y  = Bi_y.t()*Bi_y;
    mat btb_z  = Bi_z.t()*Bi_z;
    mat Bfi  = Bf.rows(ind1_y,ind2_y);
    mat Bgi  = Bg.rows(ind1_z,ind2_z);
    vec yi   = y.subvec(ind1_y,ind2_y);
    vec zi   = z.subvec(ind1_z,ind2_z);
    vec res1 = yi-Bi_y*thmu0-Bfi*alpha.row(i).t();
    vec res2 = zi-Bi_z*thnu0-Bgi*beta.row(i).t();
    mat sigaa= bar_sig(0,0,i,size(ka,ka,1));
    mat sigbb= bar_sig(ka,ka,i,size(kb,kb,1));
    mat A1=Bfi*sigaa*Bfi.t();
    double val1 = sum( A1.diag() ) + u(i)*sum(res1%res1);
    mat A2=Bgi*sigbb*Bgi.t();
    double val2 = sum( A2.diag() ) + u(i)*sum(res2%res2);
    sum_eps += val1;
    sum_xi  += val2;
    sum_ubtb_y+= u(i)*btb_y;
    sum_ubtb_z+= u(i)*btb_z;
    sum_uby += u(i)*Bi_y.t()*yi;
    sum_ubz += u(i)*Bi_z.t()*zi;
    sum_btbfua += u(i)*btb_y*thf0*alpha.row(i).t();
    sum_btbgub += u(i)*btb_z*thg0*beta.row(i).t();
  }
  double eps1 = sum_eps/sum(nobs_y);
  double xi1  = sum_xi/sum(nobs_z);
  
  // Compute theta_mu, theta_nu
  vec thmu1   = solve(sum_ubtb_y/n+lambda[0]*eps1*Omega,sum_uby/n-sum_btbfua/n);
  vec thnu1   = solve(sum_ubtb_z/n+lambda[1]*xi1 *Omega,sum_ubz/n-sum_btbgub/n);
  
  // Compute theta_f, theta_g
  cube sum_f1 = zeros<cube>(ncolB,ncolB,ka);
  cube sum_f2 = zeros<cube>(ncolB,1,ka);
  cube sum_g1 = zeros<cube>(ncolB,ncolB,kb);
  cube sum_g2 = zeros<cube>(ncolB,1,kb);
  
  for(int i=0;i<n;i++)
  {
    int ind1_y = ind_y(i);
    int ind2_y = ind_y(i+1)-1;
    int ind1_z = ind_z(i);
    int ind2_z = ind_z(i+1)-1;
    
    mat Bi_y   = By.rows(ind1_y,ind2_y);
    mat Bi_z   = Bz.rows(ind1_z,ind2_z);
    mat btb_y  = Bi_y.t()*Bi_y;
    mat btb_z  = Bi_z.t()*Bi_z;
    vec yi   = y.subvec(ind1_y,ind2_y);
    vec zi   = z.subvec(ind1_z,ind2_z);
    
    mat uab = u(i)*ab.row(i).t()*ab.row(i) + bar_sig.slice(i);
    mat uaa = uab(span(0,ka-1),span(0,ka-1));
    mat ubb = uab(span(ka,kab-1),span(ka,kab-1));
    for(int j=0;j<ka;j++)
    {
      vec tfua = thf0*uaa.col(j);
      sum_f1.slice(j) += uaa(j,j)*btb_y;
      sum_f2.slice(j) += u(i)*alpha(i,j)*Bi_y.t()*(yi-Bi_y*thmu1)+uaa(j,j)*btb_y*thf0.col(j)-btb_y*tfua;
    }
    
    for(int k=0;k<kb;k++)
    {
      vec tgub = thg0*ubb.col(k);
      sum_g1.slice(k) += ubb(k,k)*btb_z;
      sum_g2.slice(k) += u(i)*beta(i,k)*Bi_z.t()*(zi-Bi_z*thnu1)+ubb(k,k)*btb_z*thg0.col(k)-btb_z*tgub;
    }
  }
    
  mat thf1(ncolB,ka);
  for(int j=0;j<ka;j++)
  {
    mat temp = sum_f1.slice(j)/n+eps1*lambda[2]*Omega;
    vec la;
    mat V;
    eig_sym(la, V, temp);
    thf1.col(j) = V*diagmat(1/la)*V.t()*sum_f2.slice(j)/n;
  }
  
  mat thg1(ncolB,kb);
  for(int k=0;k<kb;k++)
  {
    mat temp = sum_g1.slice(k)/n+xi1*lambda[3]*Omega;
    vec la;
    mat V;
    eig_sym(la, V, temp);
    thg1.col(k) = V*diagmat(1/la)*V.t()*sum_g2.slice(k)/n;
  }
  
  // Orthogonal
  mat sum_bar = sum(bar_sig,2);
  mat hat_sigma = (ab.t()*diagmat(u)*ab+sum_bar)/n;
  mat Vaa = hat_sigma(span(0,ka-1),span(0,ka-1));
  mat Vbb = hat_sigma(span(ka,kab-1),span(ka,kab-1));
  mat Vab = hat_sigma(span(0,ka-1),span(ka,kab-1));
  List temp1 = orth_algo(thf1,Vaa);
  mat qf = temp1["Q"];
  mat Da1 = temp1["D"];
  List temp2 = orth_algo(thg1,Vbb);
  mat qg = temp2["Q"];
  mat Db1 = temp2["D"];
  mat C1  = qf.t()*thf1*Vab*thg1.t()*qg;
  mat after = zeros(kab,kab);
  after(span(0,ka-1),span(0,ka-1))     = Da1;
  after(span(0,ka-1),span(ka,kab-1))   = C1;
  after(span(ka,kab-1),span(0,ka-1))   = C1.t();
  after(span(ka,kab-1),span(ka,kab-1)) = Db1;
  
  // Change sign in order to promise the second element in each column is positive
  thf1=qf*diagmat(sign(qf.row(1)));
  thg1=qg*diagmat(sign(qg.row(1)));

  
  // Compute gamma
  double gama1;
  switch(type)
  {
  case 's':
    gama1 = -n/sum(logu);
    break;
  case 't':
    gama1 = newton_md(n,sum(logu),sum(u));
    break;
  case 'n':
    gama1=gama0;
    break;
  }
  return List::create(Named("sig_eps") = eps1,
                      Named("sig_xi")  = xi1,
                      Named("theta_mu")= thmu1,
                      Named("theta_nu")= thnu1,
                      Named("theta_f") = thf1,
                      Named("theta_g") = thg1,
                      Named("Da")      = Da1,
                      Named("Db")      = Db1,
                      Named("C")       = C1,
                      Named("gama")    = gama1,
                      Named("alpha")   = alpha,
                      Named("beta")    = beta,
                      Named("u")       = u
  );
}

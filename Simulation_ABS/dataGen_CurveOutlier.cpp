#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;

// const arma::vec mu_t(const arma::vec t);
// const arma::vec nu_t(const arma::vec t);
// const arma::vec f_y1(const arma::vec t);
// const arma::vec f_y2(const arma::vec t);
// const arma::vec f_z1(const arma::vec t);
// const arma::vec f_z2(const arma::vec t);
  
//[[Rcpp::depends(RcppArmadillo)]]
//' @keywords internal 
//[[Rcpp::export]]
const arma::vec mu_t(const arma::vec t){
    int len=t.n_elem;
    return 2.5*ones(len)+2.5*t+2.5*exp(-20*pow(t-0.6,2));
}

//' @keywords internal 
//[[Rcpp::export]]
const arma::vec nu_t(const arma::vec t){
    int len=t.n_elem;
    return 7.5*ones(len)-1.5*t-2.5*exp(-20*pow(t-0.3,2));
}

//' @keywords internal 
//[[Rcpp::export]]
const arma::vec f_y1(const arma::vec t){
    int len=t.n_elem;
    //return ones(len);
    return sqrt(15)/(1+sqrt(5))*(pow(t,2)+1/sqrt(5)*ones(len));
    //return sqrt(2)*sin(2*datum::pi*t);
}

//' @keywords internal 
//[[Rcpp::export]]
const arma::vec f_y2(const arma::vec t){
    int len=t.n_elem;
    //return -sqrt(3)*ones(len)+12*t;
    return sqrt(15)/(sqrt(5)-1)*(pow(t,2)-1/sqrt(5)*ones(len));
    //return sqrt(2)*cos(2*datum::pi*t);
}

//' @keywords internal 
//[[Rcpp::export]]
const arma::vec f_z1(const arma::vec t){
    return sqrt(2)*cos(2*datum::pi*t);
}

//' @keywords internal 
//[[Rcpp::export]]
const arma::vec f_z2(const arma::vec t){
    return sqrt(2)*sin(2*datum::pi*t);
}

//' Generate simulation data
//' 
//' @description Use this function to generate data to get familiar with this package. The simulated data have the same setting as the paper describes. 
//' 
//' @param n number of subjects
//' @param varres variance or conditional variance of residuals
//' @param gama the degree of freedoms for slash and t distribution
//' @param type model type. 'n' means normal, 't' means t distribution, 's' means slash distribution
//' @param ka number of pcs for the first response variable (Y)
//' @param kb number of pcs for the second response variable (Z)
//' 
//' @return simulated data
//' 
//' @examples 
//' rawdata = gen_data(n=100,varres=0.01, gama=2, type='t',ka=2,kb=2)
//[[Rcpp::export]]
List gen_data_curveoutlier(const int& n, const double& varres, const double& percent=0.01, const int ka=2, const int kb=2){
  // Generate latent variable
  vec u = ones(n);
  // True parameters
  vec va = "1,0.25" ;
  mat Da = diagmat(va);
  vec vb = "1.44,0.36" ;
  mat Db = diagmat(vb);
  mat C  = "-0.4,0.12;0.15,-0.05";
  mat Sigma = join_rows( join_cols(Da,C), join_cols(C.t(),Db) );
  double eps0 = varres;
  double xi0  = varres;
  mat ab = mat(n, ka+kb, fill::zeros);
  vec nknot = 1 + rbinom(n,15,0.9);
  int ntotal = sum(nknot);
  mat dataset = mat(ntotal,4,fill::zeros);
  
  // Generate Y and Z
  int k1 = 0, k2;
  vec temp0 = "0";
  for(int i=0;i<n;i++)
  {
    vec mean(ka+kb, fill::zeros);
    mat var   = Sigma/u(i);
    ab.row(i) = mvnrnd(mean,var,1).t();
    if(i<floor(n*percent)) {
      rowvec sgn = 2*round(randu(4).t())-1;
      //cout << sgn << endl;
      ab.row(i) = ab.row(i) + sgn%randu(4,distr_param(0,20)).t();
      //cout << ab.row(i) << endl;
    }
    vec temp_knot = randu<vec>(nknot(i)-1);
    vec obs_time  = sort(join_vert(temp0,temp_knot));
    k2 = k1+nknot(i)-1;
    dataset(span(k1,k2),0) = i*ones(k2-k1+1,1);
    dataset(span(k1,k2),1) = obs_time;
    vec eps = mvnrnd( zeros(nknot(i)), eps0/u(i)*eye(nknot(i),nknot(i)), 1);
    vec xi  = mvnrnd( zeros(nknot(i)), xi0/u(i)*eye(nknot(i),nknot(i)), 1);
    dataset(span(k1,k2),2) = mu_t(obs_time)+ab(i,0)*f_y1(obs_time)+ab(i,1)*f_y2(obs_time)+eps;
    dataset(span(k1,k2),3) = nu_t(obs_time)+ab(i,2)*f_z1(obs_time)+ab(i,3)*f_z2(obs_time)+xi;
    k1 = k2 + 1;
  }
  return List::create(Named("n")         = n,
                      Named("dataset")   = dataset,
                      Named("obs_times") = nknot,
                      Named("latent")    = u,
                      Named("score")     = ab
                      );
}



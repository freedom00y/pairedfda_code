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


//library('Rcpp')
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <RcppArmadillo.h>

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>
#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;






//*******************************************************************//
//                MAIN FUNC                        //
//*******************************************************************//
//' MRAID
//' 
//' @export
// [[Rcpp::export]]

SEXP MRAID_CPP(SEXP betaxin, SEXP betayin, SEXP Sigma1sin, SEXP Sigma2sin, SEXP samplen1, SEXP samplen2, SEXP Gibbsnumberin, SEXP burninproportion, SEXP initial_betain, SEXP pi_beta_shape_in, SEXP pi_beta_scale_in,
SEXP pi_c_shape_in, SEXP pi_c_scale_in, SEXP pi_1_shape_in, SEXP pi_1_scale_in,SEXP pi_0_shape_in, SEXP pi_0_scale_in, SEXP maxvarin){// *
try{
	const int Gibbs_number = Rcpp::as<int>(Gibbsnumberin);
	const double burnin_p = Rcpp::as<double>(burninproportion);
	const double maxvar = Rcpp::as<double>(maxvarin);
	const int n1 = Rcpp::as<int>(samplen1);  
    const int n2 = Rcpp::as<int>(samplen2); 
	const double lamda_beta1 = Rcpp::as<double>(pi_beta_shape_in); 
	const double lamda_beta2 = Rcpp::as<double>(pi_beta_scale_in); 
	const double lamda_c1 = Rcpp::as<double>(pi_c_shape_in); 
	const double lamda_c2 = Rcpp::as<double>(pi_c_scale_in); 
	const double lamda_21 = Rcpp::as<double>(pi_1_shape_in); 
	const double lamda_22 = Rcpp::as<double>(pi_1_scale_in); 
	const double lamda_31 = Rcpp::as<double>(pi_0_shape_in); 
	const double lamda_32 = Rcpp::as<double>(pi_0_scale_in); 
	
	const arma::vec betax = as<arma::vec>(betaxin);
    const arma::vec betay = as<arma::vec>(betayin);
    const arma::mat Sigma1s = as<arma::mat>(Sigma1sin);
    const arma::mat Sigma2s = as<arma::mat>(Sigma2sin);
    const arma::vec initial_beta = as<arma::vec>(initial_betain);
	//const double sigma2y = Rcpp::as<double>(sigma2yin);
	//const double sigma2z = Rcpp::as<double>(sigma2zin);
	//const double sigma2beta = Rcpp::as<double>(sigma2betain);
	//const double alpha_input = Rcpp::as<double>(alphainputin);
	int p1 = Sigma1s.n_cols, p2 = Sigma2s.n_cols;
    if (p1 != p2){
        perror("The dimensions of x1 and x2 are not matched");
    }
    int p = p1;
			 
vec betaxh=betax, betayh=betay;
	  
  mat Sigma1 = Sigma1s* (n1-1);
  mat Sigma2 = Sigma2s* (n2-1);				 
	//double sigma0_square = 1.0;
	double sigmaz1 = 2.0/p;
	double pi = 0.1;
	//pi means pi_beta
	double pi_2 = 0.05;
	//pi_2 mean pi_c 
	double sigma2y = 0.9;
	double y_norm = n1-1;
	double z_norm = n2-1;
    double post_sigma2y_shape = 0.5*n1-1;
	double post_sigma2z_shape = 0.5*n2-1;
	double pi_gamma_21 = 0.25;
	//pi_1 
	double pi_gamma_31 = 0.005;
	//pi_0
	double sigma2z = 0.95;	
	double sigma2gamma_1 = 1.0/p;
	//double maxvar = 100;
	double rho=0.1414;
	//double sigma2gamma_2 = 1.0/p;
	
double alpha = 0;  
  
vec latent_beta = initial_beta;  

vec pleiotropy_gamma = zeros<vec>(p);
		
    mat Identity = eye<mat>(p,p);
 
	// int Gibbs_number = 1000;
	 int burnin = ceil(burnin_p*Gibbs_number);
	vec sample_alpha(Gibbs_number);
	vec sample_pi(Gibbs_number);
	vec sample_pi_2(Gibbs_number);
	vec sample_sigma2y(Gibbs_number);
	//sigma_x in the first eqution
	vec sample_sigma2z(Gibbs_number);
	//sigma_y in the second equation
	vec sample_sigmaz1(Gibbs_number);
	//sigma_beta2
	vec sample_sigma2gamma_1(Gibbs_number);
	//vec sample_sigma2gamma_2(Gibbs_number);
	vec sample_rho(Gibbs_number);
	//vec sample_sigma2gamma_2(Gibbs_number);
	vec sample_pi_gamma_21(Gibbs_number);
	//pi_1 in the model
	vec sample_pi_gamma_31(Gibbs_number);
	//pi_0 in the model
double U_beta_1,U_gamma_1;
vec latent_gamma = zeros<vec>(p);
vec latent_gamma_pleiotropy = zeros<vec>(p);
vec I = zeros<vec>(p);
double part_1,part_2,Probability_gamma,randu_number;
double part_1_pleiotropy,part_2_pleiotropy,Probability_gamma_pleiotropy,randu_number_pleiotropy;
double part_1_I,part_2_I,Probability_I_indicator,randu_number_I;
double post_sigma2y_scale;
double post_sigma2z_scale;
double post_sigmaz1_scale;
double post_sigmaz1_shape;
double post_sigma2gamma_1_scale;
//double post_sigma2gamma_2_scale;
double post_sigma2gamma_1_shape;
//double post_sigma2gamma_2_shape;
double sigmaz1_prior_shape=p/10+1;
double sigmaz1_prior_scale=0.2;
double sigma2gamma_1_prior_shape=p/5.0+1;
double sigma2gamma_1_prior_scale=0.2;
//double sigma2gamma_2_prior_shape=p/5.0+1;
//double sigma2gamma_2_prior_scale=0.2;

mat x2left = Sigma2 - (n2-1)*Identity;
mat x1left = Sigma1 - (n1-1)*Identity;	

 for(int m=0; m<Gibbs_number; ++m){
	 
double sample_var_pleiotropy_gamma_1 = 1.0/((n2-1)/sigma2z + 1/sigma2gamma_1);

 for(int k=0; k<p; ++k){
	 
double sample_var_1 = 1.0/((n1-1)/sigma2y + (n2-1)*((alpha+rho*I(k))*(alpha+rho*I(k))/sigma2z) + 1.0/sigmaz1);
 	 
U_beta_1 = ((alpha+rho*I(k))*((n2-1)*betayh(k) - sum(Sigma2.col(k)% pleiotropy_gamma)-rho*sum(x2left.col(k)%latent_beta%I)- alpha*sum(x2left.col(k)%latent_beta))/sigma2z +

((n1-1)*betaxh(k)- sum(x1left.col(k)% latent_beta))/sigma2y)* sample_var_1;	

part_1 = exp(0.5* U_beta_1*U_beta_1/sample_var_1 + 0.5*log(sample_var_1)-0.5*log(sigmaz1)+log(pi)+ latent_gamma_pleiotropy(k)*log(pi_gamma_21)+
(1-latent_gamma_pleiotropy(k))*log(1-pi_gamma_21)+ I(k)*log(pi_2)+ (1-I(k))*log(1-pi_2));

part_2 = exp(log(1-pi)+ latent_gamma_pleiotropy(k)*log(pi_gamma_31)+ (1-latent_gamma_pleiotropy(k))*log(1-pi_gamma_31));

Probability_gamma = part_1/(part_1+part_2);

randu_number = as_scalar(randu(1));

if(randu_number <= Probability_gamma){
latent_gamma(k) = 1;
latent_beta(k) = as_scalar(randn(1)*sqrt(sample_var_1)+ U_beta_1);
} else {
latent_gamma(k) = 0;
latent_beta(k) = 0;
}

U_gamma_1 =((((n2-1)*betayh(k) -rho*sum(Sigma2.col(k)% latent_beta % I)- alpha*sum(Sigma2.col(k)% latent_beta)- sum(x2left.col(k)%pleiotropy_gamma))/sigma2z)) * sample_var_pleiotropy_gamma_1;	

part_1_pleiotropy = exp(0.5* U_gamma_1*U_gamma_1/sample_var_pleiotropy_gamma_1 + 0.5*log(sample_var_pleiotropy_gamma_1)-0.5*log(sigma2gamma_1)+
latent_gamma(k)*log(pi_gamma_21)+ (1-latent_gamma(k))*log(pi_gamma_31));

part_2_pleiotropy = exp(latent_gamma(k)*log(1-pi_gamma_21)+ (1-latent_gamma(k))*log(1-pi_gamma_31));

Probability_gamma_pleiotropy = part_1_pleiotropy/(part_1_pleiotropy+part_2_pleiotropy);

randu_number_pleiotropy = as_scalar(randu(1));

if(randu_number_pleiotropy<= Probability_gamma_pleiotropy){
latent_gamma_pleiotropy(k) = 1;
pleiotropy_gamma(k) = as_scalar(randn(1)*sqrt(sample_var_pleiotropy_gamma_1)+ U_gamma_1);
} else {
latent_gamma_pleiotropy(k) = 0;
pleiotropy_gamma(k) = 0;
}		
	
part_1_I = exp(-0.5*(rho * rho * latent_beta(k)*latent_beta(k)*(n2-1)-2*rho*latent_beta(k)*((n2-1)*betayh(k)-alpha*sum(Sigma2.col(k)% latent_beta)-
rho*sum(x2left.col(k)% latent_beta % I)- sum(Sigma2.col(k)% pleiotropy_gamma)))/sigma2z + latent_gamma(k)*log(pi_2));

part_2_I = exp(latent_gamma(k)*log(1-pi_2));

Probability_I_indicator = part_1_I/(part_1_I+part_2_I);

randu_number_I = as_scalar(randu(1));

if(randu_number_I<= Probability_I_indicator){
I(k) = 1;
} else {
I(k) = 0;
}		
}

double indicator = (as_scalar((I % latent_beta).t()* Sigma2 * (I % latent_beta))/sigma2z);

if(indicator > maxvar){
double sample_rho_variance = 1.0/(as_scalar((I % latent_beta).t()* Sigma2 * (I % latent_beta))/sigma2z);

double sample_rho_mean = (as_scalar((I % latent_beta).t()* betayh*(n2-1)-(I % latent_beta).t()*Sigma2*pleiotropy_gamma - (I % latent_beta).t()*Sigma2*latent_beta*alpha)/sigma2z)*sample_rho_variance;

sample_rho(m) = as_scalar(randn(1)*sqrt(sample_rho_variance)+ sample_rho_mean);
 
rho=sample_rho(m); 
} else {
sample_rho(m) = 0;
rho=sample_rho(m); 
}


double indicator_alpha = (as_scalar(latent_beta.t()* Sigma2 * latent_beta)/sigma2z);

if (indicator_alpha>maxvar){
  double sample_alpha_variance = 1.0/(as_scalar(latent_beta.t()* Sigma2 * latent_beta)/sigma2z);
  
  double sample_alpha_mean = (as_scalar(latent_beta.t()* betayh*(n2-1)-latent_beta.t()*Sigma2*pleiotropy_gamma-latent_beta.t()*Sigma2*(latent_beta % I)*rho)/sigma2z)*sample_alpha_variance;
  
  sample_alpha(m) = as_scalar(randn(1)*sqrt(sample_alpha_variance)+ sample_alpha_mean);
  
}else{
  sample_alpha(m) = 0;
}

double shape1 = sum(latent_gamma)+lamda_beta1;

double shape2 = sum(1.0-latent_gamma)+lamda_beta2;

sample_pi(m)= r_4beta(shape1,shape2,0,1);

pi = sample_pi(m);

double shape1_pleiotropy_21 = sum(latent_gamma % latent_gamma_pleiotropy)+lamda_21;

double shape2_pleiotropy_21 = sum(latent_gamma % (1.0-latent_gamma_pleiotropy))+lamda_22;

sample_pi_gamma_21(m)= r_4beta(shape1_pleiotropy_21,shape2_pleiotropy_21,0,1);

pi_gamma_21 = sample_pi_gamma_21(m);

double shape1_pleiotropy_31 = sum((1-latent_gamma) % latent_gamma_pleiotropy)+lamda_31;

double shape2_pleiotropy_31 = sum((1-latent_gamma) % (1.0-latent_gamma_pleiotropy))+lamda_32;

sample_pi_gamma_31(m)= r_4beta(shape1_pleiotropy_31,shape2_pleiotropy_31,0,1);

pi_gamma_31 = sample_pi_gamma_31(m);

alpha=sample_alpha(m);

double shape1_pi_2 = sum(latent_gamma % I)+lamda_c1;

double shape2_pi_2 = sum(latent_gamma % (1.0-I))+lamda_c2;

sample_pi_2(m)= r_4beta(shape1_pi_2,shape2_pi_2,0,1);

pi_2 = sample_pi_2(m);

post_sigma2y_scale = 2/(y_norm + as_scalar(latent_beta.t()* Sigma1 * latent_beta) - 2 * as_scalar(latent_beta.t()* betaxh*(n1-1)));

post_sigma2y_scale=abs(post_sigma2y_scale);

sample_sigma2y(m) = 1.0/as_scalar(randg(1, distr_param(post_sigma2y_shape,post_sigma2y_scale)));

sigma2y = sample_sigma2y(m);

post_sigma2z_scale = 2/(z_norm + as_scalar(latent_beta.t()* Sigma2 * latent_beta)*alpha*alpha + as_scalar(pleiotropy_gamma.t()* Sigma2 * pleiotropy_gamma)+ 

as_scalar((latent_beta % I).t()* Sigma2 * (latent_beta % I))*rho*rho + 2*alpha*rho*as_scalar(latent_beta.t()*Sigma2*(latent_beta % I)) + 

2*rho*as_scalar(pleiotropy_gamma.t()*Sigma2*(latent_beta % I)) +

2*alpha*as_scalar(pleiotropy_gamma.t()*Sigma2*latent_beta)- 2*rho*as_scalar((latent_beta % I).t()* betayh*(n2-1))- 2*as_scalar(pleiotropy_gamma.t()* betayh*(n2-1)) - 2 * alpha * as_scalar(latent_beta.t()* betayh*(n2-1)));

post_sigma2z_scale = abs(post_sigma2z_scale);	

sample_sigma2z(m) = 1.0/as_scalar(randg(1, distr_param(post_sigma2z_shape,post_sigma2z_scale)));

sigma2z = sample_sigma2z(m);

post_sigmaz1_scale = 2/(as_scalar(sum(latent_gamma % latent_beta % latent_beta)+2*sigmaz1_prior_scale));

post_sigmaz1_shape = 0.5*sum(latent_gamma)+ sigmaz1_prior_shape;

post_sigmaz1_scale = abs(post_sigmaz1_scale);

sample_sigmaz1(m) = 1.0/as_scalar(randg(1, distr_param(post_sigmaz1_shape,post_sigmaz1_scale)));

sigmaz1 = sample_sigmaz1(m);

post_sigma2gamma_1_scale = 2/(as_scalar(sum(latent_gamma_pleiotropy % pleiotropy_gamma % pleiotropy_gamma)+2*sigma2gamma_1_prior_scale));

post_sigma2gamma_1_shape = 0.5*sum(latent_gamma_pleiotropy)+ sigma2gamma_1_prior_shape;

post_sigma2gamma_1_scale = abs(post_sigma2gamma_1_scale);

sample_sigma2gamma_1(m) = 1.0/as_scalar(randg(1, distr_param(post_sigma2gamma_1_shape,post_sigma2gamma_1_scale)));

sigma2gamma_1 = sample_sigma2gamma_1(m);

 } 
 
 
double alpha_estimate = mean(sample_alpha.subvec((burnin-1),(Gibbs_number-1)));

double alpha_sd = stddev(sample_alpha.subvec((burnin-1),(Gibbs_number-1)));
 
double rho_estimate = mean(sample_rho.subvec((burnin-1),(Gibbs_number-1))); 

double sigma2y_estimate = mean(sample_sigma2y.subvec((burnin-1),(Gibbs_number-1)));

double sigma2z_estimate = mean(sample_sigma2z.subvec((burnin-1),(Gibbs_number-1)));

double sigmaz1_estimate = mean(sample_sigmaz1.subvec((burnin-1),(Gibbs_number-1)));

double sigma2gamma_1_estimate = mean(sample_sigma2gamma_1.subvec((burnin-1),(Gibbs_number-1)));


    return List::create(Rcpp::Named("alpha") = alpha_estimate,
	                          // Rcpp::Named("sample_alpha") = sample_alpha,
							  // Rcpp::Named("sample_pi") = sample_pi,
							  // Rcpp::Named("sample_pi_gamma_31") = sample_pi_gamma_31,
							 //Rcpp::Named("sample_pi_gamma_21") = sample_pi_gamma_21,
							  //Rcpp::Named("sample_pi_2") = sample_pi_2,
							 //Rcpp::Named("sample_rho") = sample_rho,
							 Rcpp::Named("rho") = rho_estimate,
							 Rcpp::Named("sigma2x") = sigma2y_estimate,
							 Rcpp::Named("sigma2y") = sigma2z_estimate,
							 Rcpp::Named("sigmabeta") = sigmaz1_estimate,
							 Rcpp::Named("sigmaeta") = sigma2gamma_1_estimate,
                               Rcpp::Named("sd") = alpha_sd);
                               //Rcpp::Named("sigmaX") = sigma2y,
                                //Rcpp::Named("sigmaY") = sigma2z,
                               //Rcpp::Named("sigmabeta") = sigma2beta,
                               //Rcpp::Named("loglik_seq") = loglik_out,
                               //Rcpp::Named("loglik") = loglik_max,
                               //Rcpp::Named("iteration") = Iteration-1);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}// end func

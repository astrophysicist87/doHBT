#ifndef GET_CORR_FCN_FROM_TOYS_H
#define GET_CORR_FCN_FROM_TOYS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#include <time.h>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <complex>
#include <limits>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

using namespace std;

const double PI = 3.14159265358979323846264338327950;
const double hbarC=0.197327053;  //GeV*fm
const int qnpts = 201;
const double init_q = 0.003;
const double delta_q = 0.006;
const int order = 43;
const int half_order = (order - 1)/2;
const double limit = 5.;
const double M1_4PI = 1./(4.*PI);
const double m_pion = 0.139;  //in GeV

double eps_2_bar_cmd;	//defined from 0..0.5
double eps_3_bar_cmd;	//defined from 0..0.5
double v_2_bar_cmd;	//defined from 0..0.5
double v_3_bar_cmd;	//defined from 0..0.5
double psi_2_bar_cmd;	//defined as fraction of PI
double psi_3_bar_cmd;	//defined as fraction of PI
double tau_f_cmd;
double Delta_tau_cmd;
double K_perp_cmd;
double eta_0_cmd;
double Delta_eta_cmd;
double eta_f_cmd;
double Rad_cmd;
double T0_cmd;
double flow_angle_cmd;
double Phi_K_cmd;
const double eps_2_bar_default = 0.25;
const double eps_3_bar_default = 0.;
const double v_2_bar_default = 0.;
const double v_3_bar_default = 0.;
const double psi_2_bar_default = 0.;
const double psi_3_bar_default = 0.;
const double tau_f_default = 10.;  //in fm
const double Delta_tau_default = 1.;
const double K_perp_default = 0.5;  //in GeV
const double eta_0_default = 0.;
const double Delta_eta_default = 2.3;
const double eta_f_default = 0.6;
const double Rad_default = 4.5;
const double flow_angle_default = 0.;
const double T0_default = 0.12;
const double Phi_K_default = 0.;
double r_lower=0., r_upper=4.*limit, phi_lower=-PI, phi_upper=PI;
double eta_lower=-limit, eta_upper=limit, tau_lower=0., tau_upper=4.*limit;

double my_q_out = 0., my_q_side = 0., my_q_long = 0.;

string path = "results";
//string path = ".";	//current directory
string particle_name = "Pion(+)";

typedef struct
{
   double qo, qs, ql;
   double Corr_fn, error;
}corr_fcn_data;

struct data
{
	size_t n;
	double * y;
	double * sigma;
};

struct phys_params
{
	double M_perp, T_0, eta_0, Y_rapidity;
	double Rad, eta_f, tau_f, K_perp;
	double Delta_tau, Delta_eta, beta_perp, beta_long;

	double v_3_bar, eps_3_bar, psi_3_bar, Phi_K, v_2_bar, eps_2_bar, psi_2_bar;
};

struct position
{
	double r, phi, eta, tau;
};

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

vector<int> convert_base(int number, int base)
{
	const int length = 1 + int(log(number)/log(base));
	vector<int> result(length);
	for (int diml = 1; diml <= length; diml++) result[diml-1] = (number%int(pow(base,diml)))/int(pow(base,diml-1));
	//for (int dim = 1; dim <= length; dim++) {result[length-dim] = (number%int(pow(base,dim)))/int(pow(base,dim-1));  cout << result[length-dim] << " ";}

	return (result);
}

vector<int> convert_base(int number, int base, int length)
{
	//const int length = 1 + int(log(number)/log(base));
	vector<int> result(length);
	for (int diml = 1; diml <= length; diml++) result[diml-1] = (number%int(pow(base,diml)))/int(pow(base,diml-1));
	//for (int dim = 1; dim <= length; dim++) {result[length-dim] = (number%int(pow(base,dim)))/int(pow(base,dim-1));  cout << result[length-dim] << " ";}

	return (result);
}

bool is_multiple(int x, int y)
{
	//x is number which might be a multiple of y
	return (double(x/y) == double(x)/double(y));
}

double heaviside_theta(double x)
{
	return ((x >= 0.) ? 1. : 0.);
}

struct phys_params inits;

vector<double>* xi_ptr;
vector<double>* wi_ptr;
vector<double>* xi_0pinf_ptr;
vector<double>* wi_0pinf_ptr;
vector<double>* xi_minfpinf_ptr;
vector<double>* wi_minfpinf_ptr;
vector<double>* xiq_ptr;
vector<double>* wiq_ptr;
vector<double>* xiq_0pinf_ptr;
vector<double>* wiq_0pinf_ptr;
vector<double>* xiq_minfpinf_ptr;
vector<double>* wiq_minfpinf_ptr;

vector<position>* grid_ptr;

vector<double>* lower_limits_vec;
vector<double>* upper_limits_vec;
vector<int>* integration_intervals_vec;

      //store correlation functions
      double* Correl_1D_out;
      double* Correl_1D_out_err;
      double* Correl_1D_side;
      double* Correl_1D_side_err;
      double* Correl_1D_long;
      double* Correl_1D_long_err;
      double*** Correl_3D;
      double*** Correl_3D_err;
      double* q_out;
      double* q_side;
      double* q_long;

#endif

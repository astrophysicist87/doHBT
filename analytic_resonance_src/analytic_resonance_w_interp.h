#ifndef ANALYTIC_RESONANCE_H
#define ANALYTIC_RESONANCE_H

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

#define ASSUME_ETA_SYMMETRIC 		0		// 1 means integrate only over eta_s = 0..eta_s_max, and multiply by 2 or 0 to speed up calculations
							// 0 means just integrate over given range of eta_s without worrying about symmetry

#include "../src/Arsenal.h"

typedef struct
{
   int* nbody;
   int* resonance_sign;
   double* resonance_mass;
   double* resonance_mu;
   double* resonance_gspin;
   double* resonance_Gamma;
   double* resonance_total_br;
   double** resonance_decay_masses;
}resonance_info;

const double PI = 3.14159265358979323846264338327950;
const double twopi = 2.*PI;
const double hbarC=0.197327053;  //GeV*fm
const double MeVToGeV = 0.001;
const double hbarC3=hbarC*hbarC*hbarC;
const int order = 43;
const int order2 = order;
int n_resonance;

const int nrpts = order;
const int nphipts = order;
const int netapts = order;
const int ntaupts = order;

const int nspts = order2;
const int nvpts = order2;
const int nzetapts = order2;
const int nptaupts = order2;

const double r_min = 0.0, r_max = 20.0;
const double phi_min = 0.0, phi_max = 2.0*M_PI;
const double eta_min = -5.0, eta_max = 5.0;
const double tau_min = 0.0, tau_max = 20.0;

const double M1_4PI = 1./(4.*PI);
const double m_pion = 0.13957;  //in GeV

//physical constants used
//physical constants used
const double Rad = 5.0;
const double eta0 = 0.0;
const double tau0 = 5.0;
const double Deleta = 1.2;
const double Deltau = 1.0;
const double eta_f = 0.0;
const double Jr = 0.0;
const double Tdec = 0.15;

double tau_integrated_S_prefactor, tau_integrated_Stau_prefactor, tau_integrated_Stau2_prefactor;

//physical variables
double current_KT, current_MT, current_Kphi, current_cos_Kphi, current_sin_Kphi;

const int n_K_T = 15, n_K_phi = 2;
const double K_T_min = 0.0, K_T_max = 0.7;
const double K_phi_min = 0.0, K_phi_max = 2.*M_PI;
double * K_T, * K_phi, * K_phi_weights;

string path = "results";
//string path = ".";	//current directory
string particle_name = "Pion_p";

	//vectors for gaussian integration points
	double * r_pts, * r_wts;
	double * phi_pts, * phi_wts;
	double * eta_pts, * eta_wts;
	double * tau_pts, * tau_wts;
	
	//points and weights for resonance integrals
	double * s_pts, * s_wts;
	double * v_pts, * v_wts;
	double * zeta_pts, * zeta_wts;
	double * ptau_pts, * ptau_wts;	//for proper time integration
	resonance_info resonances;

	//store spectra
	double **** resonance_spectra;
	double * PTpts, * PPhipts, * PYpts;
	const int nPTpts = 100;
	const int nPPhipts = 100;
	const int nPYpts = 100;
	const double PTmin = 0.0, PTmax = 5.0;
	const double PPhimin = 0.0, PPhimax = 2.*M_PI;
	const double PYmin = -3.0, PYmax = 3.0;

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

double place_in_range(double phi, double min, double max)
{
	while (phi < min || phi > max)
	{
		if (phi < min) phi += twopi;
		else phi -= twopi;
	}

	return (phi);
}

double g(double s);
double get_Q();
void Set_resonance_integration_points(double smin, double smax);
double Direct_contributions_to_pion_spectra(double pT, double y, double pphi);
double Resonance_decay_contributions_to_pion_spectra(double pT, double y, double pphi, int reso_idx);
double S_thermal(double r, double phi, double eta, double tau, double PT, double Y, double Phi, int reso_idx = 0);
double tauintegrated_S_thermal(double r, double phi, double eta, double PT, double Y, double Phi, int reso_idx = 0);
void Create_integrations_points_and_weights();
void Compute_direct_resonance_spectra();
void Compute_direct_pion_spectra_OLD();
int read_in_resonances(resonance_info * resonances);


#endif

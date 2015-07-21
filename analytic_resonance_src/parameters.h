#ifndef PARAMETERS_H
#define PARAMETERS_H

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

#define ASSUME_ETA_SYMMETRIC 		1		// 1 means integrate only over eta_s = 0..eta_s_max, and multiply by 2 or 0 to speed up calculations
							// 0 means just integrate over given range of eta_s without worrying about symmetry

//#include "resonance_info.h"
//#include "../src/Arsenal.h"

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
const int order = 11;
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
const double eta_min = 0.0, eta_max = 5.0;
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

const int n_K_T = 15, n_K_phi = 3;
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
	const int nPTpts = 25;
	const int nPPhipts = 25;
	const int nPYpts = 25;
	const double PTmin = 0.0, PTmax = 5.0;
	const double PPhimin = 0.0, PPhimax = 2.*M_PI;
	const double PYmin = 0.0, PYmax = 5.0;

#endif

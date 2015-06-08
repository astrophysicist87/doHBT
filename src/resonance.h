#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "gauss_quadrature.h"
//#include "gauss.h"

using namespace std;

//kinematic info
double pstar, Estar, Yp, Ym, DeltaY;
double MTbar, DeltaMT, MTp, MTm;

//pair momentum info, currently assumes pT != 0
double py, pT, pphi, mT, m;

//resonance momentum info
double PY, PT, PPhi, MT, M, PPhip, PPhim;
double m2, m3, Gamma, br, Jr;

double* Pp;
double* Pm;

const int n_tau_pts = 10;
const int n_zeta_pts = 10;
const int n_v_pts = 10;
const int n_s_pts = 10;


const double v_min = -1.;
const double v_max = 1.;
const double zeta_min = 0.;
const double zeta_max = M_PI;
const double vmin = -1.;
const double vmin = -1.;


double* tau_pts;
double* zeta_pts;
double* v_pts;
double* s_pts;
double* tau_wts;
double* zeta_wts;
double* v_wts;
double* s_wts;

tau_pts = new double [n_tau_pts];
zeta_pts = new double [n_zeta_pts];
v_pts = new double [n_v_pts];
s_pts = new double [n_s_pts];
tau_wts = new double [n_tau_pts];
zeta_wts = new double [n_zeta_pts];
v_wts = new double [n_v_pts];
s_wts = new double [n_s_pts];

Pp = new double [4];
Pm = new double [4];

//all expressions taken directly from Wiedemann & Heinz (1997)

void set_gaussian_pts_and_wts();
void set_pstar(double s);
void set_Estar();
void set_DeltaY();
void set_Ypm();
void set_MTbar();
void set_DeltaMT();
void set_MTpm();
void set_MT(double zeta);
void set_PY(double v);
void set_PPhi_vars();
void set_Ppm();
double get_Q();
double g(double s);
double s_integ();
double v_integ();
double zeta_integ();
double tau_integ();
double C(double PK[], double tau);

//End of file

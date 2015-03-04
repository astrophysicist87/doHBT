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

#include "doHBT.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "gauss.h"

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

void set_gaussian_pts_and_wts()
{
	//syntax: int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[])
	gauss_quadrature(n_tau_pts, 5, 0.0, 0.0, 0.0, Gamma, tau_pts, tau_wts);
	gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0, zeta_min, zeta_max, zeta_pts, zeta_wts);
	gauss_quadrature(n_v_pts, 1, 0.0, 0.0, v_min, v_max, v_pts, v_wts);
	gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min, s_max, s_pts, s_wts);

	return;
}

void set_pstar(double s)
{
	pstar = sqrt(((M+m)*(M+m) - s)*((M-m)*(M-m) - s)/(2.0*M));
	return;
}

void set_Estar()
{
	Estar = sqrt(m*m + pstar*pstar);
	return;
}

void set_DeltaY()
{
	double psBmT = pstar / mT;
	DeltaY = log(psBmT + sqrt(1.+psBmT*psBmT));
	return;
}

void set_Ypm()
{
	Yp = py + DeltaY(pstar);
	Ym = py - DeltaY(pstar);
	return;
}

void set_MTbar()
{
	double mT_ch_PY_py = mT*cosh(PY-py);
	double x2 = mT_ch_PY_py*mT_ch_PY_py - pT*pT;
	double num = Estar*M*mT_ch_PY_py;
	MTbar = num/x2;
	return;
}

void set_DeltaMT()
{
	double mT_ch_PY_py = mT*cosh(PY-py);
	double x2 = mT_ch_PY_py*mT_ch_PY_py - pT*pT;
	double disc = Estar*Estar - x2;
	DeltaMT = M*pT*sqrt(disc)/x2;
	return;
}

void set_MTpm()
{
	MTp = MTbar + DeltaMT;
	MTm = MTbar - DeltaMT;
	return;
}

void set_MT(double zeta)
{
	MT = MTbar + zeta*DeltaMT;
	return;
}

void set_PY(double v)
{
	PY = py + v*DeltaY;
	return;
}

void set_PPhi_vars()
{	//assumes currently that pT != 0
	double cos_PPhi_tilde = (mT*MT*cosh(PY-py) - Estar*M)/(pT*PT);
	double PPhi_tilde = arccos(cos_PPhi_tilde);
	PPhip = PPhi_tilde;
	PPhim = -PPhi_tilde;
	return;
}

void set_Ppm()
{
	Pp[0] = MT*cosh(PY);
	Pp[1] = PT*cos(PPhip);
	Pp[2] = MT*sin(PPhip);
	Pp[3] = MT*sinh(PY);
	Pm[0] = MT*cosh(PY);
	Pm[1] = MT*cos(PPhim);
	Pm[2] = MT*sin(PPhim);
	Pm[3] = MT*sinh(PY);
	return;
}

double doHBT::get_Q()
{
	double smin = (m2+m3)*(m2+m3);
	double smax = (M-m)*(M-m);
	double sum = 0.;
	
	for (int is = 0; is < n_s_pts; is++)
	{
		double sp = s_pts[is];
		double f1 = (M+m)*(M+m) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		sum += s_wts[is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double doHBT::g(double s)
{
	//assume n_body == 3...
	double pre_f = (M*br)/(2.*M_PI*s);
	double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
	double den = get_Q();
	double g_res = pre_f * num / den;

	return g_res;
}


double doHBT::s_integ()
{
	double s_sum = 0.;
	
	if (n_body == 2)
	{
		//then g(s) is delta-function, skip s-integration entirely
		//double s_loc = m2*m2;
		set_pstar(m2*m2);
		set_Estar();
		set_DeltaY();
		set_Ypm();
		s_sum = br*v_integ()/(4.*M_PI*pstar);
	}
	else if (n_body == 3)
	{
		for (is = 0; is < n_s_pts; is++)
		{
			double s_loc = s_pts[is];
			set_pstar(s_loc);
			set_Estar();
			set_DeltaY();
			set_Ypm();
	
			//g(s) only defined here for n_body == 3
			double s_factor = s_wts[is]*g(s_loc);
			s_sum += s_factor*v_integ();
		}
	}

	return s_sum;
}

double doHBT::v_integ()
{
	double v_sum = 0.;
	
	for (iv = 0; iv < n_v_pts; iv++)
	{
		double v_loc = v_pts[iv];
		set_PY();
		set_MTbar();
		set_DeltaMT();
		set_MTpm()

		double v_factor = v_wts[iv]*DeltaY/(mT*mT*cosh(v_loc*DeltaY)*cosh(v_loc*DeltaY) - pT*pT);
		v_sum += v_factor*zeta_integ();
	}
	
	return v_sum;
}

double doHBT::zeta_integ()
{
	double zeta_sum = 0.;
	
	for (izeta = 0; izeta < n_zeta_pts; izeta++)
	{
		double zeta_loc = zeta_pts[izeta];

		double zeta_factor = zeta_wts[izeta]*(MTbar + DeltaMT*cos(zeta_loc));
		zeta_sum += zeta_factor*tau_integ();
	}
}

double doHBT::tau_integ()
{
	double tau_sum = 0.;
	
	for (int itau = 0; itau < n_tau_pts; itau++)
	{
		double tau_loc = tau_pts[itau];
		double tau_factor = Gamma*exp(-Gamma*tau_loc)*tau_wts[itau];
		tau_sum += tau_factor*(C(Pp, tau_loc) + C(Pm, tau_loc));
	}
}

double doHBT::C(double PK[], double tau)
{
	//...
}

//End of file

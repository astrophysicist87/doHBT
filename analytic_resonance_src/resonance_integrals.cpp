#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <string>
#include <time.h>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <limits>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_matrix.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_linalg.h>

using namespace std;

#include "analytic_resonance.h"
#include "gauss.h"
//#include "CMDreader.h"

bool output_to_screen = true;
bool output_to_file = true;

string input_filename;
string input_filepath;

double m, Mres, Gamma, m2, m3, br, sign, Qfunc;
int n_body;

void Set_resonance_integration_points(double smin, double smax)
{
	s_pts = new double [nspts];
	s_wts = new double [nspts];
	v_pts = new double [nvpts];
	v_wts = new double [nvpts];
	zeta_pts = new double [nzetapts];
	zeta_wts = new double [nzetapts];
	ptau_pts = new double [nptaupts];
	ptau_wts = new double [nptaupts];
	gauss (nspts, 0, smin, smax, s_pts, s_wts);
	gauss (nvpts, 0, -1.0, 1.0, v_pts, v_wts);
	gauss (nzetapts, 0, 0.0, M_PI, zeta_pts, zeta_wts);
	gauss (nptaupts, 0, 0.0, 20.0, ptau_pts, ptau_wts);
	
	return;
}

double Resonance_decay_contributions_to_pion_spectra(double pT, double y, double pphi, int reso_idx)
{
	double PkT, PkPhi, PkY;
	double ssum = 0.0;
	m = m_pion;
	double mT = sqrt(m*m+pT*pT);
	Mres = resonances.resonance_mass[reso_idx-1];
	Gamma = resonances.resonance_Gamma[reso_idx-1];
	m2 = resonances.resonance_decay_masses[reso_idx-1][0];
	m3 = resonances.resonance_decay_masses[reso_idx-1][1];
	br = resonances.resonance_total_br[reso_idx-1];
	sign = (double)resonances.resonance_sign[reso_idx-1];
	double s_min = (m2 + m3)*(m2 + m3);
	double s_max = (Mres - m)*(Mres - m);
	n_body = 2;
	if (m2 > 1.e-6 && m3 > 1.e-6) n_body = 3;
	
	Set_resonance_integration_points(s_min, s_max);
	Qfunc = get_Q();
	
	if (n_body == 2)
	{
		//then g(s) is delta-function, skip s-integration entirely
		double spt = m2*m2;
		double pstar = sqrt( ((Mres+m)*(Mres+m) - spt)*((Mres-m)*(Mres-m) - spt)/(2.0*Mres) );
		double gs = br/(4.*M_PI*pstar);
		double s_factor = gs;
		double Estar = sqrt(m*m + pstar*pstar);
		double psBmT = pstar / mT;
		double DeltaY = log(psBmT + sqrt(1.+psBmT*psBmT));
		double Yp = y + DeltaY;
		double Ym = y - DeltaY;
		double vsum = 0.0;
		for (int iv = 0; iv < nvpts; iv++)
		{
			double vpt = v_pts[iv];
			double Y = y + vpt*DeltaY;
			PkY = Y;
			double mT_ch_Y_y = mT*cosh(vpt*DeltaY);
			double x2 = mT_ch_Y_y*mT_ch_Y_y - pT*pT;
			double v_factor = v_wts[iv]*DeltaY/sqrt(x2);
			double MTbar = Estar*Mres*mT_ch_Y_y/x2;
			double DeltaMT = Mres*pT*sqrt(Estar*Estar - x2)/x2;
			double MTp = MTbar + DeltaMT;
			double MTm = MTbar - DeltaMT;
			double zetasum = 0.0;
			//time (&rawtime);
			//timeinfo = localtime (&rawtime);
			//cerr << "Starting v-loop #" << iv << " at " << asctime(timeinfo);
			for (int izeta = 0; izeta < nzetapts; izeta++)
			{
				double zetapt = zeta_pts[izeta];
				double MT = MTbar + cos(zetapt)*DeltaMT;
				double zeta_factor = zeta_wts[izeta]*MT;
				double PT = sqrt(MT*MT - Mres*Mres);
				PkT = PT;
				double temp_cos_PPhi_tilde = (mT*MT*cosh(Y-y) - Estar*Mres)/(pT*PT);
				double temp_sin_PPhi_tilde = sqrt(1. - temp_cos_PPhi_tilde*temp_cos_PPhi_tilde);
				double PPhi_tilde = place_in_range( atan2(temp_sin_PPhi_tilde, temp_cos_PPhi_tilde), PPhimin, PPhimax);
				double PPhi_tilde_shift = place_in_range( pphi + PPhi_tilde, PPhimin, PPhimax);
				double PPhi_tilde_shiftFLIP = place_in_range( pphi - PPhi_tilde, PPhimin, PPhimax);
				double * Pp = new double [4];
				double * Pm = new double [4];
				//probably not the most elegant set-up, but does the job for now...
				Pp[0] = MT * cosh(Y);
				Pp[1] = PT * cos(PPhi_tilde_shift);
				Pp[2] = PT * sin(PPhi_tilde_shift);
				Pp[3] = MT * sinh(Y);
				Pm[0] = Pp[0];
				Pm[1] = PT * cos(PPhi_tilde_shiftFLIP);
				Pm[2] = PT * sin(PPhi_tilde_shiftFLIP);
				Pm[3] = Pp[3];
				double ptausum = 0.0;
				for (int iptau = 0; iptau < nptaupts; iptau++)
				{
					double Csum = 0.0;
					double ptau_factor = Gamma * exp(-Gamma*ptau_pts[iptau]) * ptau_wts[iptau];
					for (int tempidx = 1; tempidx <= 2; tempidx++)
					{
						if (tempidx != 1)
							PkPhi = PPhi_tilde_shiftFLIP;		//also takes Pp --> Pm
						//Csum += Interpolate_spectra();
						Csum += interpolate3D(PTpts, PPhipts, PYpts, resonance_spectra[reso_idx - 1],
									PkT, PkPhi, PkY, nPTpts, nPPhipts, nPYpts, 0, true, true);
					}
					ptausum += ptau_factor*Csum;
				}
				zetasum += zeta_factor*ptausum;
			}
			vsum += v_factor*zetasum;
		}
		ssum += Mres*s_factor*vsum;
	}
	else if (n_body == 3)
	{
		;
	}
	
	return (ssum);
}

double get_Q()
{
	double smin = (m2+m3)*(m2+m3);
	double smax = (Mres-m)*(Mres-m);
	double sum = 0.;
	
	for (int is = 0; is < nspts; is++)
	{
		double sp = s_pts[is];
		double f1 = (Mres+m)*(Mres+m) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		sum += s_wts[is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double g(double s)
{
	double g_res = br/(4.*M_PI);
	if (n_body == 3)
	{
		double pre_f = (Mres*br)/(2.*M_PI*s);
		double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
		double den = Qfunc;
		double g_res = pre_f * num / den;
	}

	return g_res;
}

//End of file

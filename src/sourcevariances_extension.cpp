#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "sourcevariances.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

#define VERBOSE 2			// specifies level of output - 0 is lowest (no output)
#define ASSUME_ETA_SYMMETRIC 1		// 1 means integrate only over eta_s = 0..eta_s_max, and multiply by 2 or 0 to speed up calculations
					// 0 means just integrate over given range of eta_s without worrying about symmetry

using namespace std;

void SourceVariances::Cal_dN_dypTdpTdphi_with_weights_polar_NEW(FO_surf* FOsurf_ptr, int reso_idx)
{
	//time_t rawtime;
  	//struct tm * timeinfo;
	//CURRENTLY USING WRONG DEFINITIONS OF SIGN AND DEGEN
	//NEED TO FIX BEFORE VERSION IS STABLE
	double sign = particle_sign;
	double degen = particle_gspin;
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	//assume reso_idx >= 1
	//double mu = 0.0;
	//these are constants along the FO surface,
	//so don't waste time updating them for each cell
	Tdec = (&FOsurf_ptr[0])->Tdec;
	double one_by_Tdec = 1./Tdec;
	Pdec = (&FOsurf_ptr[0])->Pdec;
	Edec = (&FOsurf_ptr[0])->Edec;
	//cerr << "DEBUG: " << Tdec << "  " << Pdec << "  " << Edec << endl;
	double deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
	
	//declare variables used below here to see if this speeds up code:
	double mu, tau, vx, vy, da0, da1, da2;
	double pi00, pi01, pi02, pi11, pi12, pi22, pi33;
	double temp_r, temp_phi, sin_temp_phi, cos_temp_phi, gammaT;
	
	double pT, sin_pphi, cos_pphi, sin_phi_m_pphi, cos_phi_m_pphi;
	double px, py, p0, pz, f0, deltaf, S_p, S_p_withweight;
	double symmetry_factor = 1.0;	//symmetry_factor accounts for the assumed reflection symmetry along eta direction
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
	mu = (&FOsurf_ptr[0])->particle_mu[particle_id];
	double temp;
	//int idx = 0;
		     
	for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
	for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
	for(int ieta=0; ieta < eta_s_npts; ieta++)
	{
		temp = 0.;
		//set momentum crap
		pT = SPinterp2_pT[ipt];
		sin_pphi = sin_SPinterp2_pphi[iphi];
		cos_pphi = cos_SPinterp2_pphi[iphi];
		px = pT*cos_pphi;
		py = pT*sin_pphi;
		p0 = SPinterp2_p0[ipt][ieta];
		pz = SPinterp2_pz[ipt][ieta];
		
	//time (&rawtime);
	//timeinfo = localtime (&rawtime);
	//if (VERBOSE > 0 && idx==0) cerr << "***Checkpoint #1 at " << asctime(timeinfo);
		temp = loop_over_FO_surface(FOsurf_ptr, p0, px, py, pz, mu);
cerr << p0 << "    " << px << "    " << py << "    " << pz << "    " << temp << endl;
	//time (&rawtime);
	//timeinfo = localtime (&rawtime);
	//if (VERBOSE > 0 && idx==999) cerr << "***Checkpoint #2 at " << asctime(timeinfo);
	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi] = temp*eta_s_weight[ieta]*symmetry_factor;
		//idx++;
	}

   return;
}

double SourceVariances::loop_over_FO_surface(FO_surf* FOsurf_ptr, double p0, double px, double py, double pz, double mu)
{
	double result = 0.0;
	double sign = particle_sign;
	double degen = particle_gspin;
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	Tdec = (&FOsurf_ptr[0])->Tdec;
	double one_by_Tdec = 1./Tdec;
	Pdec = (&FOsurf_ptr[0])->Pdec;
	Edec = (&FOsurf_ptr[0])->Edec;
	double deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
	double tau, vx, vy, da0, da1, da2;
	double pi00, pi01, pi02, pi11, pi12, pi22, pi33;
	double temp_r, temp_phi, sin_temp_phi, cos_temp_phi, gammaT;
	double pT, sin_pphi, cos_pphi, sin_phi_m_pphi, cos_phi_m_pphi;
	double f0, deltaf, S_p, S_p_withweight;
	double symmetry_factor = 1.0;	//symmetry_factor accounts for the assumed reflection symmetry along eta direction
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
	FO_surf* surf;
	for(int isurf=0; isurf<FO_length; isurf++)
	{
		surf = &FOsurf_ptr[isurf];
		//mu = surf->particle_mu[particle_id];
		tau = surf->tau;
		vx = surf->vx;
		vy = surf->vy;
		da0 = surf->da0;
		da1 = surf->da1;
		da2 = surf->da2;
		pi00 = surf->pi00;
		pi01 = surf->pi01;
		pi02 = surf->pi02;
		pi11 = surf->pi11;
		pi12 = surf->pi12;
		pi22 = surf->pi22;
		pi33 = surf->pi33;
		temp_r = surf->r;
		sin_temp_phi = surf->sin_phi;
		cos_temp_phi = surf->cos_phi;
		gammaT = surf->gammaT;

		//set space-time crap
		sin_phi_m_pphi = sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi;
		cos_phi_m_pphi = cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi;
		//temporarily specialize to just calculating spectra with sign = +1
		//generalize later with innermost loop over weight_function_index and spin-sign of particle
		//zvec[0] = tau*ch_eta_s[ieta];
		//zvec[1] = temp_r * cos_phi_m_pphi;
		//zvec[2] = temp_r * sin_phi_m_pphi;
		//zvec[3] = tau*sh_eta_s[ieta];

		//now get distribution function, emission function, etc.
		f0 = 1./(exp((gammaT*(p0*1. - px*vx - py*vy) - mu)*one_by_Tdec)+sign);	//thermal equilibrium distributions
		
		//viscous corrections
		deltaf = 0.;
		if (use_delta_f)
			deltaf = (1. - sign*f0)*(p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33)*deltaf_prefactor;

		//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
		S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);
		//ignore points where delta f is large or emission function goes negative from pdsigma
		if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol)) S_p = 0.0;
		S_p_withweight = S_p*tau;
		//temp += S_p_withweight*weight_function(zvec, wfi);
		//temporarily specialize to just calculating spectra with sign = +1
		//generalize later with innermost loop over weight_function_index and spin-sign of particle
		result += S_p_withweight;
	}
	return (result);
}

/*void SourceVariances::Cal_dN_dypTdpTdphi_with_weights_polar(FO_surf* FOsurf_ptr, int reso_idx)
{
	//CURRENTLY USING WRONG DEFINITIONS OF SIGN AND DEGEN
	//NEED TO FIX BEFORE VERSION IS STABLE
	double sign = particle_sign;
	double degen = particle_gspin;
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	//assume reso_idx >= 1
	//double mu = 0.0;
	//these are constants along the FO surface,
	//so don't waste time updating them for each cell
	Tdec = (&FOsurf_ptr[0])->Tdec;
	double one_by_Tdec = 1./Tdec;
	Pdec = (&FOsurf_ptr[0])->Pdec;
	Edec = (&FOsurf_ptr[0])->Edec;
	//cerr << "DEBUG: " << Tdec << "  " << Pdec << "  " << Edec << endl;
	double deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
	
	//declare variables used below here to see if this speeds up code:
	double mu, tau, vx, vy, da0, da1, da2;
	double pi00, pi01, pi02, pi11, pi12, pi22, pi33;
	double temp_r, temp_phi, sin_temp_phi, cos_temp_phi, gammaT;
	
	double pT, sin_pphi, cos_pphi, sin_phi_m_pphi, cos_phi_m_pphi;
	double px, py, p0, pz, f0, deltaf, S_p, S_p_withweight;
	double symmetry_factor = 1.0;
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
	FO_surf* surf;
	
	for(int isurf=0; isurf<FO_length; isurf++)
	{
		surf = &FOsurf_ptr[isurf];
		mu = surf->particle_mu[particle_id];
		tau = surf->tau;
		vx = surf->vx;
		vy = surf->vy;
		da0 = surf->da0;
		da1 = surf->da1;
		da2 = surf->da2;
		pi00 = surf->pi00;
		pi01 = surf->pi01;
		pi02 = surf->pi02;
		pi11 = surf->pi11;
		pi12 = surf->pi12;
		pi22 = surf->pi22;
		pi33 = surf->pi33;
		temp_r = surf->r;
		sin_temp_phi = surf->sin_phi;
		cos_temp_phi = surf->cos_phi;
		gammaT = surf->gammaT;
		     
		for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
			{
				pT = SPinterp2_pT[ipt];
				for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
				{
					sin_pphi = sin_SPinterp2_pphi[iphi];
					cos_pphi = cos_SPinterp2_pphi[iphi];
					px = pT*cos_pphi;
					py = pT*sin_pphi;
					sin_phi_m_pphi = sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi;
					cos_phi_m_pphi = cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi;
					zvec[1] = temp_r * cos_phi_m_pphi;
					zvec[2] = temp_r * sin_phi_m_pphi;
					for(int ieta=0; ieta < eta_s_npts; ieta++)
					{
						p0 = SPinterp2_p0[ipt][ieta];
						pz = SPinterp2_pz[ipt][ieta];
						f0 = 1./(exp((gammaT*(p0*1. - px*vx - py*vy) - mu)*one_by_Tdec)+sign);	//thermal equilibrium distributions
						
						//viscous corrections
						deltaf = 0.;
						if (use_delta_f)
						{
							deltaf = (1. - sign*f0)*(p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33)*deltaf_prefactor;
						}

						//p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
						S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);
						//ignore points where delta f is large or emission function goes negative from pdsigma
						if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol)) S_p = 0.0;
						S_p_withweight = S_p*tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor accounts for the assumed reflection symmetry along eta direction
						//S_p_withweight = 1.;
						zvec[0] = tau*ch_eta_s[ieta];
						zvec[3] = tau*sh_eta_s[ieta];
						for (int wfi = 0; wfi < n_weighting_functions; wfi++)
							dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi] += S_p_withweight*weight_function(zvec, wfi);
				}
			}
		}
	}

   return;
}*/




//End of file

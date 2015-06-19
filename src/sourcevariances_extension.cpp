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

//#define VERBOSE 2			// specifies level of output - 0 is lowest (no output)
//#define ASSUME_ETA_SYMMETRIC 1	// 1 means integrate only over eta_s = 0..eta_s_max, and multiply by 2 or 0 to speed up calculations
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
		result += S_p_withweight;
	}
	return (result);
}

void SourceVariances::interpolate_FO_loop(FO_surf* FOsurf_ptr)
{
	// use this to set point used for interpolation
	cout << "Setting points to use for freeze-out surface interpolation" << endl;
	set_FOinterp_gridpoints();
	
	// then calculate loop_over_FO_surface(...) at each of these points
	cout << "Computing values of freeze-out surface interpolation" << endl;
	fill_out_FOinterp_grid(FOsurf_ptr);
}

void SourceVariances::set_FOinterp_gridpoints()
{
	// declare arrays of points
	FOI_p0 = new double [FOI_np0pts];
	FOI_pT = new double [FOI_npTpts];
	FOI_pphi = new double [FOI_npphipts];
	FOI_sin_pphi = new double [FOI_npphipts];
	FOI_cos_pphi = new double [FOI_npphipts];
	FOI_pz = new double [FOI_npzpts];
	FOI_mu = new double [FOI_nmupts];
	
	// set limits (in GeV)
	double p0max = 100.0;
	double pTmax = 3.0;
	double pzmax = 100.0;
	double mumax = 0.6;
	double maximumMresonance = 2.5;
	
	
	//***************************************************
	cout << "--> Setting grid points" << endl;
	// fill arrays with appropriate points
	// --> p0: logspace from 0 to p0max
	//logspace(FOI_p0, 0.0, p0max, FOI_np0pts);
	//cout << "--> Finished setting p0 points" << endl;
	// --> pT: logspace from 0 to pTmax
	//logspace(FOI_pT, 0.0, pTmax, FOI_npTpts);
	FOI_pT = SPinterp2_pT;
	cout << "--> Finished setting pT points" << endl;
	// --> pphi: linspace from 0 to 2*pi
	//linspace(FOI_pphi, 0.0, twopi, FOI_npphipts);
	FOI_pphi = SPinterp2_pphi;
	cout << "--> Finished setting pphi points" << endl;
	// --> pz: stratified spacing from -pzmax to pzmax
	//stratify_npts(2.0, pzmax, (FOI_npzpts+3)/2, FOI_npzpts, FOI_pz);
	//cout << "--> Finished setting pz points" << endl;
	// --> mu: linspace from 0 to mumax
	linspace(FOI_mu, 0.0, mumax, FOI_nmupts);
	cout << "--> Finished setting mu points" << endl;
	// --> set sines and cosines of pphi as well
	for (int iFOIpphi = 0; iFOIpphi < FOI_npphipts; iFOIpphi++)
	{
		FOI_cos_pphi[iFOIpphi] = cos(FOI_pphi[iFOIpphi]);
		FOI_sin_pphi[iFOIpphi] = sin(FOI_pphi[iFOIpphi]);
	}
	FOI_eta_s = eta_s;
	cout << "--> Finished setting eta_s points" << endl;
		// --> mu: linspace from 0 to mumax
	linspace(FOI_M, 0.0, maximumMresonance, FOI_nMpts);
	cout << "--> Finished setting mu points" << endl;
	cout << "--> Finished setting all grid points" << endl;
	//***************************************************
	
	
	
	//***************************************************
	// set grid to hold all source variances from freeze-out surface interpolation
	cout << "--> Initializing source variances array" << endl;
	/*FOI_source_variances = new double ***** [n_weighting_functions];
	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		FOI_source_variances[wfi] = new double **** [FOI_np0pts];
		for (int iFOIp0 = 0; iFOIp0 < FOI_np0pts; iFOIp0++)
		{
			FOI_source_variances[wfi][iFOIp0] = new double *** [FOI_npTpts];
			for (int iFOIpT = 0; iFOIpT < FOI_npTpts; iFOIpT++)
			{
				FOI_source_variances[wfi][iFOIp0][iFOIpT] = new double ** [FOI_npphipts];
				for (int iFOIpphi = 0; iFOIpphi < FOI_npphipts; iFOIpphi++)
				{
					FOI_source_variances[wfi][iFOIp0][iFOIpT][iFOIpphi] = new double * [FOI_npzpts];
					for (int iFOIpz = 0; iFOIpz < FOI_npzpts; iFOIpz++)
					{
						FOI_source_variances[wfi][iFOIp0][iFOIpT][iFOIpphi][iFOIpz] = new double [FOI_nmupts];
						for (int iFOImu = 0; iFOImu < FOI_nmupts; iFOImu++)
							FOI_source_variances[wfi][iFOIp0][iFOIpT][iFOIpphi][iFOIpz][iFOImu] = 0.0;
					}
				}
			}
		}
	}*/
	FOI_source_variances = new double ***** [n_weighting_functions];
	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		FOI_source_variances[wfi] = new double **** [FOI_netaspts];
		for (int iFOIetas = 0; iFOIetas < FOI_netaspts; iFOIetas++)
		{
			FOI_source_variances[wfi][iFOIetas] = new double *** [FOI_npTpts];
			for (int iFOIpT = 0; iFOIpT < FOI_npTpts; iFOIpT++)
			{
				FOI_source_variances[wfi][iFOIetas][iFOIpT] = new double ** [FOI_npphipts];
				for (int iFOIpphi = 0; iFOIpphi < FOI_npphipts; iFOIpphi++)
				{
					FOI_source_variances[wfi][iFOIetas][iFOIpT][iFOIpphi] = new double * [FOI_nMpts];
					for (int iFOIM = 0; iFOIM < FOI_nMpts; iFOIM++)
					{
						FOI_source_variances[wfi][iFOIetas][iFOIpT][iFOIpphi][iFOIM] = new double [FOI_nmupts];
						for (int iFOImu = 0; iFOImu < FOI_nmupts; iFOImu++)
							FOI_source_variances[wfi][iFOIetas][iFOIpT][iFOIpphi][iFOIM][iFOImu] = 0.0;
					}
				}
			}
		}
	}
	//***************************************************
	return;
}




void SourceVariances::fill_out_FOinterp_grid(FO_surf* FOsurf_ptr)
{
	time_t rawtime;
	struct tm * timeinfo;
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
	double mu, pi00, pi01, pi02, pi11, pi12, pi22, pi33;
	double temp_r, temp_phi, sin_temp_phi, cos_temp_phi, gammaT;
	double px, py, p0, pT, pz, sin_pphi, cos_pphi, sin_phi_m_pphi, cos_phi_m_pphi;
	double f0, deltaf, S_p, S_p_withweight;
	double symmetry_factor = 1.0;	//symmetry_factor accounts for the assumed reflection symmetry along eta direction
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
	FO_surf* surf;
	double z0, z1, z2, z3;
	double pmudsmu;
	int tracker = 0;
	for(int isurf=0; isurf<FO_length; isurf++)
	{
							//cout << "made it here" << endl;
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
		for (int iFOIpT = 0; iFOIpT < FOI_npTpts; iFOIpT++)
		{
							//cout << "made it here" << endl;
			pT = FOI_pT[iFOIpT];
			for (int iFOIpphi = 0; iFOIpphi < FOI_npphipts; iFOIpphi++)
			{
							//cout << "made it here" << endl;
				//set space-time crap
				cos_pphi = FOI_cos_pphi[iFOIpphi];
				sin_pphi = FOI_sin_pphi[iFOIpphi];
				px = pT*cos_pphi;
				py = pT*sin_pphi;
				sin_phi_m_pphi = sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi;
				cos_phi_m_pphi = cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi;
				z2 = temp_r * sin_phi_m_pphi;
				for (int iFOIp0 = 0; iFOIp0 < FOI_np0pts; iFOIp0++)
				{
							//cout << "made it here" << endl;
					pmudsmu = p0*da0 + px*da1 + py*da2;
					for (int iFOImu = 0; iFOImu < FOI_nmupts; iFOImu++)
					{
							//cout << "made it here" << endl;
						mu = FOI_mu[iFOImu];
						//now get distribution function, emission function, etc.
						//assume sign = +1 for now...
						f0 = 1./(exp((gammaT*(p0*1. - px*vx - py*vy) - mu)*one_by_Tdec)+sign);	//thermal equilibrium distributions
						for (int iFOIpz = 0; iFOIpz < FOI_npzpts; iFOIpz++)
						{
							//cout << "made it here" << endl;
							pz = FOI_pz[iFOIpz];
							//viscous corrections
							deltaf = 0.;
							if (use_delta_f)
								deltaf = (1. - sign*f0)*(p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33)*deltaf_prefactor;
					
							//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
							S_p = prefactor*pmudsmu*f0*(1.+deltaf);
							//ignore points where delta f is large or emission function goes negative from pdsigma
							if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol)) S_p = 0.0;
							S_p_withweight = S_p*tau;
							FOI_source_variances[0][iFOIp0][iFOIpT][iFOIpphi][iFOIpz][iFOImu] += S_p_withweight;				//<1>
							FOI_source_variances[1][iFOIp0][iFOIpT][iFOIpphi][iFOIpz][iFOImu] += S_p_withweight*z2;				//<x_s>
							FOI_source_variances[2][iFOIp0][iFOIpT][iFOIpphi][iFOIpz][iFOImu] += S_p_withweight*z2*z2;			//<x^2_s>
							/*dN_dypTdpTdphi_moments[reso_idx][3][ipt][iphi] += S_p_withweight*z1;			//<x_o>
							dN_dypTdpTdphi_moments[reso_idx][4][ipt][iphi] += S_p_withweight*z1*z1;			//<x^2_o>
							dN_dypTdpTdphi_moments[reso_idx][5][ipt][iphi] += S_p_withweight*z3;			//<x_l>
							dN_dypTdpTdphi_moments[reso_idx][6][ipt][iphi] += S_p_withweight*z3*z3;			//<x^2_l>
							dN_dypTdpTdphi_moments[reso_idx][7][ipt][iphi] += S_p_withweight*z0;			//<t>
							dN_dypTdpTdphi_moments[reso_idx][8][ipt][iphi] += S_p_withweight*z0*z0;			//<t^2>
							dN_dypTdpTdphi_moments[reso_idx][9][ipt][iphi] += S_p_withweight*z2*z1;			//<x_s x_o>
							dN_dypTdpTdphi_moments[reso_idx][10][ipt][iphi] += S_p_withweight*z2*z3;		//<x_s x_l>
							dN_dypTdpTdphi_moments[reso_idx][11][ipt][iphi] += S_p_withweight*z2*z0;		//<x_s t>
							dN_dypTdpTdphi_moments[reso_idx][12][ipt][iphi] += S_p_withweight*z1*z3;		//<x_o x_l>
							dN_dypTdpTdphi_moments[reso_idx][13][ipt][iphi] += S_p_withweight*z1*z0;		//<x_o t>
							dN_dypTdpTdphi_moments[reso_idx][14][ipt][iphi] += S_p_withweight*z3*z0;		//<x_l t>*/
						}
					}
				}
			}
		}
		/*tracker++;
		if (tracker % 1000 == 0)
		{
			time (&rawtime);
			timeinfo = localtime (&rawtime);
			cout << "Finished " << int(double(tracker)/double(FO_length)) << "% by " << asctime(timeinfo);
		}*/
		cout << "Finished isurf = " << isurf << endl;
	}
	return;
}





//End of file

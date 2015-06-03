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

#define USE_PLANE_PSI_ORDER 0		// specifies whether to do HBT relative to flow-plane angle,
					// and at what order: 0 - use plane_psi = 0.0, !0 - use flow-plane angle
					// at given order
//#define use_delta_f 0			// indicates whether to use delta_f corrections to distribution function
					// 0 - false
#define VERBOSE 2			// specifies level of output - 0 is lowest (no output)
#define ASSUME_ETA_SYMMETRIC 1		// 1 means integrate only over eta_s = 0..eta_s_max, and multiply by 2 or 0 to speed up calculations
					// 0 means just integrate over given range of eta_s without worrying about symmetry

using namespace std;

void SourceVariances::Analyze_sourcefunction(FO_surf* FOsurf_ptr)
{
	time_t rawtime;
  	struct tm * timeinfo;

   double plane_psi = 0.0;
   if (VERBOSE > 0) *global_out_stream_ptr << "Determine nth-order plane angles..." << endl;
   //cerr << "Determine nth-order plane angles..." << endl;
   Determine_plane_angle(current_FOsurf_ptr);	//uses only thermal pions...
   int iorder = USE_PLANE_PSI_ORDER;
   if (USE_PLANE_PSI_ORDER)
   {
      if (VERBOSE > 0) *global_out_stream_ptr << "Analyzing source function w.r.t. " << iorder << " th-order participant plane angle..." << endl;
      if (VERBOSE > 0) *global_out_stream_ptr << "psi = " << plane_psi << endl;
      plane_psi = plane_angle[iorder];
   }
   else
   {
      if (VERBOSE > 0) *global_out_stream_ptr << "Analyzing source function w.r.t. psi_0 = " << plane_psi << endl;
      //cerr << "Analyzing source function w.r.t. psi_0 = " << plane_psi << endl;
   }
   global_plane_psi = plane_psi;





  for (int ir = 0; ir <= 0; ir++)	//loop over direct pions and resonances
  {
	if (VERBOSE > 0) *global_out_stream_ptr << "\t\t   * Setting spacetime moments" << endl;
	//cerr << "Setting spacetime moments" << endl;
	Set_dN_dypTdpTdphi_moments(FOsurf_ptr, ir);
	Determine_plane_angle_check(FOsurf_ptr);
   // begin HBT calculations here...
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      if (iKT == 0) continue;	//temporary, need to do something slightly different for K_T == 0...
      if (VERBOSE > 0) *global_out_stream_ptr << "   - Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
      beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
      for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
      {
		if (VERBOSE > 1) *global_out_stream_ptr << "\t --> Calculating K_phi = " << K_phi[iKphi] << " ..." << endl;
      		if (VERBOSE > 1) *global_out_stream_ptr << "\t\t ~ Calculating resonance = " << ir << endl;
		if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Reset emission data" << endl;
		//Reset_EmissionData();
		if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Loading info for this resonance..." << endl;
		Load_resonance_info(ir, K_T[iKT], K_phi[iKphi]);
		if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Loaded resonance info" << endl;
		//ir == 0 corresponds to just thermal pions
		//otherwise, gives irth resonance contribution
		if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Getting phase-space integrated spacetime moments" << endl;
		if (ir > 0)	//don't need to do any resonance integrals for thermal pions!
			Do_resonance_integrals(current_FOsurf_ptr, iKT, iKphi, ir);
		if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Finished calculations for reso_idx = " << ir << endl
			<< "\t\t --> just updating the contributions from this resonance to each source integral now..." << endl;
		Update_source_variances(iKT, iKphi, ir);
		if (VERBOSE > 1) *global_out_stream_ptr << "\t --> Finished!  Getting R^2_s..." << endl;
		Calculate_R2_side(iKT, iKphi);
		if (VERBOSE > 1) *global_out_stream_ptr << "\t --> Moving on!" << endl;
	return;
      }
   }
  }

   return;
}

void SourceVariances::Analyze_sourcefunction_check(FO_surf* FOsurf_ptr)
{
   double plane_psi = 0.0;
   global_plane_psi = plane_psi;

Set_dN_dypTdpTdphi_moments(FOsurf_ptr, 0);	//ir --> 0 for pions
Determine_plane_angle_check(FOsurf_ptr);
if (1) return;

   // begin HBT calculations here...
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      if (iKT == 0) continue;	//temporary
      double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
      beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
      for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
      {
		for (int ir = 7; ir <= 7; ir++)	//loop over direct pions and resonances
		{
			cout << "Reset emission data" << endl;
			Reset_EmissionData();
			Load_resonance_info(ir, K_T[iKT], K_phi[iKphi]);
			cout << "Loaded resonance info" << endl;
			//SetEmissionData only computes direct emission function of given particles
			//0 in 2nd argument corresponds to just thermal pions
			//otherwise, gives irth resonance contribution
			if (ir == 0)
				SetEmissionData(current_FOsurf_ptr, K_T[iKT], K_phi[iKphi]);
			else
			{
				cout << "Setting spacetime moments" << endl;
				Set_dN_dypTdpTdphi_moments(FOsurf_ptr, ir);
				//Do_resonance_integrals(current_FOsurf_ptr, K_T[iKT], K_phi[iKphi], ir);
				cout << "Getting phase-space integrated spacetime moments" << endl;
				Do_resonance_integrals(current_FOsurf_ptr, iKT, iKphi, ir);
				cout << "Finished calculations for reso_idx = " << ir << endl;
			}
			Update_source_variances(iKT, iKphi, ir);
		}
		Calculate_R2_side(iKT, iKphi);
      }
   }
   return;
}

void SourceVariances::Set_dN_dypTdpTdphi_moments(FO_surf* FOsurf_ptr, int reso_idx)
{
double mass = particle_mass;
if (reso_idx > 0) mass = current_resonance_mass;
//double mass = current_resonance_mass;
   for(int i=0; i<eta_s_npts; i++)
   {
       double local_eta_s = eta_s[i];
       double local_cosh = cosh(SP_p_y - local_eta_s);
       double local_sinh = sinh(SP_p_y - local_eta_s);
	for(int ipx=0; ipx<n_interp1_px_pts; ipx++)
	for(int ipy=0; ipy<n_interp1_py_pts; ipy++)
       {
	double mT = sqrt(mass*mass + SPinterp1_px[ipx]*SPinterp1_px[ipx] + SPinterp1_py[ipy]*SPinterp1_py[ipy]);
          SPinterp1_p0[ipx][ipy][i] = mT*local_cosh;
          SPinterp1_pz[ipx][ipy][i] = mT*local_sinh;
       }
	for(int ipt=0; ipt<n_interp2_pT_pts; ipt++)
       {
		double mT = sqrt(mass*mass + SPinterp2_pT[ipt]*SPinterp2_pT[ipt]);
		SPinterp2_p0[ipt][i] = mT*local_cosh;
		SPinterp2_pz[ipt][i] = mT*local_sinh;
       }
   }
	if (INTERPOLATION_FORMAT == 1)	//using cartesian grid for interpolation (px, py)
		Cal_dN_dypTdpTdphi_with_weights_cartesian(FOsurf_ptr, reso_idx);
	else if (INTERPOLATION_FORMAT == 2)	//using polar grid for interpolation (pT, pphi)
		Cal_dN_dypTdpTdphi_with_weights_polar(FOsurf_ptr, reso_idx);

   return;
}

void SourceVariances::Cal_dN_dypTdpTdphi_with_weights_cartesian(FO_surf* FOsurf_ptr, int reso_idx)
{
//CURRENTLY USING WRONG DEFINITIONS OF SIGN AND DEGEN
//NEED TO FIX BEFORE VERSION IS STABLE
   double sign = particle_sign;
   double degen = particle_gspin;
   double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
//assume reso_idx >= 1
//double mu = 0.0;

   for(int isurf=0; isurf<FO_length ; isurf++)
   {
      FO_surf* surf = &FOsurf_ptr[isurf];
      double mu = surf->particle_mu[particle_id];
      double tau = surf->tau;
      double vx = surf->vx;
      double vy = surf->vy;
      double Tdec = surf->Tdec;
      double Pdec = surf->Pdec;
      double Edec = surf->Edec;
      double da0 = surf->da0;
      double da1 = surf->da1;
      double da2 = surf->da2;
      double pi00 = surf->pi00;
      double pi01 = surf->pi01;
      double pi02 = surf->pi02;
      double pi11 = surf->pi11;
      double pi12 = surf->pi12;
      double pi22 = surf->pi22;
      double pi33 = surf->pi33;
	
	double x = surf->xpt;
	double y = surf->ypt;
	double temp_r = sqrt(x*x+y*y);
	double temp_phi = atan2(y,x);

      double vT = sqrt(vx*vx + vy*vy);
      double gammaT = 1./sqrt(1. - vT*vT);

      double deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
      
      //for(int ipt = 0; ipt < n_SP_pT; ipt++)
      for(int ipx = 0; ipx < n_interp1_px_pts; ipx++)
      {
      //for(int iphi = 0; iphi < n_SP_pphi; iphi++)
      for(int ipy = 0; ipy < n_interp1_py_pts; ipy++)
      {
         double px = SPinterp1_px[ipx];
         double py = SPinterp1_py[ipy];
	double pphi = atan2(py, px);
	//if (pphi < 0.) pphi += 2.*M_PI;
      for(int ieta=0; ieta < eta_s_npts; ieta++)
      {
         //double p0 = SP_p0[ipt][ieta];
         //double pz = SP_pz[ipt][ieta];
	double p0 = SPinterp1_p0[ipx][ipy][ieta];
	double pz = SPinterp1_pz[ipx][ipy][ieta];
         double expon = (gammaT*(p0*1. - px*vx - py*vy) - mu)/Tdec;
         double f0 = 1./(exp(expon)+sign);
         //if(expon > 20) f0 = 0.0e0;
         //else f0 = 1./(exp(expon)+sign);       //thermal equilibrium distributions

         //p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
         double pdsigma = p0*da0 + px*da1 + py*da2;

         //viscous corrections
         double Wfactor = p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
         double deltaf = 0.;
	 if (use_delta_f)
	 {
		deltaf = (1. - sign*f0)*Wfactor*deltaf_prefactor;
	 }

         double S_p = prefactor*pdsigma*f0*(1.+deltaf);
	double symmetry_factor = 1.0;
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
	//ignore points where delta f is large
	 if (1. + deltaf < 0.0) S_p = 0.0;
        if (flagneg == 1 && S_p < tol)
        {//neglect points where emission function goes negative from pdsigma
           S_p = 0.0e0;
        }
         double S_p_withweight = S_p*tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor accounts for the assumed reflection symmetry along eta direction
	double sin_phi = sin(temp_phi - pphi);
	double cos_phi = cos(temp_phi - pphi);
	zvec[0] = tau*ch_eta_s[ieta];
	zvec[1] = temp_r * cos_phi;
	zvec[2] = temp_r * sin_phi;
	zvec[3] = tau*sh_eta_s[ieta];
	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		dN_dypTdpTdphi_moments[reso_idx][wfi][ipx][ipy] += S_p_withweight*weight_function(zvec, wfi);
      }
      }
      }
   }
	for(int ipx = 0; ipx < n_interp1_px_pts; ipx++)
	for(int ipy = 0; ipy < n_interp1_py_pts; ipy++)
	{
	cout << SPinterp1_px[ipx] << "  " << SPinterp1_py[ipy] << "  " << dN_dypTdpTdphi_moments[reso_idx][1][ipx][ipy] << endl;
	for(int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		double temp = dN_dypTdpTdphi_moments[reso_idx][wfi][ipx][ipy];
		ln_dN_dypTdpTdphi_moments[reso_idx][wfi][ipx][ipy] = log(abs(temp));
		sign_of_dN_dypTdpTdphi_moments[reso_idx][wfi][ipx][ipy] = sgn(temp);
	}
	}
   return;
}

void SourceVariances::Cal_dN_dypTdpTdphi_with_weights_polar(FO_surf* FOsurf_ptr, int reso_idx)
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
Pdec = (&FOsurf_ptr[0])->Pdec;
Edec = (&FOsurf_ptr[0])->Edec;
//cerr << "DEBUG: " << Tdec << "  " << Pdec << "  " << Edec << endl;
double deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

   for(int isurf=0; isurf<FO_length ; isurf++)
   {
      FO_surf* surf = &FOsurf_ptr[isurf];
      double mu = surf->particle_mu[particle_id];
      double tau = surf->tau;
      double vx = surf->vx;
      double vy = surf->vy;
      double da0 = surf->da0;
      double da1 = surf->da1;
      double da2 = surf->da2;
      double pi00 = surf->pi00;
      double pi01 = surf->pi01;
      double pi02 = surf->pi02;
      double pi11 = surf->pi11;
      double pi12 = surf->pi12;
      double pi22 = surf->pi22;
      double pi33 = surf->pi33;
	double temp_r = surf->r;
	double temp_phi = surf->phi;
	double sin_temp_phi = sin(temp_phi);
	double cos_temp_phi = cos(temp_phi);
      double gammaT = surf->gammaT;
      
      for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
      {
	double pT = SPinterp2_pT[ipt];
      for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
      {
	double sin_pphi = sin_SPinterp2_pphi[iphi];
	double cos_pphi = cos_SPinterp2_pphi[iphi];
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;
	//double pphi = atan2(py, px);
	//if (pphi < 0.) pphi += 2.*M_PI;
	double sin_phi_m_pphi = sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi;
	double cos_phi_m_pphi = cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi;
	zvec[1] = temp_r * cos_phi_m_pphi;
	zvec[2] = temp_r * sin_phi_m_pphi;
      for(int ieta=0; ieta < eta_s_npts; ieta++)
      {
         double p0 = SPinterp2_p0[ipt][ieta];
         double pz = SPinterp2_pz[ipt][ieta];
         double expon = (gammaT*(p0*1. - px*vx - py*vy) - mu)/Tdec;
         double f0 = 1./(exp(expon)+sign);
         //if(expon > 20) f0 = 0.0e0;
         //else f0 = 1./(exp(expon)+sign);       //thermal equilibrium distributions

         //p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
         double pdsigma = p0*da0 + px*da1 + py*da2;

         //viscous corrections
         double Wfactor = p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
         double deltaf = 0.;
	 if (use_delta_f)
	 {
		deltaf = (1. - sign*f0)*Wfactor*deltaf_prefactor;
	 }

         double S_p = prefactor*pdsigma*f0*(1.+deltaf);
	double symmetry_factor = 1.0;
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
	//ignore points where delta f is large
	 if (1. + deltaf < 0.0) S_p = 0.0;
        if (flagneg == 1 && S_p < tol)
        {//neglect points where emission function goes negative from pdsigma
           S_p = 0.0e0;
        }
         double S_p_withweight = S_p*tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor accounts for the assumed reflection symmetry along eta direction
	zvec[0] = tau*ch_eta_s[ieta];
	zvec[3] = tau*sh_eta_s[ieta];
	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi] += S_p_withweight*weight_function(zvec, wfi);
      }
      }
      }
   }
//set log of dN_dypTdpTdphi_moments...
	for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
	for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
	for(int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		//ln_dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi] = log(dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi]);
		double temp = dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi];
		ln_dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi] = log(abs(temp));
		sign_of_dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi] = sgn(temp);
	}
   return;
}


void SourceVariances::Cal_dN_dypTdpTdphi_interpolate_cartesian_grid(double** SP_px, double** SP_py)
{
	//cout << "Starting interpolated values: " << endl;
	for(int ipt = 0; ipt < n_SP_pT; ipt++)
	for(int iphi = 0; iphi < n_SP_pphi; iphi++)
	{
		double px = SP_px[ipt][iphi];
		double py = SP_py[ipt][iphi];
		//dN_dypTdpTdphi[ipt][iphi] = exp(interpolate2D(SPinterp1_px, SPinterp1_py, dN_dypTdpTdphi_moments[0][0], px, py, n_interp1_px_pts, n_interp1_py_pts, 1, UNIFORM_SPACING));
		dN_dypTdpTdphi[ipt][iphi] = interpolate2D(SPinterp1_px, SPinterp1_py, dN_dypTdpTdphi_moments[0][0], px, py, n_interp1_px_pts, n_interp1_py_pts, 1, UNIFORM_SPACING);
		//dN_dypTdpTdphi[ipt][iphi] = exp(interpBiPolyDirect(SPinterp1_px, SPinterp1_py, ln_dN_dypTdpTdphi_moments[0][0], px, py, n_interp1_px_pts, n_interp1_py_pts));
		//dN_dypTdpTdphi[ipt][iphi] = exp(interpBiLinearDirect(SPinterp1_px, SPinterp1_py, ln_dN_dypTdpTdphi_moments[0][0], px, py, n_interp1_px_pts, n_interp1_py_pts));
		//dN_dypTdpTdphi[ipt][iphi] = exp(interpBiCubicDirect(SPinterp1_px, SPinterp1_py, ln_dN_dypTdpTdphi_moments[0][0], px, py, n_interp1_px_pts, n_interp1_py_pts));
		//dN_dypTdpTdphi[ipt][iphi] = interpBiPolyDirect(SPinterp1_px, SPinterp1_py, dN_dypTdpTdphi_moments[0][0], px, py, n_interp1_px_pts, n_interp1_py_pts);
		//dN_dypTdpTdphi[ipt][iphi] = interpBiLinearDirect(SPinterp1_px, SPinterp1_py, dN_dypTdpTdphi_moments[0][0], px, py, n_interp1_px_pts, n_interp1_py_pts);
		//dN_dypTdpTdphi[ipt][iphi] = interpBiCubicDirect(SPinterp1_px, SPinterp1_py, dN_dypTdpTdphi_moments[0][0], px, py, n_interp1_px_pts, n_interp1_py_pts);
		//cout << SP_pT[ipt] << "  " << SP_pphi[iphi] << "  " << px << "  " << py << "  " << dN_dypTdpTdphi[ipt][iphi] << endl;
	}
	//cout << endl;
	//cout << "Starting gridded values: " << endl;
	//for(int ipx = 0; ipx < n_interp1_px_pts; ipx++)
	//for(int ipy = 0; ipy < n_interp1_py_pts; ipy++)
	//	cout << SPinterp1_px[ipx] << "  " << SPinterp1_py[ipy] << "  " << dN_dypTdpTdphi_moments[0][0][ipx][ipy] << endl;

   return;
}


void SourceVariances::Cal_dN_dypTdpTdphi_interpolate_polar_grid(double* SP_pT, double* SP_pphi)
{
	//cout << "Starting interpolated values: " << endl;
	for(int ipt = 0; ipt < n_SP_pT; ipt++)
	{
	for(int iphi = 0; iphi < n_SP_pphi; iphi++)
	{
		double pT = SP_pT[ipt];
		double pphi = SP_pphi[iphi];
		dN_dypTdpTdphi[ipt][iphi] = exp(interpolate2D(SPinterp2_pT, SPinterp2_pphi, ln_dN_dypTdpTdphi_moments[0][0], pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING));
		//dN_dypTdpTdphi[ipt][iphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[0][0], pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING);
		//double temp_xs_S = sign_of_dN_dypTdpTdphi_moments[0][2][ipt][iphi]*exp(interpolate2D(SPinterp2_pT, SPinterp2_pphi, ln_dN_dypTdpTdphi_moments[0][1], pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING));
		//double temp_xs_S_alt = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[0][1], pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING);
		//double temp_x2s_S = sign_of_dN_dypTdpTdphi_moments[0][2][ipt][iphi]*exp(interpolate2D(SPinterp2_pT, SPinterp2_pphi, ln_dN_dypTdpTdphi_moments[0][2], pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING));
		//double temp_x2s_S_alt = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[0][2], pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING);
		cout << pT << "  " << pphi << "  " << pT*cos(pphi) << "  " << pT*sin(pphi) << "  " << dN_dypTdpTdphi[ipt][iphi] << endl;
	}
	cout << endl;
	}
	//cout << endl;
	//cout << "Starting gridded values: " << endl;
	//for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
	//{
	//	for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
	//		cout << SPinterp2_pT[ipt] << "  " << SPinterp2_pphi[iphi] << "  " << dN_dypTdpTdphi_moments[0][0][ipt][iphi] << endl;
	//	cout << endl;
	//}

   return;
}



void SourceVariances::Determine_plane_angle(FO_surf* FOsurf_ptr)
{
   double mass = particle_mass;
   double* mT = new double [n_SP_pT];
   double** px = new double* [n_SP_pT];
   double** py = new double* [n_SP_pT];
   double** p0 = new double* [n_SP_pT];
   double** pz = new double* [n_SP_pT];
   for(int ipt=0; ipt<n_SP_pT; ipt++)
   {
      px[ipt] = new double [n_SP_pphi];
      py[ipt] = new double [n_SP_pphi];
      p0[ipt] = new double [eta_s_npts];
      pz[ipt] = new double [eta_s_npts];
   }
   
   for(int ipt=0; ipt<n_SP_pT; ipt++)
      mT[ipt] = sqrt(mass*mass + SP_pT[ipt]*SP_pT[ipt]);
   for(int iphi = 0; iphi<n_SP_pphi; iphi++)
   {
      double cos_phi = cos(SP_pphi[iphi]);
      double sin_phi = sin(SP_pphi[iphi]);
      for(int ipt=0; ipt<n_SP_pT; ipt++)
      {
         px[ipt][iphi] = SP_pT[ipt]*cos_phi;
         py[ipt][iphi] = SP_pT[ipt]*sin_phi;
      }
   }

   for(int i=0; i<eta_s_npts; i++)
   {
       double local_eta_s = eta_s[i];
       double local_cosh = cosh(SP_p_y - local_eta_s);
       double local_sinh = sinh(SP_p_y - local_eta_s);
       for(int ipt=0; ipt<n_SP_pT; ipt++)
       {
          p0[ipt][i] = mT[ipt]*local_cosh;
          pz[ipt][i] = mT[ipt]*local_sinh;
       }
   }

Cal_dN_dypTdpTdphi(p0, px, py, pz, FOsurf_ptr);

   for(int ipt=0; ipt<n_SP_pT; ipt++)
   {
	for(int iphi=0; iphi<n_SP_pphi; iphi++)
	{
		dN_dydphi[iphi] += dN_dypTdpTdphi[ipt][iphi]*SP_pT[ipt]*SP_pT_weight[ipt];
		pTdN_dydphi[iphi] += dN_dypTdpTdphi[ipt][iphi]*SP_pT[ipt]*SP_pT[ipt]*SP_pT_weight[ipt];
		dN_dypTdpT[ipt] += dN_dypTdpTdphi[ipt][iphi]*SP_pphi_weight[iphi];
	}
   }
   
   double norm = 0.0e0;
   for(int iphi=0; iphi<n_SP_pphi; iphi++)
      norm += dN_dydphi[iphi]*SP_pphi_weight[iphi];
   for(int iorder=0; iorder < n_order; iorder++)
   {
      double cosine = 0.0e0;
      double sine = 0.0e0;
      for(int iphi=0; iphi<n_SP_pphi; iphi++)
      {
         cosine += dN_dydphi[iphi]*cos(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
         sine += dN_dydphi[iphi]*sin(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
         //get pT-differential v_n here
         for(int ipt=0; ipt<n_SP_pT; ipt++)
         {
             cosine_iorder[ipt][iorder] += dN_dypTdpTdphi[ipt][iphi]*cos(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
             sine_iorder[ipt][iorder] += dN_dypTdpTdphi[ipt][iphi]*sin(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
         }
      }
      for(int ipt=0; ipt<n_SP_pT; ipt++)
      {
         cosine_iorder[ipt][iorder] /= dN_dypTdpT[ipt];
         sine_iorder[ipt][iorder] /= dN_dypTdpT[ipt];
      }
      cosine = cosine/norm;
      sine = sine/norm;
      if( sqrt(sine*sine + cosine*cosine) < 1e-8)
         plane_angle[iorder] = 0.0e0;
      else
         plane_angle[iorder] = atan2(sine, cosine)/double(iorder);
   }

//cout << "Currently getting <p_T> stuff..." << endl;

   mean_pT = 0.;
   for(int iphi=0; iphi<n_SP_pphi; iphi++)
      mean_pT += pTdN_dydphi[iphi]*SP_pphi_weight[iphi];
   mean_pT /= norm;
   plane_angle[0] = norm;

   delete[] mT;
   for(int ipt=0; ipt<n_SP_pT; ipt++)
   {
      delete[] px[ipt];
      delete[] py[ipt];
      delete[] p0[ipt];
      delete[] pz[ipt];
   }
   delete[] px;
   delete[] py;
   delete[] p0;
   delete[] pz;

   return;
}

void SourceVariances::Determine_plane_angle_check(FO_surf* FOsurf_ptr)
{
   double mass = particle_mass;
   double* mT = new double [n_SP_pT];
   double** px = new double* [n_SP_pT];
   double** py = new double* [n_SP_pT];
   double** p0 = new double* [n_SP_pT];
   double** pz = new double* [n_SP_pT];
   for(int ipt=0; ipt<n_SP_pT; ipt++)
   {
      px[ipt] = new double [n_SP_pphi];
      py[ipt] = new double [n_SP_pphi];
      p0[ipt] = new double [eta_s_npts];
      pz[ipt] = new double [eta_s_npts];
   }
   
   for(int ipt=0; ipt<n_SP_pT; ipt++)
//{
      mT[ipt] = sqrt(mass*mass + SP_pT[ipt]*SP_pT[ipt]);
//}
   for(int iphi = 0; iphi<n_SP_pphi; iphi++)
   {
      double cos_phi = cos(SP_pphi[iphi]);
      double sin_phi = sin(SP_pphi[iphi]);
      for(int ipt=0; ipt<n_SP_pT; ipt++)
      {
         px[ipt][iphi] = SP_pT[ipt]*cos_phi;
         py[ipt][iphi] = SP_pT[ipt]*sin_phi;
      }
   }

   for(int i=0; i<eta_s_npts; i++)
   {
       double local_eta_s = eta_s[i];
       double local_cosh = cosh(SP_p_y - local_eta_s);
       double local_sinh = sinh(SP_p_y - local_eta_s);
       for(int ipt=0; ipt<n_SP_pT; ipt++)
       {
          p0[ipt][i] = mT[ipt]*local_cosh;
          pz[ipt][i] = mT[ipt]*local_sinh;
       }
   }

   //Cal_dN_dypTdpTdphi(p0, px, py, pz, FOsurf_ptr);
	if (INTERPOLATION_FORMAT == 1)	//using cartesian grid for interpolation (px, py)
		Cal_dN_dypTdpTdphi_interpolate_cartesian_grid(px, py);
	else if (INTERPOLATION_FORMAT == 2)	//using polar grid for interpolation (pT, pphi)
   		Cal_dN_dypTdpTdphi_interpolate_polar_grid(SP_pT, SP_pphi);

   for(int ipt=0; ipt<n_SP_pT; ipt++)
   {
	for(int iphi=0; iphi<n_SP_pphi; iphi++)
	{
		dN_dydphi[iphi] += dN_dypTdpTdphi[ipt][iphi]*SP_pT[ipt]*SP_pT_weight[ipt];
		pTdN_dydphi[iphi] += dN_dypTdpTdphi[ipt][iphi]*SP_pT[ipt]*SP_pT[ipt]*SP_pT_weight[ipt];
		dN_dypTdpT[ipt] += dN_dypTdpTdphi[ipt][iphi]*SP_pphi_weight[iphi];
	}
   }
   
   double norm = 0.0e0;
   for(int iphi=0; iphi<n_SP_pphi; iphi++)
      norm += dN_dydphi[iphi]*SP_pphi_weight[iphi];
   for(int iorder=0; iorder < n_order; iorder++)
   {
      double cosine = 0.0e0;
      double sine = 0.0e0;
      for(int iphi=0; iphi<n_SP_pphi; iphi++)
      {
         cosine += dN_dydphi[iphi]*cos(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
         sine += dN_dydphi[iphi]*sin(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
         //get pT-differential v_n here
         for(int ipt=0; ipt<n_SP_pT; ipt++)
         {
             cosine_iorder[ipt][iorder] += dN_dypTdpTdphi[ipt][iphi]*cos(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
             sine_iorder[ipt][iorder] += dN_dypTdpTdphi[ipt][iphi]*sin(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
         }
      }
      for(int ipt=0; ipt<n_SP_pT; ipt++)
      {
         cosine_iorder[ipt][iorder] /= dN_dypTdpT[ipt];
         sine_iorder[ipt][iorder] /= dN_dypTdpT[ipt];
      }
      cosine = cosine/norm;
      sine = sine/norm;
      if( sqrt(sine*sine + cosine*cosine) < 1e-8)
         plane_angle[iorder] = 0.0e0;
      else
         plane_angle[iorder] = atan2(sine, cosine)/double(iorder);
   }

//cout << "Currently getting <p_T> stuff..." << endl;

   mean_pT = 0.;
   for(int iphi=0; iphi<n_SP_pphi; iphi++)
      mean_pT += pTdN_dydphi[iphi]*SP_pphi_weight[iphi];
   mean_pT /= norm;
   plane_angle[0] = norm;

   delete[] mT;
   for(int ipt=0; ipt<n_SP_pT; ipt++)
   {
      delete[] px[ipt];
      delete[] py[ipt];
      delete[] p0[ipt];
      delete[] pz[ipt];
   }
   delete[] px;
   delete[] py;
   delete[] p0;
   delete[] pz;

   return;
}


//void SourceVariances::Do_resonance_integrals(FO_surf* FOsurf_ptr, double K_T_local, double K_phi_local, int reso_idx)
void SourceVariances::Do_resonance_integrals(FO_surf* FOsurf_ptr, int iKT, int iKphi, int reso_idx)
{
	//Mres = current_resonance_mass;
	//mass = particle_mass;
	//Gamma = current_resonance_Gamma;
	//br = current_resonance_total_br;
	//m2 = current_resonance_decay_masses[0];
	//m3 = current_resonance_decay_masses[1];
	//mT = sqrt(mass*mass + K_T_local*K_T_local);
	//current_K_phi = K_phi_local;
	double Kx = K_T[iKT]*cos(K_phi[iKphi]);
	double Ky = K_T[iKT]*sin(K_phi[iKphi]);
	set_surfpts();
	int idx = 0;
	if (abs(m3) <= 1.e-6)
		n_body = 2;
	else
		n_body = 3;
	//NO ETA_S LOOP --> ALREADY DID IT IN COMPUTING DN_DYPTDPTDPHI_WITH_WEIGHTS
	//for(int i=0; i<eta_s_npts; i++)
	//for(int i=0; i<1; i++)
	//{
		//double local_eta_s = eta_s[i];
		//double local_eta_s_wt = eta_s_weight[i];
		//p_y = local_eta_s;
		//ch_p_y = ch_eta_s[i];
		//sh_p_y = sh_eta_s[i];
		//double ch_localetas = cosh(local_eta_s);
		//double sh_localetas = sinh(local_eta_s);
		//double ch_localetas = ch_eta_s[i];
		//double sh_localetas = sh_eta_s[i];

		//double K0 = mT*(ch_K_y*ch_localetas - sh_K_y*sh_localetas);
		//double Kz = mT*(sh_K_y*ch_localetas - ch_K_y*sh_localetas);

		//for (int iweight = 0; iweight < n_weighting_functions; iweight++)
		//	source_variances_array[iweight] += local_eta_s_wt*Mres*do_all_integrals(iweight, i);
	do_all_integrals(iKT, iKphi, reso_idx);
	//}
//exit(1);

	return;
}

void SourceVariances::SetEmissionData(FO_surf* FOsurf_ptr, double K_T_local, double K_phi_local)
{
	double mass = particle_mass;
	double mT = sqrt(mass*mass + K_T_local*K_T_local);

	double px = K_T_local*cos(K_phi_local);
	double py = K_T_local*sin(K_phi_local);

	int idx = 0;
	for(int i=0; i<eta_s_npts; i++)
	{
		double local_eta_s = eta_s[i];
		//double ch_localetas = cosh(local_eta_s);
		//double sh_localetas = sinh(local_eta_s);
		double ch_localetas = ch_eta_s[i];
		double sh_localetas = sh_eta_s[i];

		//double p0 = mT*cosh(K_y-local_eta_s);
		//double pz = mT*sinh(K_y-local_eta_s);
		double p0 = mT*(ch_K_y*ch_localetas - sh_K_y*sh_localetas);
		double pz = mT*(sh_K_y*ch_localetas - ch_K_y*sh_localetas);

		for (int j = 0; j < FO_length; j++)
		{
			//Now that the data is loaded, cycle through it to find the freeze out surface and the emission function.
			double S_p = 0.0e0;
			S_p = Emissionfunction(p0, px, py, pz, &FOsurf_ptr[j]);
			if (flagneg == 1 && S_p < tol)
			{
				S_p = 0.0e0;
			}
			//	  else
			//        {
			double S_p_withweight = S_p*FOsurf_ptr[j].tau*eta_s_weight[i];
			(*Emissionfunction_ptr)[idx].data = S_p_withweight; 
			(*Emissionfunction_ptr)[idx].t = FOsurf_ptr[j].tau*ch_localetas;
			(*Emissionfunction_ptr)[idx].tau = FOsurf_ptr[j].tau;
			(*Emissionfunction_ptr)[idx].eta = local_eta_s;
			(*Emissionfunction_ptr)[idx].x = FOsurf_ptr[j].xpt;
			(*Emissionfunction_ptr)[idx].y = FOsurf_ptr[j].ypt;
			(*Emissionfunction_ptr)[idx].r = sqrt(FOsurf_ptr[j].xpt*FOsurf_ptr[j].xpt + FOsurf_ptr[j].ypt*FOsurf_ptr[j].ypt);
			double temp_phi = atan2(FOsurf_ptr[j].ypt, FOsurf_ptr[j].xpt);
			(*Emissionfunction_ptr)[idx].phi = temp_phi;
			(*Emissionfunction_ptr)[idx].z = FOsurf_ptr[j].tau*sh_localetas;
			idx++;
			//        }
		}
	}

	Emissionfunction_length = idx;  //only want non-zero values of emission function

	return;
}

void SourceVariances::Load_resonance_info(int reso_idx, double K_T_local, double K_phi_local)
{
	current_resonance_idx = reso_idx;
	if (reso_idx == 0)
	{
		muRES = particle_mu;
		signRES = particle_sign;
		gRES = particle_gspin;
		return;
	}
	else
	{
		current_resonance_mass = resonances.resonance_mass[reso_idx-1];
		current_resonance_Gamma = resonances.resonance_Gamma[reso_idx-1];
		current_resonance_total_br = resonances.resonance_total_br[reso_idx-1];
		current_resonance_decay_masses[0] = resonances.resonance_decay_masses[reso_idx-1][0];
		current_resonance_decay_masses[1] = resonances.resonance_decay_masses[reso_idx-1][1];

		muRES = 0.0;
		//signRES and gRES need to be set for each resonance!
		signRES = particle_sign;
		gRES = particle_gspin;

		Mres = current_resonance_mass;
		mass = particle_mass;
		Gamma = current_resonance_Gamma;
		br = current_resonance_total_br;
		m2 = current_resonance_decay_masses[0];
		m3 = current_resonance_decay_masses[1];
		mT = sqrt(mass*mass + K_T_local*K_T_local);
		pT = K_T_local;
		current_K_phi = K_phi_local;
		cos_cKphi = cos(K_phi_local);
		sin_cKphi = sin(K_phi_local);

		//VEC_tau_factor = new double [n_tau_pts];
		//for (int itau = 0; itau < n_tau_pts; itau++)
		//	VEC_tau_factor[itau] = tau_wts[reso_idx-1][itau]*Gamma*exp(-Gamma*tau_pts[reso_idx-1][itau]);

		//set up vectors of points to speed-up integrals...
		VEC_pstar = new double [n_s_pts];
		VEC_Estar = new double [n_s_pts];
		VEC_DeltaY = new double [n_s_pts];
		VEC_g_s = new double [n_s_pts];
		VEC_s_factor = new double [n_s_pts];
		VEC_v_factor = new double * [n_s_pts];
		VEC_zeta_factor = new double ** [n_s_pts];
		VEC_Yp = new double [n_s_pts];
		VEC_Ym = new double [n_s_pts];
		VEC_P_Y = new double * [n_s_pts];
		VEC_MTbar = new double * [n_s_pts];
		VEC_DeltaMT = new double * [n_s_pts];
		VEC_MTp = new double * [n_s_pts];
		VEC_MTm = new double * [n_s_pts];
		VEC_MT = new double ** [n_s_pts];
		VEC_PPhi_tilde = new double ** [n_s_pts];
		VEC_Pp = new double *** [n_s_pts];
		//cerr << "Entering loops now..." << endl;
		for (int is = 0; is < n_s_pts; is++)
		{
			//cerr << "In s loop# = " << is << endl;
			double s_loc = s_pts[reso_idx-1][is];
			double g_s_loc = g(s_loc);
			VEC_g_s[is] = g_s_loc;
			VEC_s_factor[is] = s_wts[reso_idx-1][is]*g_s_loc;
			double pstar_loc = sqrt(((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc)/(2.0*Mres));
			VEC_pstar[is] = pstar_loc;
			double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
			VEC_Estar[is] = Estar_loc;
			double psBmT = pstar_loc / mT;
			double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
			VEC_DeltaY[is] = DeltaY_loc;
			//VEC_PPhi_tilde[is] = new double ** [eta_s_npts];
			//VEC_Pp[is] = new double *** [eta_s_npts];
			/*DEBUG*///cout << Mres << "     " << mass << "     "
			/*DEBUG*///	<< ((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc)/(2.0*Mres)
			/*DEBUG*///	<< "     " << s_loc << "     " << pstar_loc << endl;
			p_y = 0.0;
			VEC_v_factor[is] = new double [n_v_pts];
			VEC_zeta_factor[is] = new double * [n_v_pts];
			VEC_Yp[is] = p_y + DeltaY_loc;
			VEC_Ym[is] = p_y - DeltaY_loc;
			VEC_P_Y[is] = new double [n_v_pts];
			VEC_MTbar[is] = new double [n_v_pts];
			VEC_DeltaMT[is] = new double [n_v_pts];
			VEC_MTp[is] = new double [n_v_pts];
			VEC_MTm[is] = new double [n_v_pts];
			VEC_MT[is] = new double * [n_v_pts];
			VEC_PPhi_tilde[is] = new double * [n_v_pts];
			VEC_Pp[is] = new double ** [n_v_pts];
			for(int iv = 0; iv < n_v_pts; iv++)
			{
				//cerr << "In v loop# = " << iv << endl;
				double v_loc = v_pts[iv];
				double P_Y_loc = p_y + v_loc*DeltaY_loc;
				VEC_P_Y[is][iv] = P_Y_loc;
				double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
				double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
				VEC_v_factor[is][iv] = v_wts[iv]*DeltaY_loc/x2;
				double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
				VEC_MTbar[is][iv] = MTbar_loc;
				double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
				VEC_DeltaMT[is][iv] = DeltaMT_loc;
				VEC_MTp[is][iv] = MTbar_loc + DeltaMT_loc;
				VEC_MTm[is][iv] = MTbar_loc - DeltaMT_loc;
				VEC_MT[is][iv] = new double [n_zeta_pts];
				VEC_PPhi_tilde[is][iv] = new double [n_zeta_pts];
				VEC_Pp[is][iv] = new double * [n_zeta_pts];
				VEC_zeta_factor[is][iv] = new double [n_zeta_pts];
				for(int izeta = 0; izeta < n_zeta_pts; izeta++)
				{
					double zeta_loc = zeta_pts[izeta];
					double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
					VEC_MT[is][iv][izeta] = MT_loc;
					VEC_zeta_factor[is][iv][izeta] = zeta_wts[izeta]*MT_loc;
					double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
					double PPhi_tilde_loc = acos( (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc) );
					VEC_PPhi_tilde[is][iv][izeta] = PPhi_tilde_loc;
					/*DEBUG*///cout << mT << "     " << pT << "     " << cosh(P_Y_loc-p_y) << "     "
					/*DEBUG*///		<< MT_loc << "     " << PT_loc << "     " << mT*MT_loc*cosh(P_Y_loc-p_y) << "     "
					/*DEBUG*///		<< Estar_loc*Mres << "     " << (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres) << "     "
					/*DEBUG*///		<< (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc) << endl;
					VEC_Pp[is][iv][izeta] = new double [4];
					VEC_Pp[is][iv][izeta][0] = MT_loc * cosh(P_Y_loc);
					VEC_Pp[is][iv][izeta][1] = PT_loc * cos(PPhi_tilde_loc);
					VEC_Pp[is][iv][izeta][2] = PT_loc * sin(PPhi_tilde_loc);	//flip sign of this component to get VEC_Pm
					VEC_Pp[is][iv][izeta][3] = MT_loc * sinh(P_Y_loc);
					/*DEBUG*///cout << VEC_Pp[is][iv][izeta][0] << "     "
					/*DEBUG*///		<< VEC_Pp[is][iv][izeta][1] << "     "
					/*DEBUG*///		<< VEC_Pp[is][iv][izeta][2] << "     "
					/*DEBUG*///		<< VEC_Pp[is][iv][izeta][3] << endl;
				}
			}
		}
	}

	return;
}


void SourceVariances::Cal_dN_dypTdpTdphi(double** SP_p0, double** SP_px, double** SP_py, double** SP_pz, FO_surf* FOsurf_ptr)
{
   double sign = particle_sign;
   double degen = particle_gspin;
   double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);

   for(int isurf=0; isurf<FO_length ; isurf++)
   {
      FO_surf* surf = &FOsurf_ptr[isurf];
      double tau = surf->tau;
      double mu = surf->particle_mu[particle_id];
//	cerr << mu << "   " << sign << "   " << degen << "   " << endl;
      double vx = surf->vx;
      double vy = surf->vy;
      double Tdec = surf->Tdec;
      double Pdec = surf->Pdec;
      double Edec = surf->Edec;
      double da0 = surf->da0;
      double da1 = surf->da1;
      double da2 = surf->da2;
      double pi00 = surf->pi00;
      double pi01 = surf->pi01;
      double pi02 = surf->pi02;
      double pi11 = surf->pi11;
      double pi12 = surf->pi12;
      double pi22 = surf->pi22;
      double pi33 = surf->pi33;

/*
cout << "isurf: " << isurf << " of " << FO_length << endl;
cout << "tau = " << surf->tau << endl;
cout << " mu = " << surf->particle_mu[particle_id] << endl;
cout << " vx = " << surf->vx << endl;
cout << " vy = " << surf->vy << endl;
cout << " Tdec = " << surf->Tdec << endl;
cout << " Pdec = " << surf->Pdec << endl;
cout << " Edec = " << surf->Edec << endl;
cout << " da0 = " << surf->da0 << endl;
cout << " da1 = " << surf->da1 << endl;
cout << " da2 = " << surf->da2 << endl;
cout << " pi00 = " << surf->pi00 << endl;
cout << " pi01 = " << surf->pi01 << endl;
cout << " pi02 = " << surf->pi02 << endl;
cout << " pi11 = " << surf->pi11 << endl;
cout << " pi12 = " << surf->pi12 << endl;
cout << " pi22 = " << surf->pi22 << endl;
cout << " pi33 = " << surf->pi33 << endl;
cout  << endl << endl << endl;
*/

      double vT = sqrt(vx*vx + vy*vy);
      double gammaT = 1./sqrt(1. - vT*vT);

      double deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
      
      for(int ipt = 0; ipt < n_SP_pT; ipt++)
      {
      for(int iphi = 0; iphi < n_SP_pphi; iphi++)
      {
         double px = SP_px[ipt][iphi];
         double py = SP_py[ipt][iphi];
      for(int ieta=0; ieta < eta_s_npts; ieta++)
      {
         double p0 = SP_p0[ipt][ieta];
         double pz = SP_pz[ipt][ieta];
         double expon = (gammaT*(p0*1. - px*vx - py*vy) - mu)/Tdec;
         double f0;
         if(expon > 20) f0 = 0.0e0;
         else f0 = 1./(exp(expon)+sign);       //thermal equilibrium distributions

         //p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
         double pdsigma = p0*da0 + px*da1 + py*da2;

         //viscous corrections
         double Wfactor = p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
         double deltaf = 0.;
	 if (use_delta_f)
	 {
		deltaf = (1. - sign*f0)*Wfactor*deltaf_prefactor;
	 }

         double S_p = prefactor*pdsigma*f0*(1.+deltaf);
	double symmetry_factor = 1.0;
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
	 if (1. + deltaf < 0.0) S_p = 0.0;
         double S_p_withweight = S_p*tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor accounts for the assumed reflection symmetry along eta direction
//cout << "(ipt, iphi, ieta) = (" << ipt << ", " << iphi << ", " << ieta << "): " << "dN_dypTdpTdphi[ipt][iphi] = " << dN_dypTdpTdphi[ipt][iphi] << endl;
         dN_dypTdpTdphi[ipt][iphi] += S_p_withweight;
      }
//cout << "dN_dydphi[" << ipt << "][" << iphi << "] = " << dN_dypTdpTdphi[ipt][iphi] << endl;
      }
      }
   }
   return;
}

double SourceVariances::weight_function(double zvec[], int weight_function_index)
{
	double result = 0.0;
	switch(weight_function_index)
	{
		case 0:
			result = 1.;					//<1>
			break;
		case 1:
			result = zvec[2];				//<x_s>
			break;
		case 2:
			result = zvec[2]*zvec[2];			//<x^2_s>
			break;
		case 3:
			result = zvec[1];				//<x_o>
			break;
		case 4:
			result = zvec[1]*zvec[1];			//<x^2_o>
			break;
		case 5:
			result = zvec[3];				//<x_l>
			break;
		case 6:
			result = zvec[3]*zvec[3];			//<x^2_l>
			break;
		case 7:
			result = zvec[0];				//<t>
			break;
		case 8:
			result = zvec[0]*zvec[0];			//<t^2>
			break;
		case 9:
			result = zvec[2]*zvec[1];			//<x_s x_o>
			break;
		case 10:
			result = zvec[2]*zvec[3];			//<x_s x_l>
			break;
		case 11:
			result = zvec[2]*zvec[0];			//<x_s t>
			break;
		case 12:
			result = zvec[1]*zvec[3];			//<x_o x_l>
			break;
		case 13:
			result = zvec[1]*zvec[0];			//<x_o t>
			break;
		case 14:
			result = zvec[3]*zvec[0];			//<x_l t>
			break;
	}

	return result;
}


double SourceVariances::Emissionfunction(double p0, double px, double py, double pz, FO_surf* surf)
{
	//double mu = 0.0;
	//if (current_resonance_idx == 0)
	//   mu = surf->particle_mu[particle_id];
   //double sign = particle_sign;
   //double degen = particle_gspin;

   double vx = surf->vx;
   double vy = surf->vy;
   //double Tdec = surf->Tdec;
   //double Pdec = surf->Pdec;
   //double Edec = surf->Edec;
   double da0 = surf->da0;
   double da1 = surf->da1;
   double da2 = surf->da2;
   double pi00 = surf->pi00;
   double pi01 = surf->pi01;
   double pi02 = surf->pi02;
   double pi11 = surf->pi11;
   double pi12 = surf->pi12;
   double pi22 = surf->pi22;
   double pi33 = surf->pi33;

   //double vT = sqrt(vx*vx + vy*vy);
   //double gammaT = 1./sqrt(1. - vT*vT);
	double gammaT = surf->gammaT;
	

   double expon = (gammaT*(p0*1. - px*vx - py*vy) - muRES)/Tdec;
   double f0 = 1./(exp(expon)+signRES);       //thermal equilibrium distributions

   //p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
   double pdsigma = p0*da0 + px*da1 + py*da2;

   //viscous corrections
   double Wfactor = p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
   double deltaf = 0.;
   if (use_delta_f)
   {
      deltaf = (1. - signRES*f0)*Wfactor/(2.0*Tdec*Tdec*(Edec+Pdec));
   }

   //double dN_dyd2pTdphi = 1.0*gRES/(8.0*(M_PI*M_PI*M_PI))*pdsigma*f0*(1.+deltaf)/hbarC/hbarC/hbarC;
   double dN_dyd2pTdphi = gRES*pdsigma*f0*(1.+deltaf);  //put prefactor outside EmissionFunction
   if (1. + deltaf < 0.0) dN_dyd2pTdphi = 0.0;
   //out << "Spectral funct = " << dN_dyd2pTdphi << endl;

   return (dN_dyd2pTdphi);
}

void SourceVariances::Update_source_variances(int iKT, int iKphi, int reso_idx)
{
	if (reso_idx == 0)
	{
		double phi_K = K_phi[iKphi];
		double KT = K_T[iKT];
		//S_func[iKT][iKphi] = exp(interpolate2D(SPinterp2_pT, SPinterp2_pphi, ln_dN_dypTdpTdphi_moments[0][0], KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING));
		//xs_S[iKT][iKphi] = exp(interpolate2D(SPinterp2_pT, SPinterp2_pphi, ln_dN_dypTdpTdphi_moments[0][1], KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING));
		//xs2_S[iKT][iKphi] = exp(interpolate2D(SPinterp2_pT, SPinterp2_pphi, ln_dN_dypTdpTdphi_moments[0][2], KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING));
		S_func[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[0][0], KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING);
		xs_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[0][1], KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING);
		xs2_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[0][2], KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING);
		//cerr << "DEBUG: " << S_func[iKT][iKphi] << "  " << xs_S[iKT][iKphi] << "  " << xs2_S[iKT][iKphi] << endl;
	}
	else
	{
		//WRONG, JUST WANT TIMING RIGHT NOW
		//SHOULD BE LINEAR COMBINATION OF INTEGRATED_SPACETIME_MOMENTS
		S_func[iKT][iKphi] += integrated_spacetime_moments[reso_idx-1][0][iKT][iKphi];
		xs_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx-1][1][iKT][iKphi];
		xs2_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx-1][2][iKT][iKphi];
	}

return;
}

void SourceVariances::Calculate_R2_side(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xs2_S[iKT][iKphi];
   double term2 = xs_S[iKT][iKphi];

   R2_side[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   cerr << "R^2_s(KT = " << K_T[iKT] << ", Kphi = " << K_phi[iKphi] << ") = " << R2_side[iKT][iKphi] << endl;
   return;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

void SourceVariances::Output_SVdN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpTdphi;
	filename_stream_dN_dypTdpTdphi << global_path << "/SV_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpTdphi(filename_stream_dN_dypTdpTdphi.str().c_str());

	for(int iphi=0; iphi<n_SP_pphi; iphi++)
	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpTdphi << SP_pT[ipt] << "   " << SP_pphi[iphi] << "   " << dN_dypTdpTdphi[ipt][iphi] << endl;

	output_dN_dypTdpTdphi.close();

	return;
}

void SourceVariances::Output_SVdN_dypTdpT(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpT;
	filename_stream_dN_dypTdpT << global_path << "/SV_dN_dypTdpT_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpT(filename_stream_dN_dypTdpT.str().c_str());

	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpT << SP_pT[ipt] << "   " << dN_dypTdpT[ipt] << endl;

	output_dN_dypTdpT.close();

	return;
}

//End of file

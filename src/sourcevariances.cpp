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

using namespace std;

//need to define some variables for quick evaluation of S_direct
double Sdir_Y = 0.0, Sdir_R = 5., Sdir_Deltau = 1., Sdir_Deleta = 1.2;
double Sdir_eta0 = 0.0, Sdir_tau0 = 5., Sdir_etaf = 0.0, Sdir_T = 0.15;
double Sdir_prefactor, Sdir_rt, Sdir_H, Sdir_etat, Sdir_ch_Y_m_eta, Sdir_expon;
double Sdir_term1, Sdir_term2, Sdir_term3;

double S_direct(double r, double eta, double tau, double MT, double PT, double cos_phi_m_Phi)
{
	//Sdir_prefactor = (2.*Sdir_Jr + 1.)/(twopi*twopi*twopi);
	Sdir_rt = r/Sdir_R;
	Sdir_term1 = 0.5*Sdir_rt*Sdir_rt;
	Sdir_term2 = 0.5*(eta-Sdir_eta0)*(eta-Sdir_eta0)/(Sdir_Deleta*Sdir_Deleta);
	Sdir_term3 = 0.5*(tau-Sdir_tau0)*(tau-Sdir_tau0)/(Sdir_Deltau*Sdir_Deltau);
	//Sdir_term2 = 0.0;
	//Sdir_term3 = 0.0;
	Sdir_H = exp( - Sdir_term1 - Sdir_term2 - Sdir_term3 ) / (M_PI * Sdir_Deltau);
	//Sdir_H = exp( - Sdir_term1 - Sdir_term2 - Sdir_term3 );
	Sdir_ch_Y_m_eta = cosh(Sdir_Y - eta);
	Sdir_expon = -(MT/Sdir_T)*Sdir_ch_Y_m_eta*cosh(Sdir_etaf*Sdir_rt)
			+ (PT/Sdir_T)*sinh(Sdir_etaf*Sdir_rt)*cos_phi_m_Phi;
	return (MT * Sdir_ch_Y_m_eta * Sdir_H * exp(Sdir_expon));
	//return (1.);
	//return((r > Sdir_R) ? 0. : 1.);
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

void SourceVariances::Analyze_sourcefunction(FO_surf* FOsurf_ptr)
{
	time_t rawtime;
  	struct tm * timeinfo;

   //double local_prefactor = 1.0/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
//cerr << "local_prefactor * S_direct = " << local_prefactor * S_direct(1.0, 0.0, 7.5, sqrt(1.0 + 0.77*0.77), 1.0, -1.0) << endl;
//cin >> local_prefactor;
//if(1) exit;

	double plane_psi = 0.0;
	int iorder = USE_PLANE_PSI_ORDER;
	if (USE_PLANE_PSI_ORDER)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "Determine nth-order plane angles..." << endl;
		//cerr << "Determine nth-order plane angles..." << endl;
		Determine_plane_angle(current_FOsurf_ptr, 0, true);	//uses only thermal pions...
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

	//int resonance_loop_cutoff = 0;				//loop over direct pions only
	//int resonance_loop_cutoff = n_resonance;			//loop over direct pions and resonances
	int resonance_loop_cutoff = 1;					//other
	
	for (int ir = 0; ir <= resonance_loop_cutoff; ir++)
	{
		if (ir != 0 && ir != 1)
		{
			cout << "ir = " << ir << ": skipping." << endl;
			continue;
		}
		time (&rawtime);
		timeinfo = localtime (&rawtime);
		//if (VERBOSE > 0) *global_out_stream_ptr << "\t\t   * Started setting spacetime moments at " << asctime(timeinfo);
		//cerr << "Setting spacetime moments" << endl;
		Set_dN_dypTdpTdphi_moments(FOsurf_ptr, ir);
		time (&rawtime);
		timeinfo = localtime (&rawtime);
		//if (VERBOSE > 0) *global_out_stream_ptr << "\t\t   * Finished setting spacetime moments at " << asctime(timeinfo);
		Determine_plane_angle(FOsurf_ptr, ir);	//for ir==0, should give same answer as first call to Determine_plane_angle above...
		if (INTERPOLATION_FORMAT == 0) return;		//if calculating everything exactly, just do dN_d3p stuff
		//return;
		// begin HBT calculations here...
		for(int iKT = 0; iKT < n_localp_T; iKT++)
		{
			//if (iKT == 0) continue;	//temporary, need to do something slightly different for K_T == 0...
			if (abs(K_T[iKT]) < 1.e-10) continue;
			time (&rawtime);
			timeinfo = localtime (&rawtime);
			//if (VERBOSE > 0) *global_out_stream_ptr << "   - Calculating K_T = " << K_T[iKT] << " GeV at " << asctime(timeinfo);
			for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
			{
				//if (VERBOSE > 1) *global_out_stream_ptr << "\t --> Calculating K_phi = " << K_phi[iKphi] << " ..." << endl;
				//if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Loading info for this resonance..." << endl;
				Load_resonance_info(ir, K_T[iKT], K_phi[iKphi]);
				//if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Loaded resonance info" << endl;
				//ir == 0 corresponds to just thermal pions
				//otherwise, gives irth resonance contribution
				//if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Getting phase-space integrated spacetime moments" << endl;
				if (ir > 0)	//don't need to do any resonance integrals for thermal pions!
					Do_resonance_integrals(iKT, iKphi, ir);
				//if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Finished calculations for reso_idx = " << ir << endl
				//	<< "\t\t --> just updating the contributions from this resonance to each source integral now..." << endl;
				Update_source_variances(iKT, iKphi, ir);
				//if (VERBOSE > 1) *global_out_stream_ptr << "\t\t   * Moving to next K_phi loop..." << endl;
			}
		}
	}
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		//if (iKT == 0) continue;	//temporary, need to do something slightly different for K_T == 0...
		if (abs(K_T[iKT]) < 1.e-10) continue;
		double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
		beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
		for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			if (VERBOSE > 1) *global_out_stream_ptr << "\t --> Finished!  Getting R^2_ij..." << endl;
			Calculate_R2_side(iKT, iKphi);
			Calculate_R2_out(iKT, iKphi);
			Calculate_R2_long(iKT, iKphi);
			Calculate_R2_outside(iKT, iKphi);
			Calculate_R2_sidelong(iKT, iKphi);
			Calculate_R2_outlong(iKT, iKphi);
			if (VERBOSE > 1) *global_out_stream_ptr << "\t --> Moving on!" << endl;
		}
		R2_Fourier_transform(iKT, plane_psi);
	}

   return;
}

void SourceVariances::Set_dN_dypTdpTdphi_moments(FO_surf* FOsurf_ptr, int reso_idx)
{
double localmass = particle_mass;
if (reso_idx > 0) localmass = resonances.resonance_mass[reso_idx-1];
   for(int i=0; i<eta_s_npts; i++)
   {
       double local_eta_s = eta_s[i];
       double local_cosh = cosh(SP_p_y - local_eta_s);
       double local_sinh = sinh(SP_p_y - local_eta_s);
	for(int ipx=0; ipx<n_interp1_px_pts; ipx++)
	for(int ipy=0; ipy<n_interp1_py_pts; ipy++)
       {
	double mT = sqrt(localmass*localmass + SPinterp1_px[ipx]*SPinterp1_px[ipx] + SPinterp1_py[ipy]*SPinterp1_py[ipy]);
          SPinterp1_p0[ipx][ipy][i] = mT*local_cosh;
          SPinterp1_pz[ipx][ipy][i] = mT*local_sinh;
       }
	for(int ipt=0; ipt<n_interp2_pT_pts; ipt++)
       {
		double mT = sqrt(localmass*localmass + SPinterp2_pT[ipt]*SPinterp2_pT[ipt]);
		SPinterp2_p0[ipt][i] = mT*local_cosh;
		SPinterp2_pz[ipt][i] = mT*local_sinh;
		//cerr << SPinterp2_p0[ipt][i] << "    " << SPinterp2_pz[ipt][i] << "    " << SPinterp2_pT[ipt] << endl;
       }
   }
//if (1) exit(1);
//cerr << "Made it inside Set_dN_dypTdpTdphi_moments" << endl;
	if (INTERPOLATION_FORMAT == 1)	//using cartesian grid for interpolation (px, py)
		Cal_dN_dypTdpTdphi_with_weights_cartesian(FOsurf_ptr, reso_idx);
	else if (INTERPOLATION_FORMAT == 2)	//using polar grid for interpolation (pT, pphi)
		Cal_dN_dypTdpTdphi_with_weights_polar(FOsurf_ptr, reso_idx);
		//Cal_dN_dypTdpTdphi_with_weights_polar_NEW(FOsurf_ptr, reso_idx);

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

   for(int isurf=0; isurf<FO_length; isurf++)
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
	//time_t rawtime;
  	//struct tm * timeinfo;
	//time (&rawtime);
	//timeinfo = localtime (&rawtime);
	//if (VERBOSE > 0) *global_out_stream_ptr << "***Checkpoint #1 at " << asctime(timeinfo);
//cout << "***Entered Cal_dN_dypTdpTdphi_with_weights_polar" << endl;
double z0, z1, z2, z3;
double sign, degen, localmass;
if (reso_idx == 0)
{
	sign = particle_sign;
	degen = particle_gspin;
	localmass = particle_mass;
}
else
{
	sign = resonances.resonance_sign[reso_idx - 1];
	degen = resonances.resonance_gspin[reso_idx - 1];
	localmass = resonances.resonance_mass[reso_idx - 1];
}
   double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
//these are constants along the FO surface,
//so don't waste time updating them for each cell
Tdec = (&FOsurf_ptr[0])->Tdec;
double one_by_Tdec = 1./Tdec;
Pdec = (&FOsurf_ptr[0])->Pdec;
Edec = (&FOsurf_ptr[0])->Edec;
//cerr << "DEBUG: " << Tdec << "  " << Pdec << "  " << Edec << endl;
double deltaf_prefactor = 0.;
if (use_delta_f) deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

//declare variables used below here to see if this speeds up code:
double mu, tau, vx, vy, da0, da1, da2;
double pi00, pi01, pi02, pi11, pi12, pi22, pi33;
double temp_r, temp_phi, xpt, ypt, zpt, tpt, sin_temp_phi, cos_temp_phi, gammaT, expon;

double pT, pphi, sin_pphi, cos_pphi, sin_phi_m_pphi, cos_phi_m_pphi;
double px, py, p0, pz, f0, deltaf, S_p, S_p_withweight;
      FO_surf* surf;
   double eta_odd_factor = 1.0, eta_even_factor = 1.0;
   if (ASSUME_ETA_SYMMETRIC)
   {
	eta_odd_factor = 0.0;
	eta_even_factor = 2.0;
   }
	if (CHECKING_RESONANCE_CALC && USE_ANALYTIC_S)
		mu = 0.0;
	else
		mu = FOsurf_ptr[0].particle_mu[particle_id];

//cerr << "***Checkpoint #2" << endl;
   for(int isurf=0; isurf<FO_length ; isurf++)
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
	temp_phi = place_in_range(surf->phi, 0.0, 2.0*M_PI);
	xpt = surf->xpt;
	ypt = surf->ypt;
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
	pphi = SPinterp2_pphi[iphi];
//if (fabs(pphi - M_PI) < 1.e-3 && fabs(pT - 0.05) < 1.e-6)
//	cerr << "isurf = " << isurf << ": made it to pT = " << pT << " and pphi = " << pphi << endl;
	px = pT*cos_pphi;
	py = pT*sin_pphi;
	sin_phi_m_pphi = sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi;
	cos_phi_m_pphi = cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi;
	z1 = temp_r * cos_phi_m_pphi;
	z2 = temp_r * sin_phi_m_pphi;
      for(int ieta=0; ieta < eta_s_npts; ieta++)
      {
         p0 = SPinterp2_p0[ipt][ieta];
         pz = SPinterp2_pz[ipt][ieta];
	
	zpt = tau*sh_eta_s[ieta];
	tpt = tau*ch_eta_s[ieta];
	
	//now get distribution function, emission function, etc.
	if (TRUNCATE_COOPER_FRYE)
	{
		expon = (gammaT*(p0*1. - px*vx - py*vy) - mu)*one_by_Tdec;
		if (expon > 20.) continue;
		f0 = 1./(exp(expon)+sign);	//thermal equilibrium distributions
	}
	else
		f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
	//viscous corrections
	deltaf = 0.;
	if (use_delta_f)
		deltaf = deltaf_prefactor*(1. - sign*f0)*(p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);
	if (CHECKING_RESONANCE_CALC && USE_ANALYTIC_S)
	{
		S_p = prefactor * S_direct(temp_r, eta_s[ieta], tau, sqrt(pT*pT + localmass*localmass), pT, cos(temp_phi - 0.00342436));
		//S_p = S_direct(temp_r, eta_s[ieta], tau, sqrt(pT*pT + localmass*localmass), pT, cos_phi_m_pphi);
		/*if (reso_idx == 1 && fabs(pphi - M_PI) < 1.e-6 && fabs(pT - 0.05) < 1.e-6)
		{
			cerr << reso_idx << "   " << setw(8) << setprecision(10) << temp_r << "   " << temp_phi << "   " << eta_s[ieta] << "   " << tau
				<< "   " << sqrt(pT*pT + localmass*localmass) << "   " << pT << "   " << cos(temp_phi - 0.00342436)
				<< "   " << S_p << endl;
		}*/
		S_p = prefactor * S_direct(temp_r, eta_s[ieta], tau, sqrt(pT*pT + localmass*localmass), pT, cos_phi_m_pphi);
		//if ((1. + deltaf < 0.0) || ((p0*da0 + px*da1 + py*da2) < 0.0) || (flagneg == 1 && S_p < tol))
		//	S_p = 0.0;
	}
	else
	{
		//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
		S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);
		//ignore points where delta f is large or emission function goes negative from pdsigma
		if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol)) S_p = 0.0;
	}
        if (flagneg == 1 && S_p < tol)
        {
           S_p = 0.0e0;
        }

	S_p_withweight = S_p*tau*eta_s_weight[ieta]; //eta_even_factor accounts for the assumed reflection symmetry along eta direction
	z0 = tau*ch_eta_s[ieta];
	z3 = tau*sh_eta_s[ieta];
	dN_dypTdpTdphi_moments[reso_idx][0][ipt][iphi] += eta_even_factor*S_p_withweight;			//<1>
	dN_dypTdpTdphi_moments[reso_idx][1][ipt][iphi] += eta_even_factor*S_p_withweight*z2;			//<x_s>
	dN_dypTdpTdphi_moments[reso_idx][2][ipt][iphi] += eta_even_factor*S_p_withweight*z2*z2;			//<x^2_s>
	dN_dypTdpTdphi_moments[reso_idx][3][ipt][iphi] += eta_even_factor*S_p_withweight*z1;			//<x_o>
	dN_dypTdpTdphi_moments[reso_idx][4][ipt][iphi] += eta_even_factor*S_p_withweight*z1*z1;			//<x^2_o>
	dN_dypTdpTdphi_moments[reso_idx][5][ipt][iphi] += eta_odd_factor*S_p_withweight*z3;			//<x_l>
	dN_dypTdpTdphi_moments[reso_idx][6][ipt][iphi] += eta_even_factor*S_p_withweight*z3*z3;			//<x^2_l>
	dN_dypTdpTdphi_moments[reso_idx][7][ipt][iphi] += eta_even_factor*S_p_withweight*z0;			//<t>
	dN_dypTdpTdphi_moments[reso_idx][8][ipt][iphi] += eta_even_factor*S_p_withweight*z0*z0;			//<t^2>
	dN_dypTdpTdphi_moments[reso_idx][9][ipt][iphi] += eta_even_factor*S_p_withweight*z2*z1;			//<x_s x_o>
	dN_dypTdpTdphi_moments[reso_idx][10][ipt][iphi] += eta_odd_factor*S_p_withweight*z2*z3;			//<x_s x_l>
	dN_dypTdpTdphi_moments[reso_idx][11][ipt][iphi] += eta_even_factor*S_p_withweight*z2*z0;		//<x_s t>
	dN_dypTdpTdphi_moments[reso_idx][12][ipt][iphi] += eta_odd_factor*S_p_withweight*z1*z3;			//<x_o x_l>
	dN_dypTdpTdphi_moments[reso_idx][13][ipt][iphi] += eta_even_factor*S_p_withweight*z1*z0;		//<x_o t>
	dN_dypTdpTdphi_moments[reso_idx][14][ipt][iphi] += eta_odd_factor*S_p_withweight*z3*z0;			//<x_l t>
      }
      }
      }
   }
double temp;
//set log of dN_dypTdpTdphi_moments...
	for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
	for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
	for(int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		temp = dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi];
		ln_dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi] = log(abs(temp));
		sign_of_dN_dypTdpTdphi_moments[reso_idx][wfi][ipt][iphi] = sgn(temp);
		if (wfi == 0 || wfi == 5 || wfi == 6)
			cout << "-1   " << reso_idx << "   " << setw(8) << setprecision(15) << wfi << "   "
				<< SPinterp2_pT[ipt] << "   " << SPinterp2_pphi[iphi] << "   " << temp << endl;
	}
//cerr << "***Checkpoint #3" << endl;
//cout << "***Exited Cal_dN_dypTdpTdphi_with_weights_polar" << endl;
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
		//dN_dypTdpTdphi[ipt][iphi] = exp(interpolate2D(SPinterp1_px, SPinterp1_py, dN_dypTdpTdphi_moments[0][0],
		//					px, py, n_interp1_px_pts, n_interp1_py_pts, INTERPOLATION_KIND, UNIFORM_SPACING));
		dN_dypTdpTdphi[ipt][iphi] = interpolate2D(SPinterp1_px, SPinterp1_py, dN_dypTdpTdphi_moments[0][0],
							px, py, n_interp1_px_pts, n_interp1_py_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
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
		//dN_dypTdpTdphi[ipt][iphi] = exp(interpolate2D(SPinterp2_pT, SPinterp2_pphi, ln_dN_dypTdpTdphi_moments[0][0],
		//					pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING));
		dN_dypTdpTdphi[ipt][iphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[0][0],
							pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		cout << "-2   " << pT << "  " << pphi << "  " << pT*cos(pphi) << "  " << pT*sin(pphi) << "  " << dN_dypTdpTdphi[ipt][iphi] << endl;
	}
	//cout << endl;
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



void SourceVariances::Determine_plane_angle(FO_surf* FOsurf_ptr, int reso_idx, bool thermal_particles_only /*= false*/)
{
   double localmass = particle_mass;
   if (reso_idx > 0) return;	//don't really care about thermal resonance distributions at the moment,
				//so just skip this part
   //if (reso_idx > 0) localmass = resonances.resonance_mass[reso_idx-1];
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
      mT[ipt] = sqrt(localmass*localmass + SP_pT[ipt]*SP_pT[ipt]);
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

	if (INTERPOLATION_FORMAT == 0 || thermal_particles_only)	//no interpolation at all; calculate everything exactly
	{
		Cal_dN_dypTdpTdphi(p0, px, py, pz, FOsurf_ptr);
		if (VERBOSE > 0) *global_out_stream_ptr << "  --> No interpolation!  Calculating everything exactly..." << endl;
	}
	else if (INTERPOLATION_FORMAT == 1)				//using cartesian grid for interpolation (px, py)
	{
		Cal_dN_dypTdpTdphi_interpolate_cartesian_grid(px, py);
		if (VERBOSE > 0) *global_out_stream_ptr << "  --> Interpolating on a Cartesian grid..." << endl;
	}
	else if (INTERPOLATION_FORMAT == 2)				//using polar grid for interpolation (pT, pphi)
	{
   		Cal_dN_dypTdpTdphi_interpolate_polar_grid(SP_pT, SP_pphi);
		if (VERBOSE > 0) *global_out_stream_ptr << "  --> Interpolating on a polar grid..." << endl;
	}
	else
	{
		if (VERBOSE > 0) cerr << "Determine_plane_angle: Options not supported!  Exiting..." << endl;
		exit(1);
	}

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

   for(int ipt=0; ipt<n_SP_pT; ipt++)
	dN_dypTdpT[ipt] /= (2.*M_PI);
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

void SourceVariances::Load_resonance_info(int reso_idx, double K_T_local, double K_phi_local)
{
	current_resonance_idx = reso_idx;
	//cerr << "Entering Load_resonance_info for reso_idx = " << reso_idx << endl;
	if (reso_idx == 0)
	{
		muRES = particle_mu;
		signRES = particle_sign;
		gRES = particle_gspin;
		
		/*if (DEBUG)
		{
			cerr << "Working on resonance # = " << current_resonance_idx << ":" << endl
				<< "  --> muRES = " << muRES << endl
				<< "  --> signRES = " << signRES << endl
				<< "  --> gRES = " << gRES << endl;
		}*/
		
		return;
	}
	else
	{
		//cerr << "Entering loop in Load_resonance_info for reso_idx = " << reso_idx << endl;
		current_resonance_mass = resonances.resonance_mass[reso_idx-1];
		current_resonance_Gamma = resonances.resonance_Gamma[reso_idx-1];
		current_resonance_total_br = resonances.resonance_total_br[reso_idx-1];
		current_resonance_decay_masses[0] = resonances.resonance_decay_masses[reso_idx-1][0];
		current_resonance_decay_masses[1] = resonances.resonance_decay_masses[reso_idx-1][1];

		muRES = resonances.resonance_mu[reso_idx-1];
		signRES = resonances.resonance_sign[reso_idx-1];
		gRES = resonances.resonance_gspin[reso_idx-1];

		Mres = current_resonance_mass;
		mass = particle_mass;
		Gamma = current_resonance_Gamma;
		double one_by_Gamma_Mres = hbarC/(Gamma*Mres);
		br = current_resonance_total_br;
		m2 = current_resonance_decay_masses[0];
		m3 = current_resonance_decay_masses[1];
		if (DEBUG)
		{
			cerr << "Working on resonance # = " << current_resonance_idx << ":" << endl
				<< "  --> muRES = " << muRES << endl
				<< "  --> signRES = " << signRES << endl
				<< "  --> gRES = " << gRES << endl
				<< "  --> Mres = " << Mres << endl
				<< "  --> mass = " << mass << endl
				<< "  --> Gamma = " << Gamma << endl
				<< "  --> br = " << br << endl
				<< "  --> m2 = " << m2 << endl
				<< "  --> m3 = " << m3 << endl << endl;
		}
		if (abs(m3) <= 1.e-6)
			n_body = 2;
		else
			n_body = 3;
		mT = sqrt(mass*mass + K_T_local*K_T_local);
		pT = K_T_local;
		current_K_phi = K_phi_local;
		cos_cKphi = cos(K_phi_local);
		sin_cKphi = sin(K_phi_local);

		if (n_body == 2)
		{
			//set up vectors of points to speed-up integrals...
			double s_loc = m2*m2;
			double pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc) )/(2.0*Mres);
			VEC_n2_pstar = pstar_loc;
			double g_s_loc = g(s_loc)/pstar_loc;	//for n_body == 2, doesn't actually use s_loc since result is just a factor * delta(...); just returns factor
			VEC_n2_g_s = g_s_loc;
			//VEC_n2_s_factor = br/(4.*M_PI*VEC_n2_pstar);
			VEC_n2_s_factor = VEC_n2_g_s;
			double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
			VEC_n2_Estar = Estar_loc;
			double psBmT = pstar_loc / mT;
			double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
			VEC_n2_DeltaY = DeltaY_loc;
			p_y = 0.0;
			VEC_n2_v_factor = new double [n_v_pts];
			VEC_n2_zeta_factor = new double * [n_v_pts];
			VEC_n2_Yp = p_y + DeltaY_loc;
			VEC_n2_Ym = p_y - DeltaY_loc;
			VEC_n2_P_Y = new double [n_v_pts];
			VEC_n2_MTbar = new double [n_v_pts];
			VEC_n2_DeltaMT = new double [n_v_pts];
			VEC_n2_MTp = new double [n_v_pts];
			VEC_n2_MTm = new double [n_v_pts];
			VEC_n2_MT = new double * [n_v_pts];
			VEC_n2_PPhi_tilde = new double * [n_v_pts];
			VEC_n2_PPhi_tildeFLIP = new double * [n_v_pts];
			VEC_n2_PT = new double * [n_v_pts];
			VEC_n2_Pp = new double ** [n_v_pts];
			VEC_n2_alpha = new double ** [n_v_pts];
			VEC_n2_Pm = new double ** [n_v_pts];
			VEC_n2_alpha_m = new double ** [n_v_pts];
			for(int iv = 0; iv < n_v_pts; iv++)
			{
				//cerr << "In v loop# = " << iv << endl;
				double v_loc = v_pts[iv];
				double P_Y_loc = p_y + v_loc*DeltaY_loc;
				VEC_n2_P_Y[iv] = P_Y_loc;
				double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
				double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
				VEC_n2_v_factor[iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);
				double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
				VEC_n2_MTbar[iv] = MTbar_loc;
				double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
				VEC_n2_DeltaMT[iv] = DeltaMT_loc;
				VEC_n2_MTp[iv] = MTbar_loc + DeltaMT_loc;
				VEC_n2_MTm[iv] = MTbar_loc - DeltaMT_loc;
				VEC_n2_MT[iv] = new double [n_zeta_pts];
				VEC_n2_PPhi_tilde[iv] = new double [n_zeta_pts];
				VEC_n2_PPhi_tildeFLIP[iv] = new double [n_zeta_pts];
				VEC_n2_PT[iv] = new double [n_zeta_pts];
				VEC_n2_Pp[iv] = new double * [n_zeta_pts];
				VEC_n2_alpha[iv] = new double * [n_zeta_pts];
				VEC_n2_Pm[iv] = new double * [n_zeta_pts];
				VEC_n2_alpha_m[iv] = new double * [n_zeta_pts];
				VEC_n2_zeta_factor[iv] = new double [n_zeta_pts];
				for(int izeta = 0; izeta < n_zeta_pts; izeta++)
				{
					double zeta_loc = zeta_pts[izeta];
					double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
					VEC_n2_MT[iv][izeta] = MT_loc;
					VEC_n2_zeta_factor[iv][izeta] = zeta_wts[izeta]*MT_loc;
					double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
					double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
					//assume that PPhi_tilde is +ve in next step...
					double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
					double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), interp2_pphi_min, interp2_pphi_max);
					VEC_n2_PPhi_tilde[iv][izeta] = place_in_range( K_phi_local + PPhi_tilde_loc, interp2_pphi_min, interp2_pphi_max);
					VEC_n2_PPhi_tildeFLIP[iv][izeta] = place_in_range( K_phi_local - PPhi_tilde_loc, interp2_pphi_min, interp2_pphi_max);
					VEC_n2_PT[iv][izeta] = PT_loc;
					/*DEBUG*///cout << mT << "     " << pT << "     " << cosh(P_Y_loc-p_y) << "     "
					/*DEBUG*///		<< MT_loc << "     " << PT_loc << "     " << mT*MT_loc*cosh(P_Y_loc-p_y) << "     "
					/*DEBUG*///		<< Estar_loc*Mres << "     " << (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres) << "     "
					/*DEBUG*///		<< (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc) << endl;
					VEC_n2_Pp[iv][izeta] = new double [4];
					VEC_n2_alpha[iv][izeta] = new double [4];
					VEC_n2_Pm[iv][izeta] = new double [4];
					VEC_n2_alpha_m[iv][izeta] = new double [4];
					//probably not the most elegant set-up, but does the job for now...
					VEC_n2_Pp[iv][izeta][0] = MT_loc * cosh(P_Y_loc);
					VEC_n2_Pp[iv][izeta][1] = PT_loc * cos(K_phi_local + PPhi_tilde_loc);
					VEC_n2_Pp[iv][izeta][2] = PT_loc * sin(K_phi_local + PPhi_tilde_loc);
					VEC_n2_Pp[iv][izeta][3] = MT_loc * sinh(P_Y_loc);
					VEC_n2_Pm[iv][izeta][0] = VEC_n2_Pp[iv][izeta][0];
					VEC_n2_Pm[iv][izeta][1] = PT_loc * cos(K_phi_local - PPhi_tilde_loc);
					VEC_n2_Pm[iv][izeta][2] = PT_loc * sin(K_phi_local - PPhi_tilde_loc);
					VEC_n2_Pm[iv][izeta][3] = VEC_n2_Pp[iv][izeta][3];
					for (int ii=0; ii<4; ii++)
					{
						VEC_n2_alpha[iv][izeta][ii] = one_by_Gamma_Mres * VEC_n2_Pp[iv][izeta][ii];
						VEC_n2_alpha_m[iv][izeta][ii] = one_by_Gamma_Mres * VEC_n2_Pm[iv][izeta][ii];
					}
					/*DEBUG*///cout << VEC_n2_Pp[iv][izeta][0] << "     "
					/*DEBUG*///		<< VEC_n2_Pp[iv][izeta][1] << "     "
					/*DEBUG*///		<< VEC_n2_Pp[iv][izeta][2] << "     "
					/*DEBUG*///		<< VEC_n2_Pp[iv][izeta][3] << endl;
				}
			}
		}
		else
		{
			Qfunc = get_Q();
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
			VEC_PPhi_tildeFLIP = new double ** [n_s_pts];
			VEC_PT = new double ** [n_s_pts];
			VEC_Pp = new double *** [n_s_pts];
			VEC_alpha = new double *** [n_s_pts];
			VEC_Pm = new double *** [n_s_pts];
			VEC_alpha_m = new double *** [n_s_pts];
			//cerr << "Entering loops now..." << endl;
			for (int is = 0; is < n_s_pts; is++)
			{
				//cerr << "In s loop# = " << is << endl;
				double s_loc = s_pts[reso_idx-1][is];
				double g_s_loc = g(s_loc);
				VEC_g_s[is] = g_s_loc;
				VEC_s_factor[is] = s_wts[reso_idx-1][is]*g_s_loc;
				double pstar_loc = sqrt(((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc))/(2.0*Mres);
				VEC_pstar[is] = pstar_loc;
				double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
				VEC_Estar[is] = Estar_loc;
				double psBmT = pstar_loc / mT;
				double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
				VEC_DeltaY[is] = DeltaY_loc;
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
				VEC_PPhi_tildeFLIP[is] = new double * [n_v_pts];
				VEC_PT[is] = new double * [n_v_pts];
				VEC_Pp[is] = new double ** [n_v_pts];
				VEC_alpha[is] = new double ** [n_v_pts];
				VEC_Pm[is] = new double ** [n_v_pts];
				VEC_alpha_m[is] = new double ** [n_v_pts];
				for(int iv = 0; iv < n_v_pts; iv++)
				{
					//cerr << "In v loop# = " << iv << endl;
					double v_loc = v_pts[iv];
					double P_Y_loc = p_y + v_loc*DeltaY_loc;
					VEC_P_Y[is][iv] = P_Y_loc;
					double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
					double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
					VEC_v_factor[is][iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);
					double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
					VEC_MTbar[is][iv] = MTbar_loc;
					double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
					VEC_DeltaMT[is][iv] = DeltaMT_loc;
					VEC_MTp[is][iv] = MTbar_loc + DeltaMT_loc;
					VEC_MTm[is][iv] = MTbar_loc - DeltaMT_loc;
					VEC_MT[is][iv] = new double [n_zeta_pts];
					VEC_PPhi_tilde[is][iv] = new double [n_zeta_pts];
					VEC_PPhi_tildeFLIP[is][iv] = new double [n_zeta_pts];
					VEC_PT[is][iv] = new double [n_zeta_pts];
					VEC_Pp[is][iv] = new double * [n_zeta_pts];
					VEC_alpha[is][iv] = new double * [n_zeta_pts];
					VEC_Pm[is][iv] = new double * [n_zeta_pts];
					VEC_alpha_m[is][iv] = new double * [n_zeta_pts];
					VEC_zeta_factor[is][iv] = new double [n_zeta_pts];
					for(int izeta = 0; izeta < n_zeta_pts; izeta++)
					{
						double zeta_loc = zeta_pts[izeta];
						double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
						VEC_MT[is][iv][izeta] = MT_loc;
						VEC_zeta_factor[is][iv][izeta] = zeta_wts[izeta]*MT_loc;
						double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
						double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
						//assume that PPhi_tilde is +ve in next step...
						double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
						double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), interp2_pphi_min, interp2_pphi_max);
						VEC_PPhi_tilde[is][iv][izeta] = place_in_range( K_phi_local + PPhi_tilde_loc, interp2_pphi_min, interp2_pphi_max);
						VEC_PPhi_tildeFLIP[is][iv][izeta] = place_in_range( K_phi_local - PPhi_tilde_loc, interp2_pphi_min, interp2_pphi_max);
						VEC_PT[is][iv][izeta] = PT_loc;
						/*DEBUG*///cout << mT << "     " << pT << "     " << cosh(P_Y_loc-p_y) << "     "
						/*DEBUG*///		<< MT_loc << "     " << PT_loc << "     " << mT*MT_loc*cosh(P_Y_loc-p_y) << "     "
						/*DEBUG*///		<< Estar_loc*Mres << "     " << (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres) << "     "
						/*DEBUG*///		<< (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc) << endl;
						VEC_Pp[is][iv][izeta] = new double [4];
						VEC_alpha[is][iv][izeta] = new double [4];
						VEC_Pm[is][iv][izeta] = new double [4];
						VEC_alpha_m[is][iv][izeta] = new double [4];
						VEC_Pp[is][iv][izeta][0] = MT_loc * cosh(P_Y_loc);
						VEC_Pp[is][iv][izeta][1] = PT_loc * cos(K_phi_local + PPhi_tilde_loc);
						VEC_Pp[is][iv][izeta][2] = PT_loc * sin(K_phi_local + PPhi_tilde_loc);
						VEC_Pp[is][iv][izeta][3] = MT_loc * sinh(P_Y_loc);
						VEC_Pm[is][iv][izeta][0] = VEC_Pp[is][iv][izeta][0];
						VEC_Pm[is][iv][izeta][1] = PT_loc * cos(K_phi_local - PPhi_tilde_loc);
						VEC_Pm[is][iv][izeta][2] = PT_loc * sin(K_phi_local - PPhi_tilde_loc);
						VEC_Pm[is][iv][izeta][3] = VEC_Pp[is][iv][izeta][3];
						for (int ii=0; ii<4; ii++)
						{
							VEC_alpha[is][iv][izeta][ii] = one_by_Gamma_Mres * VEC_Pp[is][iv][izeta][ii];
							VEC_alpha_m[is][iv][izeta][ii] = one_by_Gamma_Mres * VEC_Pm[is][iv][izeta][ii];
						}
						/*DEBUG*///cout << VEC_Pp[is][iv][izeta][0] << "     "
						/*DEBUG*///		<< VEC_Pp[is][iv][izeta][1] << "     "
						/*DEBUG*///		<< VEC_Pp[is][iv][izeta][2] << "     "
						/*DEBUG*///		<< VEC_Pp[is][iv][izeta][3] << endl;
					}
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
	double localmass = particle_mass;

   for(int isurf=0; isurf<FO_length ; isurf++)
   {
      FO_surf* surf = &FOsurf_ptr[isurf];
      double tau = surf->tau;
	double mu;
	if (CHECKING_RESONANCE_CALC && USE_ANALYTIC_S)
		mu = 0.0;
	else
		mu = FOsurf_ptr[0].particle_mu[particle_id];
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
      double vT = sqrt(vx*vx + vy*vy);
      double gammaT = 1./sqrt(1. - vT*vT);
	double temp_r = surf->r;
	double temp_phi = surf->phi;
	//double sin_temp_phi = surf->sin_phi;
	//double cos_temp_phi = surf->cos_phi;

      double deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
      
      for(int ipt = 0; ipt < n_SP_pT; ipt++)
      {
	double pT = SP_pT[ipt];
      for(int iphi = 0; iphi < n_SP_pphi; iphi++)
      {
         double px = SP_px[ipt][iphi];
         double py = SP_py[ipt][iphi];
	double cos_phi_m_pphi = cos(temp_phi - SP_pphi[iphi]);
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

	double S_p;

	if (CHECKING_RESONANCE_CALC && USE_ANALYTIC_S)
	{
		S_p = prefactor * S_direct(temp_r, eta_s[ieta], tau, sqrt(pT*pT + localmass*localmass), pT, cos_phi_m_pphi);
		//S_p = S_direct(temp_r, eta_s[ieta], tau, sqrt(pT*pT + localmass*localmass), pT, cos_phi_m_pphi);
		//if (VERBOSE > 0 && isurf == 0) *global_out_stream_ptr << "  --> Cal_dN_dypTdpTdphi(): computed emission function with S_direct()..." << endl;
		//if (isurf == 0)
		//{
		//	cerr << SP_pphi[iphi] << "   " << SP_pT[ipt] << endl;
		//}
		//if (/*isurf == 0 && */fabs(SP_pphi[iphi] - 0.00386099) < 1.e-3 && fabs(pT - 0.0180112) < 1.e-6)
		//if (isurf == 0)
		/*{
			cerr << temp_r << "   " << temp_phi << "   " << eta_s[ieta] << "   " << tau
				<< "   " << sqrt(pT*pT + localmass*localmass) << "   " << pT << "   " << cos_phi_m_pphi
				<< "   " << S_p*tau << endl;
		}*/
		//if ((1. + deltaf < 0.0) || ((p0*da0 + px*da1 + py*da2) < 0.0) || (flagneg == 1 && S_p < tol))
		//	S_p = 0.0;
	}
	else
	{
		//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
		S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);
		//ignore points where delta f is large or emission function goes negative from pdsigma
		//if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol)) S_p = 0.0;
        if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
        //if (1. + deltaf < 0.0)
		S_p = 0.0e0;
	}


         //double S_p = prefactor*pdsigma*f0*(1.+deltaf);
	double symmetry_factor = 1.0;
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
	 //if (1. + deltaf < 0.0) S_p = 0.0;
         double S_p_withweight = S_p*tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor accounts for the assumed reflection symmetry along eta direction
//cout << "(ipt, iphi, ieta) = (" << ipt << ", " << iphi << ", " << ieta << "): " << "dN_dypTdpTdphi[ipt][iphi] = " << dN_dypTdpTdphi[ipt][iphi] << endl;
//JUST COMPARING ANALYTIC RESONANCES WITH SV.S CALCULATIONS: WRONG OTHERWISE
	//S_p_withweight = (1.)*eta_s_weight[ieta]*symmetry_factor;
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
	switch(weight_function_index)
	{
		case 0:
			return (1.);					//<1>
		case 1:
			return (zvec[2]);				//<x_s>
		case 2:
			return (zvec[2]*zvec[2]);			//<x^2_s>
		case 3:
			return (zvec[1]);				//<x_o>
		case 4:
			return (zvec[1]*zvec[1]);			//<x^2_o>
		case 5:
			return (zvec[3]);				//<x_l>
		case 6:
			return (zvec[3]*zvec[3]);			//<x^2_l>
		case 7:
			return (zvec[0]);				//<t>
		case 8:
			return (zvec[0]*zvec[0]);			//<t^2>
		case 9:
			return (zvec[2]*zvec[1]);			//<x_s x_o>
		case 10:
			return (zvec[2]*zvec[3]);			//<x_s x_l>
		case 11:
			return (zvec[2]*zvec[0]);			//<x_s t>
		case 12:
			return (zvec[1]*zvec[3]);			//<x_o x_l>
		case 13:
			return (zvec[1]*zvec[0]);			//<x_o t>
		case 14:
			return (zvec[3]*zvec[0]);			//<x_l t>
	}
}

void SourceVariances::Update_source_variances(int iKT, int iKphi, int reso_idx)
{
	if (VERBOSE > 1) *global_out_stream_ptr << "\t\t\t   () Successfully entered Update_source_variances()..." << endl;
	//if (reso_idx == 0 || reso_idx == 1)
	if (reso_idx == 0)
	{
		double phi_K = K_phi[iKphi];
		double KT = K_T[iKT];
		double temp_factor = 0.0;
		S_func[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][0],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xs_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][1],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xs2_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][2],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xo_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][3],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xo2_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][4],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xl_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][5],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xl2_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][6],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);		
		t_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][7],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		t2_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][8],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);		
		xo_xs_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][9],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xl_xs_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][10],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xs_t_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][11],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xo_xl_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][12],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xo_t_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][13],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xl_t_S[iKT][iKphi] = temp_factor * interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][14],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		if (DEBUG)		
			cerr << "DEBUG: reso_idx = " << reso_idx << "  " << S_func[iKT][iKphi] << "  " << xs_S[iKT][iKphi] << "  " << xs2_S[iKT][iKphi] << endl;
	}
	else
	{
		S_func[iKT][iKphi] += integrated_spacetime_moments[reso_idx][0][iKT][iKphi];
		//cerr << "reso_idx = 1: " << integrated_spacetime_moments[reso_idx-1][0][iKT][iKphi] << endl;
		xs_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][1][iKT][iKphi];
		xs2_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][2][iKT][iKphi];
		xo_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][3][iKT][iKphi];
		xo2_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][4][iKT][iKphi];
		xl_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][5][iKT][iKphi];
		xl2_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][6][iKT][iKphi];
		t_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][7][iKT][iKphi];
		t2_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][8][iKT][iKphi];
		xo_xs_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][9][iKT][iKphi];
		xl_xs_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][10][iKT][iKphi];
		xs_t_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][11][iKT][iKphi];
		xo_xl_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][12][iKT][iKphi];
		xo_t_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][13][iKT][iKphi];
		xl_t_S[iKT][iKphi] += integrated_spacetime_moments[reso_idx][14][iKT][iKphi];
		if (DEBUG)		
			cerr << "DEBUG: reso_idx = " << reso_idx << "  " << integrated_spacetime_moments[reso_idx][0][iKT][iKphi]
				<< "  " << integrated_spacetime_moments[reso_idx][1][iKT][iKphi]
				<< "  " << integrated_spacetime_moments[reso_idx][2][iKT][iKphi] << endl;
	}

	if (VERBOSE > 1) *global_out_stream_ptr << "\t\t\t   () Successfully completed Update_source_variances()..." << endl;
	return;
}

void SourceVariances::Calculate_R2_side(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xs2_S[iKT][iKphi];
   double term2 = xs_S[iKT][iKphi];

   R2_side[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   if (DEBUG)
	cerr << "DEBUG: R^2_s(KT = " << K_T[iKT] << ", Kphi = " << K_phi[iKphi] << ") = " << R2_side[iKT][iKphi] << endl;
   return;
}

void SourceVariances::Calculate_R2_out(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xo2_S[iKT][iKphi] - 2.*beta_perp*xo_t_S[iKT][iKphi] + beta_perp*beta_perp*t2_S[iKT][iKphi];
   double term2 = xo_S[iKT][iKphi] - beta_perp*t_S[iKT][iKphi];

   R2_out[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   if (DEBUG)
	cerr << "DEBUG: R^2_o(KT = " << K_T[iKT] << ", Kphi = " << K_phi[iKphi] << ") = " << R2_out[iKT][iKphi] << endl;
   return;
}

void SourceVariances::Calculate_R2_outside(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xo_xs_S[iKT][iKphi] - beta_perp*xs_t_S[iKT][iKphi];
   double term2 = xo_S[iKT][iKphi] - beta_perp*t_S[iKT][iKphi];
   double term3 = xs_S[iKT][iKphi];

   R2_outside[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
   if (DEBUG)
	cerr << "DEBUG: R^2_os(KT = " << K_T[iKT] << ", Kphi = " << K_phi[iKphi] << ") = " << R2_outside[iKT][iKphi] << endl;
   return;
}

void SourceVariances::Calculate_R2_long(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xl2_S[iKT][iKphi] - 2.*beta_l*xl_t_S[iKT][iKphi] + beta_l*beta_l*t2_S[iKT][iKphi];
   double term2 = xl_S[iKT][iKphi] - beta_l*t_S[iKT][iKphi];

   R2_long[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   if (DEBUG)
	cerr << "DEBUG: R^2_l(KT = " << K_T[iKT] << ", Kphi = " << K_phi[iKphi] << ") = " << R2_long[iKT][iKphi] << endl;
   return;
}

void SourceVariances::Calculate_R2_outlong(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xo_xl_S[iKT][iKphi] - beta_perp*xl_t_S[iKT][iKphi] - beta_l*xo_t_S[iKT][iKphi] + beta_perp*beta_l*t2_S[iKT][iKphi];
   double term2 = xo_S[iKT][iKphi] - beta_perp*t_S[iKT][iKphi];
   double term3 = xl_S[iKT][iKphi] - beta_l*t_S[iKT][iKphi];

   R2_outlong[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
   if (DEBUG)
	cerr << "DEBUG: R^2_ol(KT = " << K_T[iKT] << ", Kphi = " << K_phi[iKphi] << ") = " << R2_outlong[iKT][iKphi] << endl;
   return;
}

void SourceVariances::Calculate_R2_sidelong(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xl_xs_S[iKT][iKphi] - beta_l*xs_t_S[iKT][iKphi];
   double term2 = xs_S[iKT][iKphi];
   double term3 = xl_S[iKT][iKphi] - beta_l*t_S[iKT][iKphi];

   R2_sidelong[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
   if (DEBUG)
	cerr << "DEBUG: R^2_sl(KT = " << K_T[iKT] << ", Kphi = " << K_phi[iKphi] << ") = " << R2_sidelong[iKT][iKphi] << endl;
   return;
}

void SourceVariances::R2_Fourier_transform(int iKT, double plane_psi)
{
	for(int Morder=0; Morder<n_order; Morder++)
	{
		double cos_mK_phi[n_localp_phi], sin_mK_phi[n_localp_phi];
		for(int i=0; i<n_localp_phi; i++)
		{
			cos_mK_phi[i] = cos(Morder*(K_phi[i] - plane_psi));
			sin_mK_phi[i] = sin(Morder*(K_phi[i] - plane_psi));
		}
		double temp_sum_side_cos = 0.0e0;
		double temp_sum_side_sin = 0.0e0;
		double temp_sum_out_cos = 0.0e0;
		double temp_sum_out_sin = 0.0e0;
		double temp_sum_outside_cos = 0.0e0;
		double temp_sum_outside_sin = 0.0e0;
		double temp_sum_long_cos = 0.0e0;
		double temp_sum_long_sin = 0.0e0;
		double temp_sum_sidelong_cos = 0.0e0;
		double temp_sum_sidelong_sin = 0.0e0;
		double temp_sum_outlong_cos = 0.0e0;
		double temp_sum_outlong_sin = 0.0e0;
		for(int i=0; i<n_localp_phi; i++)
		{
			temp_sum_side_cos += R2_side[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_side_sin += R2_side[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_out_cos += R2_out[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_out_sin += R2_out[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_outside_cos += R2_outside[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_outside_sin += R2_outside[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_long_cos += R2_long[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_long_sin += R2_long[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_sidelong_cos += R2_sidelong[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_sidelong_sin += R2_sidelong[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
			temp_sum_outlong_cos += R2_outlong[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
			temp_sum_outlong_sin += R2_outlong[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
		}
		R2_side_C[iKT][Morder] = temp_sum_side_cos/(2*M_PI);
		R2_side_S[iKT][Morder] = temp_sum_side_sin/(2*M_PI);
		R2_out_C[iKT][Morder] = temp_sum_out_cos/(2*M_PI);
		R2_out_S[iKT][Morder] = temp_sum_out_sin/(2*M_PI);
		R2_outside_C[iKT][Morder] = temp_sum_outside_cos/(2*M_PI);
		R2_outside_S[iKT][Morder] = temp_sum_outside_sin/(2*M_PI);
		R2_long_C[iKT][Morder] = temp_sum_long_cos/(2*M_PI);
		R2_long_S[iKT][Morder] = temp_sum_long_sin/(2*M_PI);
		R2_sidelong_C[iKT][Morder] = temp_sum_sidelong_cos/(2*M_PI);
		R2_sidelong_S[iKT][Morder] = temp_sum_sidelong_sin/(2*M_PI);
		R2_outlong_C[iKT][Morder] = temp_sum_outlong_cos/(2*M_PI);
		R2_outlong_S[iKT][Morder] = temp_sum_outlong_sin/(2*M_PI);
	}
	return;
}

//**********************************************************************************************
bool SourceVariances::particles_are_the_same(int idx1, int idx2)
{
    if (all_particles[idx1].sign!=all_particles[idx2].sign) return false;
    if (abs(all_particles[idx1].mass-all_particles[idx2].mass) / (all_particles[idx2].mass+1e-30) > PARTICLE_DIFF_TOLERANCE) return false;
    for (long l=0; l<FO_length; l++)
    {
        double chem1 = current_FOsurf_ptr[l].particle_mu[idx1], chem2 = current_FOsurf_ptr[l].particle_mu[idx2];
        if (abs(chem1-chem2)/(chem2+1e-30) > PARTICLE_DIFF_TOLERANCE)
        {
          return false;
        }

    }

    return true;
}

//End of file

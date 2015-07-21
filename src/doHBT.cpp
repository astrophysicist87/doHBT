#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<algorithm>
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

using namespace std;

//need to define some variables for quick evaluation of S_direct
double Sdir_Y = 0.0, Sdir_R = 5., Sdir_Deltau = 1., Sdir_Deleta = 1.2;
double Sdir_eta0 = 0.0, Sdir_tau0 = 5., Sdir_etaf = 0.0, Sdir_T = 0.15;
double Sdir_prefactor, Sdir_rt, Sdir_H, Sdir_etat, Sdir_ch_Y_m_eta, Sdir_expon;
double Sdir_term1, Sdir_term2, Sdir_term3;
double S_prefactor = 1.0/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);

double S_direct(double r, double eta, double tau, double MT, double PT, double cos_phi_m_Phi)
{
	//Sdir_prefactor = (2.*Sdir_Jr + 1.)/(twopi*twopi*twopi*hbarC3);
	//Sdir_prefactor = 1./(twopi*twopi*twopi*hbarC3);
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
	//double result = MT * Sdir_ch_Y_m_eta * Sdir_H * exp(Sdir_expon);
	//cerr << r << "   " << eta << "   " << tau << "   " << MT << "   " << PT << "   " << cos_phi_m_Phi << "   " << result << endl;
	return (MT * Sdir_ch_Y_m_eta * Sdir_H * exp(Sdir_expon));
	//return (exp(-r*r));
	//return((r > Sdir_R) ? 0. : 1.);
}

void doHBT::Analyze_sourcefunction(FO_surf* FOsurf_ptr)
{
   double plane_psi = 0.0;
   bool includezeroes = false;
   if (VERBOSE > 0) *global_out_stream_ptr << "Determine nth-order plane angles..." << endl;
   Determine_plane_angle(FOsurf_ptr);
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
   }
   global_plane_psi = plane_psi;

   // begin HBT calculations here...
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      //if (iKT == 0) continue;	//temporary
      //cout << "Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      if (VERBOSE > 0) *global_out_stream_ptr << "   - Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
      beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
      for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
      {
         if (VERBOSE > 1) *global_out_stream_ptr << "\t --> Calculating K_phi = " << K_phi[iKphi] << " ..." << endl;
         Reset_EmissionData();
         SetEmissionData(FOsurf_ptr, K_T[iKT], K_phi[iKphi], includezeroes);
         Get_source_variances(iKT, iKphi);
	 //Update_avgSource_function(iKT, iKphi);
         Calculate_R2_side(iKT, iKphi);
//return;
         Calculate_R2_out(iKT, iKphi);
         Calculate_R2_outside(iKT, iKphi);
         Calculate_R2_long(iKT, iKphi);
         Calculate_R2_sidelong(iKT, iKphi);
         Calculate_R2_outlong(iKT, iKphi);
      }
      R2_Fourier_transform(iKT, plane_psi);
	Svars_Fourier_transform(iKT, global_plane_psi);
   }
   return;
}

void doHBT::dummy_function(int i /*= -1*/, int j /*= -1*/)
{
	if (i < 0 && j < 0) cout << "You ran with the defaults!" << endl;
	else if (i < 0 && j >= 0) cout << "You changed the first parameter!" << endl;
	else if (i >= 0 && j < 0) cout << "You changed the second parameter!" << endl;
	else if (i >= 0 && j >= 0) cout << "You changed both parameters!" << endl;
	return;
}

void doHBT::quick_Analyze_sourcefunction()
{
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      //cout << "Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      R2_Fourier_transform(iKT, global_plane_psi);
   }
   return;
}

void doHBT::Get_azimuthally_averaged_multiplicities()
{
	for (int iKT = 0; iKT < n_localp_T; iKT++)
	{
		azavg_S_func[iKT] = 0;
	
		for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
			azavg_S_func[iKT] += S_func[iKT][iKphi]*K_phi_weight[iKphi];
		
		azavg_S_func[iKT] /= (2.*M_PI);
	}
	
	return;
}

void doHBT::Get_azimuthally_averaged_EBE_HBTradii()
{
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "   - Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
		double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
		beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
		Calculate_azimuthally_averaged_squared_S_func(iKT);
		Calculate_azimuthally_averaged_R2_side(iKT);
		Calculate_azimuthally_averaged_R2_out(iKT);
		Calculate_azimuthally_averaged_R2_long(iKT);
		Calculate_azimuthally_averaged_R2_outside(iKT);
		Calculate_azimuthally_averaged_R2_outlong(iKT);
		Calculate_azimuthally_averaged_R2_sidelong(iKT);
		if (VERBOSE > 0) *global_out_stream_ptr << "   - Finished calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
	}

	return;
}


void doHBT::quick_Analyze_sourcefunction_vars()
{
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      //cout << "Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      Svars_Fourier_transform(iKT, global_plane_psi);
   }
   return;
}

void doHBT::Analyze_AVG_sourcefunction()
{
   //Calculate_avgSource_function();  //redundant
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      if (VERBOSE > 1) *global_out_stream_ptr << "     --> Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
      //double m_perp = sqrt(K_T[iKT]*K_T[iKT] + 0.13957*0.13957);
      beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
      for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
      {
	 //Calculate_avgSource_function(iKT, iKphi);  //comment out if running with Cavg calculations...
         Calculate_avgR2_side(iKT, iKphi);
         Calculate_avgR2_out(iKT, iKphi);
         Calculate_avgR2_outside(iKT, iKphi);
         Calculate_avgR2_long(iKT, iKphi);
         Calculate_avgR2_sidelong(iKT, iKphi);
         Calculate_avgR2_outlong(iKT, iKphi);
      }
      avgR2_Fourier_transform(iKT, 0.);
   }

   return;
}

void doHBT::Analyze_CAVG_sourcefunction()
{
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      if (VERBOSE > 1) *global_out_stream_ptr << "     --> Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
      beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
      for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
      {
         Calculate_CavgR2_side(iKT, iKphi);
         Calculate_CavgR2_out(iKT, iKphi);
         Calculate_CavgR2_outside(iKT, iKphi);
         Calculate_CavgR2_long(iKT, iKphi);
         Calculate_CavgR2_sidelong(iKT, iKphi);
         Calculate_CavgR2_outlong(iKT, iKphi);
      }
      CavgR2_Fourier_transform(iKT, 0.);
   }

   return;
}

void doHBT::Analyze_azimuthally_averaged_AVG_sourcefunction()
{
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      if (VERBOSE > 1) *global_out_stream_ptr << "     --> Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
      beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
      Calculate_azimuthally_averaged_squared_avgS_func(iKT);
      Calculate_azimuthally_averaged_avgR2_side(iKT);
      Calculate_azimuthally_averaged_avgR2_out(iKT);
      Calculate_azimuthally_averaged_avgR2_outside(iKT);
      Calculate_azimuthally_averaged_avgR2_long(iKT);
      Calculate_azimuthally_averaged_avgR2_sidelong(iKT);
      Calculate_azimuthally_averaged_avgR2_outlong(iKT);
   }

   return;
}

void doHBT::Analyze_azimuthally_averaged_CAVG_sourcefunction()
{
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      if (VERBOSE > 1) *global_out_stream_ptr << "     --> Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
      beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
      Calculate_azimuthally_averaged_Cavg_squared_S_func(iKT);
      Calculate_azimuthally_averaged_CavgR2_side(iKT);
      Calculate_azimuthally_averaged_CavgR2_out(iKT);
      Calculate_azimuthally_averaged_CavgR2_outside(iKT);
      Calculate_azimuthally_averaged_CavgR2_long(iKT);
      Calculate_azimuthally_averaged_CavgR2_sidelong(iKT);
      Calculate_azimuthally_averaged_CavgR2_outlong(iKT);
   }

   return;
}

void doHBT::Determine_plane_angle(FO_surf* FOsurf_ptr)
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
         anisotropic_flows_pTdiff[ipt][iorder] = sqrt(sine_iorder[ipt][iorder]*sine_iorder[ipt][iorder] + cosine_iorder[ipt][iorder]*cosine_iorder[ipt][iorder]);
         if( anisotropic_flows_pTdiff[ipt][iorder] < 1e-8)
            anisotropic_flows_pTdiff_psin[ipt][iorder] = 0.0e0;
         else
            anisotropic_flows_pTdiff_psin[ipt][iorder] = atan2(sine_iorder[ipt][iorder], cosine_iorder[ipt][iorder])/double(iorder);
      }
      cosine = cosine/norm;
      sine = sine/norm;
      anisotropic_flows[iorder] = sqrt(sine*sine + cosine*cosine);
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

void doHBT::Get_EdNd3p_cfs(FO_surf* FOsurf_ptr, bool determine_plane_angle /* = true*/)
{
	if (determine_plane_angle)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "Fourier-expanding single-particle spectra..." << endl;
		Determine_plane_angle(FOsurf_ptr);
	}

	if (VERBOSE > 0) *global_out_stream_ptr << "Obtaining Fourier coefficients..." << endl;
	for(int ipt=0; ipt<n_SP_pT; ipt++)
	{
		//double n = dN_dypTdpT[ipt]/(2.*M_PI);
		double n = dN_dypTdpT[ipt];	//already dividing by 2*pi above, so no need to do it again
		EdNd3p_cfs[ipt][0] = n;
		EdNd3p_phases[ipt][0] = 0.0;
		for(int iorder=1; iorder<n_order; iorder++)
		{
			double c = cosine_iorder[ipt][iorder];
			double s = sine_iorder[ipt][iorder];
			EdNd3p_cfs[ipt][iorder] = n*sqrt(c*c + s*s);
			EdNd3p_phases[ipt][iorder] = atan2(s, c)/double(iorder);
		}
	}

   return;
}


void doHBT::SetEmissionData(FO_surf* FOsurf_ptr, double K_T_local, double K_phi_local, bool includezeroes)
{
  double mass = particle_mass;
//mass = 0.001;			//DUMMY TEST
  double mT = sqrt(mass*mass + K_T_local*K_T_local);
  double px = K_T_local*cos(K_phi_local);
  double py = K_T_local*sin(K_phi_local);

  int idx = 0;
  double CDF = 0.0;
  for(int i=0; i<eta_s_npts; i++)
  {
      double local_eta_s = eta_s[i];
      double ch_localetas = ch_eta_s[i];
      double sh_localetas = sh_eta_s[i];

      double p0 = mT*cosh(K_y-local_eta_s);
      double pz = mT*sinh(K_y-local_eta_s);

      for (int j = 0; j < FO_length; j++)
	{
	  //Now that the data is loaded, cycle through it to find the freeze out surface and the emission function.
	  double S_p = 0.0e0;
	  double temp_r = sqrt(FOsurf_ptr[j].xpt*FOsurf_ptr[j].xpt + FOsurf_ptr[j].ypt*FOsurf_ptr[j].ypt);
          double temp_phi = atan2(FOsurf_ptr[j].ypt, FOsurf_ptr[j].xpt);
		if (temp_phi < 0.0)
			temp_phi += 2.*M_PI;
	if (USE_ANALYTIC_S)
	{
		S_p = S_prefactor * S_direct(temp_r, local_eta_s, FOsurf_ptr[j].tau, mT, K_T_local, cos(temp_phi - K_phi_local));
		//if (j == 0)
		//{
		//	cerr << K_phi_local << "   " << K_T_local << endl;
		//}
		if (j==0 && fabs(K_phi_local - M_PI) < 1.e-3 && fabs(K_T_local - 0.05) < 1.e-6)
			cerr << temp_r << "   " << temp_phi << "   " << local_eta_s << "   " << FOsurf_ptr[j].tau
				<< "   " << mT << "   " << K_T_local << "   " << cos(temp_phi - K_phi_local)
				<< "   " << S_p*FOsurf_ptr[j].tau << endl;
	}
	else
	        S_p = Emissionfunction(p0, px, py, pz, &FOsurf_ptr[j]);
        if (flagneg == 1 && S_p < tol)
        {
           S_p = 0.0e0;
        }
	  else
        {
           double S_p_withweight = S_p*FOsurf_ptr[j].tau*eta_s_weight[i];
           (*Emissionfunction_ptr)[idx].data = S_p_withweight; 
           (*Emissionfunction_ptr)[idx].t = FOsurf_ptr[j].tau*ch_localetas;
           (*Emissionfunction_ptr)[idx].tau = FOsurf_ptr[j].tau;
           (*Emissionfunction_ptr)[idx].eta = local_eta_s;
           (*Emissionfunction_ptr)[idx].x = FOsurf_ptr[j].xpt;
           (*Emissionfunction_ptr)[idx].y = FOsurf_ptr[j].ypt;
           (*Emissionfunction_ptr)[idx].r = temp_r;
           (*Emissionfunction_ptr)[idx].phi = temp_phi;
           (*Emissionfunction_ptr)[idx].z = FOsurf_ptr[j].tau*sh_localetas;
           CDF += S_p_withweight;
           (*Emissionfunction_ptr)[idx].CDF_value = CDF;
           idx++;
        }
      }
  }
//cerr << "Emissionfunction_length in SetEmissionData (before) for K_T = " << K_T_local << " and K_phi = " << K_phi_local << " is " << Emissionfunction_length << endl;
if (!includezeroes) Emissionfunction_length = idx;  //only want non-zero values of emission function
//cerr << "Emissionfunction_length in SetEmissionData (after) for K_T = " << K_T_local << " and K_phi = " << K_phi_local << " is " << Emissionfunction_length << endl;

  //nomalize CDF to unity
  //ofstream CDF_check("CDF_check.dat");
  for(int i=0; i<Emissionfunction_length; i++)
  {
     (*Emissionfunction_ptr)[i].CDF_value = (*Emissionfunction_ptr)[i].CDF_value / CDF;
     //CDF_check << scientific << setw(15) << setprecision(8)
     //          << i << "   " << Emissionfunction_Data_CDF[i] << endl;
  }
  //CDF_check.close();
  return;
}

//average over Phi_K, eta_s --> leave x, y, tau, K_T (K_Y == 0)
double doHBT::Average_Emission_Function_on_FOsurface(FO_surf* FOsurf_ptr, int FOcell, int iKT)
{
    double mass = particle_mass;
    double K_T_local = K_T[iKT];
    double mT = sqrt(mass*mass + K_T_local*K_T_local);

    double sum = 0.;

    for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
    {
	double tempsum = 0.;
	double K_phi_local = K_phi[iKphi];
	double px = K_T_local*cos(K_phi_local);
	double py = K_T_local*sin(K_phi_local);

	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		double local_eta_s = eta_s[ieta];
		//double ch_localetas = cosh(local_eta_s);
		//double sh_localetas = sinh(local_eta_s);

		double p0 = mT*cosh(K_y-local_eta_s);
		double pz = mT*sinh(K_y-local_eta_s);
		double S_p = Emissionfunction(p0, px, py, pz, &FOsurf_ptr[FOcell]);
		if (S_p < tol) S_p = 0.0e0;
		double symmetry_factor = 1.0;
		if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
		tempsum += S_p*FOsurf_ptr[FOcell].tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor count for the assumed reflection symmetry along eta
	}

	sum += tempsum*K_phi_weight[iKphi];
    }

    return sum;
}

//average over just eta_s --> leave x, y, tau, Phi_K, K_T (K_Y == 0)
double doHBT::Average_Emission_Function_on_FOsurface(FO_surf* FOsurf_ptr, int FOcell, int iKT, int iKphi)
{
    double mass = particle_mass;
    double K_T_local = K_T[iKT];
    double mT = sqrt(mass*mass + K_T_local*K_T_local);

    double sum = 0.;

    //for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
    //{
	double tempsum = 0.;
	double K_phi_local = K_phi[iKphi];
	double px = K_T_local*cos(K_phi_local);
	double py = K_T_local*sin(K_phi_local);

	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		double local_eta_s = eta_s[ieta];
		//double ch_localetas = cosh(local_eta_s);
		//double sh_localetas = sinh(local_eta_s);

		double p0 = mT*cosh(K_y-local_eta_s);
		double pz = mT*sinh(K_y-local_eta_s);
		double S_p = Emissionfunction(p0, px, py, pz, &FOsurf_ptr[FOcell]);
		if (S_p < tol) S_p = 0.0e0;
		double symmetry_factor = 1.0;
		if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
		tempsum += S_p*FOsurf_ptr[FOcell].tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor count for the assumed reflection symmetry along eta
	}

	//sum += tempsum*K_phi_weight[iKphi];
	sum = tempsum;
    //}

    return sum;
}

void doHBT::Average_sourcefunction_on_FOsurface(FO_surf* FOsurf_ptr)
{
	//avgFOsurf_ptr = new FOsurf[FO_length*n_localp_T];
	int idx = 0;

	for (int iKT = 0; iKT < n_localp_T; iKT++)
	for (int iFOcell = 0; iFOcell < FO_length; iFOcell++)
	{
		if (iFOcell == 0 && VERBOSE > 0) *global_out_stream_ptr << "Averaging over Phi_K, eta_s for K_T = " << K_T[iKT] << endl;
		(*avgFOsurf_ptr)[idx].tau = FOsurf_ptr[iFOcell].tau;
		(*avgFOsurf_ptr)[idx].x = FOsurf_ptr[iFOcell].xpt;
		(*avgFOsurf_ptr)[idx].y = FOsurf_ptr[iFOcell].ypt;
		(*avgFOsurf_ptr)[idx].data = Average_Emission_Function_on_FOsurface(FOsurf_ptr, iFOcell, iKT);
		idx++;
	}

	return;
}

void doHBT::Average_sourcefunction_on_FOsurface(FO_surf* FOsurf_ptr, int iKphi)
{
	//avgFOsurf_ptr = new FOsurf[FO_length*n_localp_T];
	int idx = 0;

	for (int iKT = 0; iKT < n_localp_T; iKT++)
	for (int iFOcell = 0; iFOcell < FO_length; iFOcell++)
	{
		if (iFOcell == 0 && VERBOSE > 0) *global_out_stream_ptr << "Averaging over eta_s for K_T = " << K_T[iKT] << " and K_phi = " << K_phi[iKphi] << endl;
		(*avgFOsurf_ptr)[idx].tau = FOsurf_ptr[iFOcell].tau;
		(*avgFOsurf_ptr)[idx].x = FOsurf_ptr[iFOcell].xpt;
		(*avgFOsurf_ptr)[idx].y = FOsurf_ptr[iFOcell].ypt;
		(*avgFOsurf_ptr)[idx].data = Average_Emission_Function_on_FOsurface(FOsurf_ptr, iFOcell, iKT, iKphi);
		idx++;
	}

	return;
}

void doHBT::Cal_dN_dypTdpTdphi(double** SP_p0, double** SP_px, double** SP_py, double** SP_pz, FO_surf* FOsurf_ptr)
{
   double sign = particle_sign;
   double degen = particle_gspin;
   double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double localmass = particle_mass;

   for(int isurf=0; isurf<FO_length ; isurf++)
   {
//cout << "isurf = " << isurf << endl;
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
	double temp_r = surf->r;
	double temp_phi = surf->phi;

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
	//double f0 = 1./(exp(expon)+sign);       //thermal equilibrium distributions

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

	if (USE_ANALYTIC_S)
	{
		S_p = prefactor * S_direct(temp_r, eta_s[ieta], tau, sqrt(pT*pT + localmass*localmass), pT, cos_phi_m_pphi);
		//if (VERBOSE > 0 && isurf == 0) *global_out_stream_ptr << "  --> Cal_dN_dypTdpTdphi(): computed emission function with S_direct()..." << endl;
		if (/*isurf == 0 && */fabs(SP_pphi[iphi] - 0.00386099) < 1.e-3 && fabs(pT - 0.0180112) < 1.e-6)
		//if (isurf==0)
			cerr << temp_r << "   " << temp_phi << "   " << eta_s[ieta] << "   " << FOsurf_ptr[isurf].tau
				<< "   " << sqrt(pT*pT + localmass*localmass) << "   " << pT << "   " << cos_phi_m_pphi
				<< "   " << S_p*FOsurf_ptr[isurf].tau << endl;
	}
	else
	{
		S_p = prefactor*pdsigma*f0*(1.+deltaf);
	        if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
	        //if (1. + deltaf < 0.0)
			S_p = 0.0e0;
	}
//cout << "S_p = " << S_p << endl;
	double symmetry_factor = 1.0;
	if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
         double S_p_withweight = S_p*tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor count for the assumed reflection symmetry along eta direction
//cout << "(ipt, iphi, ieta) = (" << ipt << ", " << iphi << ", " << ieta << "): " << "dN_dypTdpTdphi[ipt][iphi] = " << dN_dypTdpTdphi[ipt][iphi] << endl;
         dN_dypTdpTdphi[ipt][iphi] += S_p_withweight;
      }
//cout << "dN_dypTdpTdphi[" << ipt << "][" << iphi << "] = " << dN_dypTdpTdphi[ipt][iphi] << endl;
      }
      }
   }
   return;
}

double doHBT::Emissionfunction(double p0, double px, double py, double pz, FO_surf* surf)
{
   double mu = surf->particle_mu[particle_id];
   double sign = particle_sign;
   double degen = particle_gspin;

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

   double expon = (gammaT*(p0*1. - px*vx - py*vy) - mu)/Tdec;
   double f0 = 1./(exp(expon)+sign);       //thermal equilibrium distributions

   //p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
   double pdsigma = p0*da0 + px*da1 + py*da2;

   //viscous corrections
   double Wfactor = p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
   double deltaf = 0.;
   if (use_delta_f)
   {
      deltaf = (1. - sign*f0)*Wfactor/(2.0*Tdec*Tdec*(Edec+Pdec));
   }

   double dN_dyd2pTdphi = 1.0*degen/(8.0*(M_PI*M_PI*M_PI))*pdsigma*f0*(1.+deltaf)/hbarC/hbarC/hbarC;
   if (1. + deltaf < 0.0) dN_dyd2pTdphi = 0.0;
   //out << "Spectral funct = " << dN_dyd2pTdphi << endl;

   return (dN_dyd2pTdphi);
}

void doHBT::Get_source_variances(int iKT, int iKphi)
{
   double phi_K = K_phi[iKphi];
   double eta_odd_factor = 1.0, eta_even_factor = 1.0;
   if (ASSUME_ETA_SYMMETRIC)
   {
	eta_odd_factor = 0.0;
	eta_even_factor = 2.0;
   }
   for(int i=0; i<Emissionfunction_length; i++)
   {
     double r = (*Emissionfunction_ptr)[i].r;
     double phi = (*Emissionfunction_ptr)[i].phi;
     double t = (*Emissionfunction_ptr)[i].t;
     double z = (*Emissionfunction_ptr)[i].z;
     double S_x_K = (*Emissionfunction_ptr)[i].data;
     double sin_phi = sin(phi - phi_K);
     double cos_phi = cos(phi - phi_K);

	S_func[iKT][iKphi] += eta_even_factor*S_x_K;
	xs_S[iKT][iKphi] += eta_even_factor*S_x_K*r*sin_phi;
	xo_S[iKT][iKphi] += eta_even_factor*S_x_K*r*cos_phi;
	xl_S[iKT][iKphi] += eta_odd_factor*S_x_K*z;
	t_S[iKT][iKphi] += eta_even_factor*S_x_K*t;
	xs_t_S[iKT][iKphi] += eta_even_factor*S_x_K*r*t*sin_phi;
	xo_t_S[iKT][iKphi] += eta_even_factor*S_x_K*r*t*cos_phi;
	xl_t_S[iKT][iKphi] += eta_odd_factor*S_x_K*z*t;
	xo_xs_S[iKT][iKphi] += eta_even_factor*S_x_K*r*r*sin_phi*cos_phi;
	xl_xs_S[iKT][iKphi] += eta_odd_factor*S_x_K*z*r*sin_phi;
	xo_xl_S[iKT][iKphi] += eta_odd_factor*S_x_K*z*r*cos_phi;
	xs2_S[iKT][iKphi] += eta_even_factor*S_x_K*r*r*sin_phi*sin_phi;
	xo2_S[iKT][iKphi] += eta_even_factor*S_x_K*r*r*cos_phi*cos_phi;
	xl2_S[iKT][iKphi] += eta_even_factor*S_x_K*z*z;
	t2_S[iKT][iKphi] += eta_even_factor*S_x_K*t*t;
   }

//cerr << "DEBUG: " << S_func[iKT][iKphi] << "  " << xs_S[iKT][iKphi] << "  " << xs2_S[iKT][iKphi] << endl;

return;
}

void doHBT::Get_HBTradii_from_Cbar_and_Cavg()
{
//N.B. - f(-1,-1) means do it for all iKT and iKphi

  for (int i = initial_event; i < initial_event + n_events; i++)
  {
	if (VERBOSE > 0) *global_out_stream_ptr << "Reading in results from event = " << i << endl;
	Set_path(global_runfolder + "/" + global_resultsfolder_stem + "-" + patch::to_string(i));
	Readin_results(i);
	Update_avgSource_function(-1,-1);
	Update_CavgSource_function(-1,-1);
  }

  Calculate_avgSource_function(-1,-1);
  Calculate_CavgSource_function(-1,-1);

	Determine_avgplane_angle();
	//Determine_Cavgplane_angle();

  Analyze_AVG_sourcefunction();
  Analyze_CAVG_sourcefunction();

  return;
}

void doHBT::Get_HBTradii_from_azimuthally_averaged_Cavg_and_Cbar()
{
//N.B. - f(-1) means do it for all iKT

  for (int i = initial_event; i < initial_event + n_events; i++)
  {
	if (VERBOSE > 0) *global_out_stream_ptr << "Reading in results from event = " << i << endl;
	Set_path(global_runfolder + "/" + global_resultsfolder_stem + "-" + patch::to_string(i));
	Readin_results(i);
	Update_avgSource_function(-1, -1);
	Update_CavgSource_function(-1, -1);
  }

  Calculate_avgSource_function(-1, -1);
  Calculate_CavgSource_function(-1, -1);

  //Determine_avgplane_angle();
  Analyze_AVG_sourcefunction();	//azimuthally averaged calculation most convenient in terms of \bar{R}^2_{ij},
				//so need to calculate those first

  Analyze_azimuthally_averaged_AVG_sourcefunction();
  Analyze_azimuthally_averaged_CAVG_sourcefunction();

  return;
}

void doHBT::Simulate_subensemble_averaging()
{
	for (int iM = 0; iM < MmaxSEA; iM++)
	{
		// using built-in random generator:
		random_shuffle ( eventvector.begin(), eventvector.end() );
		
		int nb = total_Nev/n_events;
		for (int ibin = 0; ibin < nb; ibin++)
		{
			Reset_source_variances_and_HBT_radii();
			Get_HBTradii_from_Cbar_and_Cavg_random(ibin);
			Update_subensemble_indexfile(iM, ibin);
			Output_CAVG_random_results(iM, ibin);
		}
	}
}

void doHBT::Simulate_azimuthally_averaged_subensemble_averaging()
{
	for (int i=1; i<=total_Nev; ++i) eventvector.push_back(i);

	for (int iM = 0; iM < MmaxSEA; iM++)
	{
		// using built-in random generator:
		random_shuffle ( eventvector.begin(), eventvector.end() );
		if (VERBOSE > 0) *global_out_stream_ptr << "Beginning binning iteration = " << iM << endl;
		
		int nb = total_Nev/n_events;
		for (int ibin = 0; ibin < nb; ibin++)
		{
			if (VERBOSE > 1) *global_out_stream_ptr << "   - binning #" << iM << ": bin = " << ibin << endl;
			Reset_source_variances_and_HBT_radii();
			Get_HBTradii_from_azimuthally_averaged_Cavg_random(ibin);
			Update_subensemble_indexfile(iM, ibin);
			Output_azimuthally_averaged_CAVG_random_results(iM, ibin);
		}
		if (VERBOSE > 0) *global_out_stream_ptr << "Finished binning iteration = " << iM << endl << endl;
	}
}


void doHBT::Get_HBTradii_from_Cbar_and_Cavg_random(int ibin)
{
  // for i in given bin of events
  for (int i = ibin*n_events; i < (ibin+1)*n_events; i++)
  {
	int ith_event = eventvector[i];
	if (VERBOSE > 2) *global_out_stream_ptr << "Reading in results from event = " << ith_event << endl;
	Set_path(global_runfolder + "/" + global_resultsfolder_stem + "-" + patch::to_string(ith_event));
	Readin_results(ith_event);
	//Update_avgSource_function(-1,-1);
	Update_CavgSource_function(-1,-1);
  }

  //Calculate_avgSource_function(-1,-1);
  Calculate_CavgSource_function(-1,-1);

  //Determine_avgplane_angle();
  //Determine_Cavgplane_angle();

  //Analyze_AVG_sourcefunction();
  Analyze_CAVG_sourcefunction();

  return;
}

void doHBT::Get_HBTradii_from_azimuthally_averaged_Cavg_random(int ibin)
{
  // for i in given bin of events
  for (int i = ibin*n_events; i < (ibin+1)*n_events; i++)
  {
	int ith_event = eventvector[i];
	if (VERBOSE > 2) *global_out_stream_ptr << "Reading in results from event = " << ith_event << endl;
	Set_path(global_runfolder + "/" + global_resultsfolder_stem + "-" + patch::to_string(ith_event));
	Readin_azimuthally_averaged_results(ith_event);
	Update_azavgCavgSource_function(-1);
  }

  //with running sums all accumulated, divide by number of events in sub-ensemble
  Calculate_azavgCavgSource_function(-1);

  //effective averages calculated, use to extract corresponding radii
  Analyze_azimuthally_averaged_CAVG_sourcefunction();

  return;
}

void doHBT::Get_HBTradii_from_Cbar_and_Cavg_fixed_events(vector<int> list_of_events)
{
  // for i in given bin of events
  for (int i = 0; i < (int)list_of_events.size(); i++)
  {
	int ith_event = list_of_events[i];
	if (VERBOSE > 0) *global_out_stream_ptr << "Reading in results from event = " << ith_event << endl;
	Set_path(global_runfolder + "/" + global_resultsfolder_stem + "-" + patch::to_string(ith_event));
	Readin_results(ith_event);
	Update_avgSource_function(-1,-1);
	Update_CavgSource_function(-1,-1);
  }

  Calculate_avgSource_function(-1,-1);
  Calculate_CavgSource_function(-1,-1);

  Determine_avgplane_angle();
  //Determine_Cavgplane_angle();

  Analyze_AVG_sourcefunction();
  Analyze_CAVG_sourcefunction();

  return;
}

void doHBT::Determine_avgplane_angle()
{
	for (int Morder=0; Morder < n_order; Morder++)
	{
		double cos_mK_phi[n_localp_phi], sin_mK_phi[n_localp_phi];
		double vnCOS, vnSIN;
		double denominator = 0.0;
		double numeratorCOS = 0.0;
		double numeratorSIN = 0.0;
		for (int i = 0; i < n_localp_phi; i++)
		{
			cos_mK_phi[i] = cos(Morder*K_phi[i]);
			sin_mK_phi[i] = sin(Morder*K_phi[i]);
		}

		for (int iKT = 0; iKT < n_localp_T; iKT++)
		for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			denominator += K_T[iKT]*K_phi_weight[iKphi]*avgS_func[iKT][iKphi];
			numeratorCOS += K_T[iKT]*cos_mK_phi[iKphi]*K_phi_weight[iKphi]*avgS_func[iKT][iKphi];
			numeratorSIN += K_T[iKT]*sin_mK_phi[iKphi]*K_phi_weight[iKphi]*avgS_func[iKT][iKphi];
		}

		if (Morder == 0)
			avgplane_angle[Morder] = denominator;
		else
		{
			vnCOS = numeratorCOS / denominator;
			vnSIN = numeratorSIN / denominator;
			avgplane_angle[Morder] = atan2(vnSIN, vnCOS) / (double)Morder;
		}
	}

	return;
}

void doHBT::Determine_Cavgplane_angle()
{
	for (int Morder=0; Morder < n_order; Morder++)
	{
		double cos_mK_phi[n_localp_phi], sin_mK_phi[n_localp_phi];
		double vnCOS, vnSIN;
		double denominator = 0.0;
		double numeratorCOS = 0.0;
		double numeratorSIN = 0.0;
		for (int i = 0; i < n_localp_phi; i++)
		{
			cos_mK_phi[i] = cos(Morder*K_phi[i]);
			sin_mK_phi[i] = sin(Morder*K_phi[i]);
		}

		for (int iKT = 0; iKT < n_localp_T; iKT++)
		for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			denominator += K_T[iKT]*K_phi_weight[iKphi]*CavgS_func_squared[iKT][iKphi];
			numeratorCOS += K_T[iKT]*cos_mK_phi[iKphi]*K_phi_weight[iKphi]*CavgS_func_squared[iKT][iKphi];
			numeratorSIN += K_T[iKT]*sin_mK_phi[iKphi]*K_phi_weight[iKphi]*CavgS_func_squared[iKT][iKphi];
		}

		if (Morder == 0)
			Cavgplane_angle[Morder] = denominator;
		else
		{
			vnCOS = numeratorCOS / denominator;
			vnSIN = numeratorSIN / denominator;
			Cavgplane_angle[Morder] = atan2(vnSIN, vnCOS) / (double)Morder;
		}
	}

	return;
}


void doHBT::R2_Fourier_transform(int iKT, double plane_psi)
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
//      cout << K_T[iKT] << "  " << Morder << "  " << R2_side_C[iKT][Morder] << "   " << R2_side_S[iKT][Morder] << "  " << R2_out_C[iKT][Morder] << "  " << R2_out_S[iKT][Morder]
//		<< "  " << R2_outside_C[iKT][Morder] << "   " << R2_outside_S[iKT][Morder] << "  " << R2_long_C[iKT][Morder] << "  " << R2_long_S[iKT][Morder]
//		<< "  " << R2_sidelong_C[iKT][Morder] << "   " << R2_sidelong_S[iKT][Morder] << "  " << R2_outlong_C[iKT][Morder] << "  " << R2_outlong_S[iKT][Morder] << endl;
   }
   return;
}

void doHBT::Svars_Fourier_transform(int iKT, double plane_psi)
{
   for(int Morder=0; Morder<n_order; Morder++)
   {
      double cos_mK_phi[n_localp_phi], sin_mK_phi[n_localp_phi];
      for(int i=0; i<n_localp_phi; i++)
      {
         cos_mK_phi[i] = cos(Morder*(K_phi[i] - plane_psi));
         sin_mK_phi[i] = sin(Morder*(K_phi[i] - plane_psi));
      }
      double tempsum_os_cos = 0.0e0;
      double tempsum_os_sin = 0.0e0;
      double tempsum_sl_cos = 0.0e0;
      double tempsum_sl_sin = 0.0e0;
      double tempsum_ol_cos = 0.0e0;
      double tempsum_ol_sin = 0.0e0;
      double tempsum_s2_cos = 0.0e0;
      double tempsum_s2_sin = 0.0e0;
      double tempsum_o2_cos = 0.0e0;
      double tempsum_o2_sin = 0.0e0;
      double tempsum_l2_cos = 0.0e0;
      double tempsum_l2_sin = 0.0e0;
      double tempsum_t2_cos = 0.0e0;
      double tempsum_t2_sin = 0.0e0;
      double tempsum_st_cos = 0.0e0;
      double tempsum_st_sin = 0.0e0;
      double tempsum_ot_cos = 0.0e0;
      double tempsum_ot_sin = 0.0e0;
      double tempsum_lt_cos = 0.0e0;
      double tempsum_lt_sin = 0.0e0;

      for(int i=0; i<n_localp_phi; i++)
      {
         double norm = S_func[iKT][i];
         tempsum_os_cos += (xo_xs_S[iKT][i]/norm - (xs_S[iKT][i]/norm)*(xo_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_os_sin += (xo_xs_S[iKT][i]/norm - (xs_S[iKT][i]/norm)*(xo_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_sl_cos += (xl_xs_S[iKT][i]/norm - (xl_S[iKT][i]/norm)*(xs_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_sl_sin += (xl_xs_S[iKT][i]/norm - (xl_S[iKT][i]/norm)*(xs_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_ol_cos += (xo_xl_S[iKT][i]/norm - (xl_S[iKT][i]/norm)*(xo_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_ol_sin += (xo_xl_S[iKT][i]/norm - (xl_S[iKT][i]/norm)*(xo_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_ot_cos += (xo_t_S[iKT][i]/norm - (t_S[iKT][i]/norm)*(xo_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_ot_sin += (xo_t_S[iKT][i]/norm - (t_S[iKT][i]/norm)*(xo_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_st_cos += (xs_t_S[iKT][i]/norm - (t_S[iKT][i]/norm)*(xs_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_st_sin += (xs_t_S[iKT][i]/norm - (t_S[iKT][i]/norm)*(xs_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_lt_cos += (xl_t_S[iKT][i]/norm - (xl_S[iKT][i]/norm)*(t_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_lt_sin += (xl_t_S[iKT][i]/norm - (xl_S[iKT][i]/norm)*(t_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_s2_cos += (xs2_S[iKT][i]/norm - (xs_S[iKT][i]/norm)*(xs_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_s2_sin += (xs2_S[iKT][i]/norm - (xs_S[iKT][i]/norm)*(xs_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_o2_cos += (xo2_S[iKT][i]/norm - (xo_S[iKT][i]/norm)*(xo_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_o2_sin += (xo2_S[iKT][i]/norm - (xo_S[iKT][i]/norm)*(xo_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_l2_cos += (xl2_S[iKT][i]/norm - (xl_S[iKT][i]/norm)*(xl_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_l2_sin += (xl2_S[iKT][i]/norm - (xl_S[iKT][i]/norm)*(xl_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
         tempsum_t2_cos += (t2_S[iKT][i]/norm - (t_S[iKT][i]/norm)*(t_S[iKT][i]/norm))*cos_mK_phi[i]*K_phi_weight[i];
         tempsum_t2_sin += (t2_S[iKT][i]/norm - (t_S[iKT][i]/norm)*(t_S[iKT][i]/norm))*sin_mK_phi[i]*K_phi_weight[i];
      }
      xs_t_cos[iKT][Morder] = tempsum_st_cos/(2*M_PI);
      xo_t_cos[iKT][Morder] = tempsum_ot_cos/(2*M_PI);
      xl_t_cos[iKT][Morder] = tempsum_lt_cos/(2*M_PI);
      xo_xs_cos[iKT][Morder] = tempsum_os_cos/(2*M_PI);
      xl_xs_cos[iKT][Morder] = tempsum_sl_cos/(2*M_PI);
      xo_xl_cos[iKT][Morder] = tempsum_ol_cos/(2*M_PI);
      xs2_cos[iKT][Morder] = tempsum_s2_cos/(2*M_PI);
      xo2_cos[iKT][Morder] = tempsum_o2_cos/(2*M_PI);
      xl2_cos[iKT][Morder] = tempsum_l2_cos/(2*M_PI);
      t2_cos[iKT][Morder] = tempsum_t2_cos/(2*M_PI);
      xs_t_sin[iKT][Morder] = tempsum_st_sin/(2*M_PI);
      xo_t_sin[iKT][Morder] = tempsum_ot_sin/(2*M_PI);
      xl_t_sin[iKT][Morder] = tempsum_lt_sin/(2*M_PI);
      xo_xs_sin[iKT][Morder] = tempsum_os_sin/(2*M_PI);
      xl_xs_sin[iKT][Morder] = tempsum_sl_sin/(2*M_PI);
      xo_xl_sin[iKT][Morder] = tempsum_ol_sin/(2*M_PI);
      xs2_sin[iKT][Morder] = tempsum_s2_sin/(2*M_PI);
      xo2_sin[iKT][Morder] = tempsum_o2_sin/(2*M_PI);
      xl2_sin[iKT][Morder] = tempsum_l2_sin/(2*M_PI);
      t2_sin[iKT][Morder] = tempsum_t2_sin/(2*M_PI);
   }
   return;
}

void doHBT::avgR2_Fourier_transform(int iKT, double plane_psi)
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
         temp_sum_side_cos += avgR2_side[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_side_sin += avgR2_side[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_out_cos += avgR2_out[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_out_sin += avgR2_out[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_outside_cos += avgR2_outside[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_outside_sin += avgR2_outside[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_long_cos += avgR2_long[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_long_sin += avgR2_long[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_sidelong_cos += avgR2_sidelong[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_sidelong_sin += avgR2_sidelong[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_outlong_cos += avgR2_outlong[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_outlong_sin += avgR2_outlong[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
      }
      avgR2_side_C[iKT][Morder] = temp_sum_side_cos/(2*M_PI);
      avgR2_side_S[iKT][Morder] = temp_sum_side_sin/(2*M_PI);
      avgR2_out_C[iKT][Morder] = temp_sum_out_cos/(2*M_PI);
      avgR2_out_S[iKT][Morder] = temp_sum_out_sin/(2*M_PI);
      avgR2_outside_C[iKT][Morder] = temp_sum_outside_cos/(2*M_PI);
      avgR2_outside_S[iKT][Morder] = temp_sum_outside_sin/(2*M_PI);
      avgR2_long_C[iKT][Morder] = temp_sum_long_cos/(2*M_PI);
      avgR2_long_S[iKT][Morder] = temp_sum_long_sin/(2*M_PI);
      avgR2_sidelong_C[iKT][Morder] = temp_sum_sidelong_cos/(2*M_PI);
      avgR2_sidelong_S[iKT][Morder] = temp_sum_sidelong_sin/(2*M_PI);
      avgR2_outlong_C[iKT][Morder] = temp_sum_outlong_cos/(2*M_PI);
      avgR2_outlong_S[iKT][Morder] = temp_sum_outlong_sin/(2*M_PI);
      if (VERBOSE > 3) *global_out_stream_ptr << K_T[iKT] << "  " << Morder << "  " << avgR2_side_C[iKT][Morder] << "   " << avgR2_side_S[iKT][Morder] << "  " << avgR2_out_C[iKT][Morder] << "  " << avgR2_out_S[iKT][Morder]
	<< "  " << avgR2_outside_C[iKT][Morder] << "   " << avgR2_outside_S[iKT][Morder] << "  " << avgR2_long_C[iKT][Morder] << "  " << avgR2_long_S[iKT][Morder]
		<< "  " << avgR2_sidelong_C[iKT][Morder] << "   " << avgR2_sidelong_S[iKT][Morder] << "  " << avgR2_outlong_C[iKT][Morder] << "  " << avgR2_outlong_S[iKT][Morder] << endl;
   }
   return;
}

void doHBT::CavgR2_Fourier_transform(int iKT, double plane_psi)
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
         temp_sum_side_cos += CavgR2_side[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_side_sin += CavgR2_side[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_out_cos += CavgR2_out[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_out_sin += CavgR2_out[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_outside_cos += CavgR2_outside[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_outside_sin += CavgR2_outside[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_long_cos += CavgR2_long[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_long_sin += CavgR2_long[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_sidelong_cos += CavgR2_sidelong[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_sidelong_sin += CavgR2_sidelong[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
         temp_sum_outlong_cos += CavgR2_outlong[iKT][i]*cos_mK_phi[i]*K_phi_weight[i];
         temp_sum_outlong_sin += CavgR2_outlong[iKT][i]*sin_mK_phi[i]*K_phi_weight[i];
      }
      CavgR2_side_C[iKT][Morder] = temp_sum_side_cos/(2*M_PI);
      CavgR2_side_S[iKT][Morder] = temp_sum_side_sin/(2*M_PI);
      CavgR2_out_C[iKT][Morder] = temp_sum_out_cos/(2*M_PI);
      CavgR2_out_S[iKT][Morder] = temp_sum_out_sin/(2*M_PI);
      CavgR2_outside_C[iKT][Morder] = temp_sum_outside_cos/(2*M_PI);
      CavgR2_outside_S[iKT][Morder] = temp_sum_outside_sin/(2*M_PI);
      CavgR2_long_C[iKT][Morder] = temp_sum_long_cos/(2*M_PI);
      CavgR2_long_S[iKT][Morder] = temp_sum_long_sin/(2*M_PI);
      CavgR2_sidelong_C[iKT][Morder] = temp_sum_sidelong_cos/(2*M_PI);
      CavgR2_sidelong_S[iKT][Morder] = temp_sum_sidelong_sin/(2*M_PI);
      CavgR2_outlong_C[iKT][Morder] = temp_sum_outlong_cos/(2*M_PI);
      CavgR2_outlong_S[iKT][Morder] = temp_sum_outlong_sin/(2*M_PI);
      if (VERBOSE > 3) *global_out_stream_ptr << K_T[iKT] << "  " << Morder << "  " << CavgR2_side_C[iKT][Morder] << "   " << CavgR2_side_S[iKT][Morder] << "  " << CavgR2_out_C[iKT][Morder] << "  " << CavgR2_out_S[iKT][Morder]
		<< "  " << CavgR2_outside_C[iKT][Morder] << "   " << CavgR2_outside_S[iKT][Morder] << "  " << CavgR2_long_C[iKT][Morder] << "  " << CavgR2_long_S[iKT][Morder]
		<< "  " << CavgR2_sidelong_C[iKT][Morder] << "   " << CavgR2_sidelong_S[iKT][Morder] << "  " << CavgR2_outlong_C[iKT][Morder] << "  " << CavgR2_outlong_S[iKT][Morder] << endl;
   }
   return;
}

//End of file

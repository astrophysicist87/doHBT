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

using namespace std;

void doHBT::Analyze_sourcefunction(FO_surf* FOsurf_ptr)
{
   double plane_psi = 0.0;
   bool includezeroes = false;
   *global_out_stream_ptr << "Determine nth-order plane angles..." << endl;
   Determine_plane_angle(FOsurf_ptr);
   int iorder = 3;
   *global_out_stream_ptr << "Analyzing source function w.r.t " << iorder << " th-order participant plane angle..." << endl;
   *global_out_stream_ptr << "psi = " << plane_angle[iorder] << endl;
   plane_psi = plane_angle[iorder];
   global_plane_psi = plane_psi;
   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      //cout << "Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      *global_out_stream_ptr << "   - Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
      beta_perp = K_T[iKT]/(m_perp*cosh(K_y));
      for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
      {
         *global_out_stream_ptr << "\t --> Calculating K_phi = " << K_phi[iKphi] << " ..." << endl;
         Reset_EmissionData();
         SetEmissionData(FOsurf_ptr, K_T[iKT], K_phi[iKphi], includezeroes);
         Get_source_variances(iKT, iKphi);
	 Update_avgSource_function(iKT, iKphi);
         Calculate_R2_side(iKT, iKphi);
         Calculate_R2_out(iKT, iKphi);
         Calculate_R2_outside(iKT, iKphi);
         Calculate_R2_long(iKT, iKphi);
         Calculate_R2_sidelong(iKT, iKphi);
         Calculate_R2_outlong(iKT, iKphi);
      }
      R2_Fourier_transform(iKT, plane_psi);
   }
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
      cout << "Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
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
      double ch_localetas = cosh(local_eta_s);
      double sh_localetas = sinh(local_eta_s);

      double p0 = mT*cosh(K_y-local_eta_s);
      double pz = mT*sinh(K_y-local_eta_s);

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
           CDF += S_p_withweight;
           (*Emissionfunction_ptr)[idx].CDF_value = CDF;
           idx++;
//        }
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

void doHBT::Output_Emission_Function(int iKT, int iKphi, int folderindex)
{
int coords = 2;
//coords:	0 - output in (x,y,z,t,data) format
//		1 - output in (r,phi,eta,tau,data) format
//		2 - output in (x,y,eta,tau,data) format

	ostringstream filename_stream_S;

	filename_stream_S << path << folderindex << "/Emissionfunction_S_spacetime_kT_" 
			<< fixed << setprecision(2) << setw(3) << K_T[iKT] << "_kphi_"
			<< fixed << setprecision(5) << setw(6) << K_phi[iKphi] << ".dat";
	cout << "Output Emissionfunction to " << filename_stream_S.str() << endl;

	ofstream output_S(filename_stream_S.str().c_str());

cerr << "Emissionfunction_length in Output_Emission_Function is " << Emissionfunction_length << endl;

	for(int i=0; i<Emissionfunction_length; i++)
	{
		if (coords == 2) {
		output_S << scientific << setprecision(8) << setw(15)
				<< (*Emissionfunction_ptr)[i].eta << "   "
				<< (*Emissionfunction_ptr)[i].tau << "   "
				<< (*Emissionfunction_ptr)[i].x << "   "
				<< (*Emissionfunction_ptr)[i].y << "   "
				<< (*Emissionfunction_ptr)[i].data << "   ";
		}
		else if (coords == 1) {
		output_S << scientific << setprecision(8) << setw(15)
				<< (*Emissionfunction_ptr)[i].r << "   "
				<< (*Emissionfunction_ptr)[i].phi << "   "
				<< (*Emissionfunction_ptr)[i].eta << "   "
				<< (*Emissionfunction_ptr)[i].tau << "   "
				<< (*Emissionfunction_ptr)[i].data << "   ";
		}
		else {
		output_S << scientific << setprecision(8) << setw(15)
				<< (*Emissionfunction_ptr)[i].x << "   "
				<< (*Emissionfunction_ptr)[i].y << "   "
				<< (*Emissionfunction_ptr)[i].z << "   "
				<< (*Emissionfunction_ptr)[i].t << "   "
				<< (*Emissionfunction_ptr)[i].data << "   ";
		}
		output_S << endl;
	}

	output_S.close();

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
		tempsum += S_p*FOsurf_ptr[FOcell].tau*eta_s_weight[ieta]*2.0; //2.0 count for the assumed reflection symmetry along eta
	}

	sum += tempsum*K_phi_weight[iKphi];
    }

    return sum;
}

void doHBT::Average_sourcefunction_on_FOsurface(FO_surf* FOsurf_ptr)
{
	//avgFOsurf_ptr = new FOsurf[FO_length*n_localp_T];
	int idx = 0;

	for (int iKT = 0; iKT < n_localp_T; iKT++)
	for (int iFOcell = 0; iFOcell < FO_length; iFOcell++)
	{
		if (iFOcell == 0) *global_out_stream_ptr << "Averaging over Phi_K for K_T = " << K_T[iKT] << endl;
		(*avgFOsurf_ptr)[idx].tau = FOsurf_ptr[iFOcell].tau;
		(*avgFOsurf_ptr)[idx].x = FOsurf_ptr[iFOcell].xpt;
		(*avgFOsurf_ptr)[idx].y = FOsurf_ptr[iFOcell].ypt;
		(*avgFOsurf_ptr)[idx].data = Average_Emission_Function_on_FOsurface(FOsurf_ptr, iFOcell, iKT);
		idx++;
	}

	return;
}

void doHBT::Cal_dN_dypTdpTdphi(double** SP_p0, double** SP_px, double** SP_py, double** SP_pz, FO_surf* FOsurf_ptr)
{
   double sign = particle_sign;
   double degen = particle_gspin;
   double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);

   for(int isurf=0; isurf<FO_length ; isurf++)
   {
//cout << "isurf = " << isurf << endl;
      FO_surf* surf = &FOsurf_ptr[isurf];
      double tau = surf->tau;
      double mu = surf->particle_mu[particle_id];
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
         double deltaf = (1. - sign*f0)*Wfactor*deltaf_prefactor;
//deltaf=0.;	//doing this temporarily for understanding source variance fluctuations

         double S_p = prefactor*pdsigma*f0*(1.+deltaf);
	 if (1. + deltaf < 0.0) S_p = 0.0;
//cout << "S_p = " << S_p << endl;
         double S_p_withweight = S_p*tau*eta_s_weight[ieta]*2.0; //2.0 count for the assumed reflection symmetry along eta direction
//cout << "(ipt, iphi, ieta) = (" << ipt << ", " << iphi << ", " << ieta << "): " << "dN_dypTdpTdphi[ipt][iphi] = " << dN_dypTdpTdphi[ipt][iphi] << endl;
         dN_dypTdpTdphi[ipt][iphi] += S_p_withweight;
      }
//cout << "dN_dydphi[" << ipt << "][" << iphi << "] = " << dN_dypTdpTdphi[ipt][iphi] << endl;
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
   double deltaf = (1. - sign*f0)*Wfactor/(2.0*Tdec*Tdec*(Edec+Pdec));
//deltaf=0.;	//doing this temporarily for understanding source variance fluctuations

   double dN_dyd2pTdphi = 1.0*degen/(8.0*(M_PI*M_PI*M_PI))*pdsigma*f0*(1.+deltaf)/hbarC/hbarC/hbarC;
   if (1. + deltaf < 0.0) dN_dyd2pTdphi = 0.0;
   //out << "Spectral funct = " << dN_dyd2pTdphi << endl;

   return (dN_dyd2pTdphi);
}

void doHBT::Get_source_variances(int iKT, int iKphi)
{
   double phi_K = K_phi[iKphi];
   for(int i=0; i<Emissionfunction_length; i++)
   {
     double r = (*Emissionfunction_ptr)[i].r;
     double phi = (*Emissionfunction_ptr)[i].phi;
     double t = (*Emissionfunction_ptr)[i].t;
     double z = (*Emissionfunction_ptr)[i].z;
     double S_x_K = (*Emissionfunction_ptr)[i].data;
     double sin_phi = sin(phi - phi_K);
     double cos_phi = cos(phi - phi_K);

	S_func[iKT][iKphi] += S_x_K;
	xs_S[iKT][iKphi] += S_x_K*r*sin_phi;
	xo_S[iKT][iKphi] += S_x_K*r*cos_phi;
	xl_S[iKT][iKphi] += S_x_K*z;
	t_S[iKT][iKphi] += S_x_K*t;
	xs_t_S[iKT][iKphi] += S_x_K*r*t*sin_phi;
	xo_t_S[iKT][iKphi] += S_x_K*r*t*cos_phi;
	xl_t_S[iKT][iKphi] += S_x_K*z*t;
	xo_xs_S[iKT][iKphi] += S_x_K*r*r*sin_phi*cos_phi;
	xl_xs_S[iKT][iKphi] += S_x_K*z*r*sin_phi;
	xo_xl_S[iKT][iKphi] += S_x_K*z*r*cos_phi;
	xs2_S[iKT][iKphi] += S_x_K*r*r*sin_phi*sin_phi;
	xo2_S[iKT][iKphi] += S_x_K*r*r*cos_phi*cos_phi;
	xl2_S[iKT][iKphi] += S_x_K*z*z;
	t2_S[iKT][iKphi] += S_x_K*t*t;
   }

return;
}

void doHBT::Get_HBTradii_from_C_ev()
{

  for (int i = 0; i < n_localp_T; i++)
  {
    for (int j = 0; j < n_localp_phi; j++)
    {
	cout << "Getting HBT radii from <C>_{ev} at K_T = " << K_T[i] << ", K_phi = " << K_phi[j] << endl;
	Calculate_HBTradii_from_C_ev(i, j);
    }
    CavgR2_Fourier_transform(i, 0.);
  }

  Output_CAVG_results();
  Analyze_AVG_sourcefunction();
  Output_AVG_results();

  return;
}

void doHBT::Calculate_HBTradii_from_C_ev(int iKT, int iKphi)
{
double sumside = 0.;
double sumout = 0.;
double sumlong = 0.;
double sumoutside = 0.;
double sumsidelong = 0.;
double sumoutlong = 0.;

double fluctuations_term = 0.;
double fluctuations_term_alt = 0.;

//Readin_AVG_results();

//for (int local_folderindex = 1; local_folderindex <= n_events; local_folderindex++)
for (int local_folderindex = initial_event; local_folderindex < initial_event + n_events; local_folderindex++)
{
   Readin_results(local_folderindex);
   Update_avgSource_function(iKT, iKphi);
}

Calculate_avgSource_function(iKT, iKphi);

//for (int local_folderindex = 1; local_folderindex <= n_events; local_folderindex++)
for (int local_folderindex = initial_event; local_folderindex < initial_event + n_events; local_folderindex++)
{
   Readin_results(local_folderindex);

//compute correction factor to cancel out multiplicity fluctuations, to get "true" radii
fluctuations_term += (S_func[iKT][iKphi] - avgS_func[iKT][iKphi]) * (S_func[iKT][iKphi] - avgS_func[iKT][iKphi]) / (avgS_func[iKT][iKphi] * avgS_func[iKT][iKphi]);
fluctuations_term_alt += (S_func[iKT][iKphi] * S_func[iKT][iKphi]) / (avgS_func[iKT][iKphi] * avgS_func[iKT][iKphi]);

//cout << "debug: S_func[" << iKT << "][" << iKphi << "] = " << S_func[iKT][iKphi] << " and avgS_func[" << iKT << "][" << iKphi << "] = " << avgS_func[iKT][iKphi] << endl;
   double prefactor = (S_func[iKT][iKphi] * S_func[iKT][iKphi]) / (avgS_func[iKT][iKphi] * avgS_func[iKT][iKphi]);
   sumside += prefactor * R2_side[iKT][iKphi];
   sumout += prefactor * R2_out[iKT][iKphi];
   sumlong += prefactor * R2_long[iKT][iKphi];
   sumoutside += prefactor * R2_outside[iKT][iKphi];
   sumsidelong += prefactor * R2_sidelong[iKT][iKphi];
   sumoutlong += prefactor * R2_outlong[iKT][iKphi];
}

fluctuations_term /= double(n_events);
fluctuations_term_alt /= double(n_events);

//divide by this correcting normalization factor to get normalized curvature of correlator at q --> 0
double correction_factor = 1 + fluctuations_term;
cout << "correction_factor = " << correction_factor << endl;
cout << "other correction factor (as a check) = " << 1./fluctuations_term_alt << endl;

CavgR2_side[iKT][iKphi] = sumside / double(n_events) / correction_factor;
CavgR2_out[iKT][iKphi] = sumout / double(n_events) / correction_factor;
CavgR2_long[iKT][iKphi] = sumlong / double(n_events) / correction_factor;
CavgR2_outside[iKT][iKphi] = sumoutside / double(n_events) / correction_factor;
CavgR2_outlong[iKT][iKphi] = sumoutlong / double(n_events) / correction_factor;
CavgR2_sidelong[iKT][iKphi] = sumsidelong / double(n_events) / correction_factor;

return;
}

void doHBT::Update_avgSource_function(int iKT, int iKphi)
{
//N.B. - avgs. only contains sums, have not actually been averaged yet
//	for(int iKT = 0; iKT < n_localp_T; iKT++)
//	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
//	{
		avgS_func[iKT][iKphi] += S_func[iKT][iKphi];
		avgxs_S[iKT][iKphi] += xs_S[iKT][iKphi];
		avgxo_S[iKT][iKphi] += xo_S[iKT][iKphi];
		avgxl_S[iKT][iKphi] += xl_S[iKT][iKphi];
		avgt_S[iKT][iKphi] += t_S[iKT][iKphi];
		avgxs_t_S[iKT][iKphi] += xs_t_S[iKT][iKphi];
		avgxo_t_S[iKT][iKphi] += xo_t_S[iKT][iKphi];
		avgxl_t_S[iKT][iKphi] += xl_t_S[iKT][iKphi];
		avgxo_xs_S[iKT][iKphi] += xo_xs_S[iKT][iKphi];
		avgxl_xs_S[iKT][iKphi] += xl_xs_S[iKT][iKphi];
		avgxo_xl_S[iKT][iKphi] += xo_xl_S[iKT][iKphi];
		avgxs2_S[iKT][iKphi] += xs2_S[iKT][iKphi];
		avgxo2_S[iKT][iKphi] += xo2_S[iKT][iKphi];
		avgxl2_S[iKT][iKphi] += xl2_S[iKT][iKphi];
		avgt2_S[iKT][iKphi] += t2_S[iKT][iKphi];
//	}
return;
}

void doHBT::Update_avgSource_function()
{
//N.B. - avgs. only contains sums, have not actually been averaged yet
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		avgS_func[iKT][iKphi] += S_func[iKT][iKphi];
		avgxs_S[iKT][iKphi] += xs_S[iKT][iKphi];
		avgxo_S[iKT][iKphi] += xo_S[iKT][iKphi];
		avgxl_S[iKT][iKphi] += xl_S[iKT][iKphi];
		avgt_S[iKT][iKphi] += t_S[iKT][iKphi];
		avgxs_t_S[iKT][iKphi] += xs_t_S[iKT][iKphi];
		avgxo_t_S[iKT][iKphi] += xo_t_S[iKT][iKphi];
		avgxl_t_S[iKT][iKphi] += xl_t_S[iKT][iKphi];
		avgxo_xs_S[iKT][iKphi] += xo_xs_S[iKT][iKphi];
		avgxl_xs_S[iKT][iKphi] += xl_xs_S[iKT][iKphi];
		avgxo_xl_S[iKT][iKphi] += xo_xl_S[iKT][iKphi];
		avgxs2_S[iKT][iKphi] += xs2_S[iKT][iKphi];
		avgxo2_S[iKT][iKphi] += xo2_S[iKT][iKphi];
		avgxl2_S[iKT][iKphi] += xl2_S[iKT][iKphi];
		avgt2_S[iKT][iKphi] += t2_S[iKT][iKphi];
	}
return;
}

void doHBT::Calculate_avgSource_function(int iKT, int iKphi)
{
//N.B. - avgs. only contains sums, doing averaging here
//	for(int iKT = 0; iKT < n_localp_T; iKT++)
//	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
//	{
		avgS_func[iKT][iKphi] /= double(n_events);
		avgxs_S[iKT][iKphi] /= double(n_events);
		avgxo_S[iKT][iKphi] /= double(n_events);
		avgxl_S[iKT][iKphi] /= double(n_events);
		avgt_S[iKT][iKphi] /= double(n_events);
		avgxs_t_S[iKT][iKphi] /= double(n_events);
		avgxo_t_S[iKT][iKphi] /= double(n_events);
		avgxl_t_S[iKT][iKphi] /= double(n_events);
		avgxo_xs_S[iKT][iKphi] /= double(n_events);
		avgxl_xs_S[iKT][iKphi] /= double(n_events);
		avgxo_xl_S[iKT][iKphi] /= double(n_events);
		avgxs2_S[iKT][iKphi] /= double(n_events);
		avgxo2_S[iKT][iKphi] /= double(n_events);
		avgxl2_S[iKT][iKphi] /= double(n_events);
		avgt2_S[iKT][iKphi] /= double(n_events);
//	}
return;
}

void doHBT::Calculate_avgSource_function()
{
//N.B. - avgs. only contains sums, doing averaging here
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		avgS_func[iKT][iKphi] /= double(n_events);
		avgxs_S[iKT][iKphi] /= double(n_events);
		avgxo_S[iKT][iKphi] /= double(n_events);
		avgxl_S[iKT][iKphi] /= double(n_events);
		avgt_S[iKT][iKphi] /= double(n_events);
		avgxs_t_S[iKT][iKphi] /= double(n_events);
		avgxo_t_S[iKT][iKphi] /= double(n_events);
		avgxl_t_S[iKT][iKphi] /= double(n_events);
		avgxo_xs_S[iKT][iKphi] /= double(n_events);
		avgxl_xs_S[iKT][iKphi] /= double(n_events);
		avgxo_xl_S[iKT][iKphi] /= double(n_events);
		avgxs2_S[iKT][iKphi] /= double(n_events);
		avgxo2_S[iKT][iKphi] /= double(n_events);
		avgxl2_S[iKT][iKphi] /= double(n_events);
		avgt2_S[iKT][iKphi] /= double(n_events);
	}
return;
}

void doHBT::Calculate_R2_side(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xs2_S[iKT][iKphi];
   double term2 = xs_S[iKT][iKphi];

   R2_side[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   return;
}

void doHBT::Calculate_R2_out(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xo2_S[iKT][iKphi] - 2.*beta_perp*xo_t_S[iKT][iKphi] + beta_perp*beta_perp*t2_S[iKT][iKphi];
   double term2 = xo_S[iKT][iKphi] - beta_perp*t_S[iKT][iKphi];

   R2_out[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   return;
}

void doHBT::Calculate_R2_outside(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xo_xs_S[iKT][iKphi] - beta_perp*xs_t_S[iKT][iKphi];
   double term2 = xo_S[iKT][iKphi] - beta_perp*t_S[iKT][iKphi];
   double term3 = xs_S[iKT][iKphi];

   R2_outside[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
   return;
}

void doHBT::Calculate_R2_long(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xl2_S[iKT][iKphi] - 2.*beta_l*xl_t_S[iKT][iKphi] + beta_l*beta_l*t2_S[iKT][iKphi];
   double term2 = xl_S[iKT][iKphi] - beta_l*t_S[iKT][iKphi];

   R2_long[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   return;
}

void doHBT::Calculate_R2_outlong(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xo_xl_S[iKT][iKphi] - beta_perp*xl_t_S[iKT][iKphi] - beta_l*xo_t_S[iKT][iKphi] + beta_perp*beta_l*t2_S[iKT][iKphi];
   double term2 = xo_S[iKT][iKphi] - beta_perp*t_S[iKT][iKphi];
   double term3 = xl_S[iKT][iKphi] - beta_l*t_S[iKT][iKphi];

   R2_outlong[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
   return;
}

void doHBT::Calculate_R2_sidelong(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xl_xs_S[iKT][iKphi] - beta_l*xs_t_S[iKT][iKphi];
   double term2 = xs_S[iKT][iKphi];
   double term3 = xl_S[iKT][iKphi] - beta_l*t_S[iKT][iKphi];

   R2_sidelong[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
   return;
}

void doHBT::Calculate_avgR2_side(int iKT, int iKphi)
{
   double norm = avgS_func[iKT][iKphi];
   double term1 = avgxs2_S[iKT][iKphi];
   double term2 = avgxs_S[iKT][iKphi];

   avgR2_side[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
//debug
//   cout << "avgR2_side[" << K_T[iKT] << "][" << K_phi[iKphi] << "] = " << avgR2_side[iKT][iKphi] << endl;
   return;
}

void doHBT::Calculate_avgR2_out(int iKT, int iKphi)
{
   double norm = avgS_func[iKT][iKphi];
   double term1 = avgxo2_S[iKT][iKphi] - 2.*beta_perp*avgxo_t_S[iKT][iKphi] + beta_perp*beta_perp*avgt2_S[iKT][iKphi];
   double term2 = avgxo_S[iKT][iKphi] - beta_perp*avgt_S[iKT][iKphi];

   avgR2_out[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
//debug
//   cout << "avgR2_out[" << K_T[iKT] << "][" << K_phi[iKphi] << "] = " << avgR2_out[iKT][iKphi] << endl;
//   cout << "avgxo2_S[" << K_T[iKT] << "][" << K_phi[iKphi] << "] = " << avgxo2_S[iKT][iKphi] << endl;
//   cout << "avgxo_t_S[" << K_T[iKT] << "][" << K_phi[iKphi] << "] = " << avgxo_t_S[iKT][iKphi] << endl;
//   cout << "avgt2_S[" << K_T[iKT] << "][" << K_phi[iKphi] << "] = " << avgt2_S[iKT][iKphi] << endl;
//   cout << "beta_perp = " << beta_perp << endl;
   return;
}

void doHBT::Calculate_avgR2_outside(int iKT, int iKphi)
{
   double norm = avgS_func[iKT][iKphi];
   double term1 = avgxo_xs_S[iKT][iKphi] - beta_perp*avgxs_t_S[iKT][iKphi];
   double term2 = avgxo_S[iKT][iKphi] - beta_perp*avgt_S[iKT][iKphi];
   double term3 = avgxs_S[iKT][iKphi];

   avgR2_outside[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
   return;
}

void doHBT::Calculate_avgR2_long(int iKT, int iKphi)
{
   double norm = avgS_func[iKT][iKphi];
   double term1 = avgxl2_S[iKT][iKphi] - 2.*beta_l*avgxl_t_S[iKT][iKphi] + beta_l*beta_l*avgt2_S[iKT][iKphi];
   double term2 = avgxl_S[iKT][iKphi] - beta_l*avgt_S[iKT][iKphi];

   avgR2_long[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   return;
}

void doHBT::Calculate_avgR2_outlong(int iKT, int iKphi)
{
   double norm = avgS_func[iKT][iKphi];
   double term1 = avgxo_xl_S[iKT][iKphi] - beta_perp*avgxl_t_S[iKT][iKphi] - beta_l*avgxo_t_S[iKT][iKphi] + beta_perp*beta_l*avgt2_S[iKT][iKphi];
   double term2 = avgxo_S[iKT][iKphi] - beta_perp*avgt_S[iKT][iKphi];
   double term3 = avgxl_S[iKT][iKphi] - beta_l*avgt_S[iKT][iKphi];

   avgR2_outlong[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
   return;
}

void doHBT::Calculate_avgR2_sidelong(int iKT, int iKphi)
{
   double norm = avgS_func[iKT][iKphi];
   double term1 = avgxl_xs_S[iKT][iKphi] - beta_l*avgxs_t_S[iKT][iKphi];
   double term2 = avgxs_S[iKT][iKphi];
   double term3 = avgxl_S[iKT][iKphi] - beta_l*avgt_S[iKT][iKphi];

   avgR2_sidelong[iKT][iKphi] = term1/norm - term2*term3/(norm*norm);
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
      cout << K_T[iKT] << "  " << Morder << "  " << avgR2_side_C[iKT][Morder] << "   " << avgR2_side_S[iKT][Morder] << "  " << avgR2_out_C[iKT][Morder] << "  " << avgR2_out_S[iKT][Morder]
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
      cout << K_T[iKT] << "  " << Morder << "  " << CavgR2_side_C[iKT][Morder] << "   " << CavgR2_side_S[iKT][Morder] << "  " << CavgR2_out_C[iKT][Morder] << "  " << CavgR2_out_S[iKT][Morder]
		<< "  " << CavgR2_outside_C[iKT][Morder] << "   " << CavgR2_outside_S[iKT][Morder] << "  " << CavgR2_long_C[iKT][Morder] << "  " << CavgR2_long_S[iKT][Morder]
		<< "  " << CavgR2_sidelong_C[iKT][Morder] << "   " << CavgR2_sidelong_S[iKT][Morder] << "  " << CavgR2_outlong_C[iKT][Morder] << "  " << CavgR2_outlong_S[iKT][Morder] << endl;
   }
   return;
}

//End of file

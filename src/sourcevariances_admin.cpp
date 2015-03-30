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

#include "sourcevariances.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

SourceVariances::SourceVariances(particle_info* particle)
{
	//particle information
	particle_name = particle->name;
	particle_mass = particle->mass;
	particle_sign = particle->sign;
	particle_gspin = particle->gspin;
	//particle_id = particle_idx;
	S_prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;

	//just need this for various dummy momentum calculations
	P_eval = new double [4];
	Pp = new double [4];
	Pm = new double [4];

	n_tau_pts = 10;
	n_zeta_pts = 10;
	n_v_pts = 10;
	n_s_pts = 10;
	v_min = -1.;
	v_max = 1.;
	zeta_min = 0.;
	zeta_max = M_PI;
	
   //default: use delta_f in calculations
   use_delta_f = true;
   no_df_stem = "";
   n_resonance = 0;
   Emissionfunction_ptr = new vector<Emissionfunction_data> (1);
//cerr << "made it inside!" << endl;
	//n_weighting_functions = 15;
	n_weighting_functions = 3;  //just doing R^2_s to start
	zvec = new double [4];

	ifstream tempresonancefile("/home/plumberg.1/HBTPlumberg/EOS/temporary_resonance_data.dat");
	tempresonancefile >> n_resonance;
	resonances.resonance_mass = new double [n_resonance];
	resonances.resonance_Gamma = new double [n_resonance];
	resonances.resonance_total_br = new double [n_resonance];
	resonances.resonance_decay_masses = new double* [n_resonance];
	current_resonance_mass = 0.0;
	current_resonance_Gamma = 0.0;
	current_resonance_total_br = 0.0;
	current_resonance_decay_masses = new double [2];
	current_resonance_decay_masses[0] = 0.0;
	current_resonance_decay_masses[1] = 0.0;
	for (int ir=0; ir<n_resonance; ir++)
	{
		resonances.resonance_decay_masses[ir] = new double [2];
		resonances.resonance_decay_masses[ir][0] = 0.0;
		resonances.resonance_decay_masses[ir][1] = 0.0;
	}
	while (!tempresonancefile.eof())
	{
		int row_index = 0;
		tempresonancefile >> row_index;
		//note that we have to convert given table values to GeV
		tempresonancefile >> resonances.resonance_mass[row_index-1];
		tempresonancefile >> resonances.resonance_decay_masses[row_index-1][0];
		tempresonancefile >> resonances.resonance_decay_masses[row_index-1][1];
		tempresonancefile >> resonances.resonance_Gamma[row_index-1];
		tempresonancefile >> resonances.resonance_total_br[row_index-1];
	}
	tempresonancefile.close();
//cerr << "finished reading and processing resonances file..." << endl;
	const double MeVToGeV = 0.001;
	tau_pts = new double * [n_resonance];
	s_pts = new double * [n_resonance];
	tau_wts = new double * [n_resonance];
	s_wts = new double * [n_resonance];
	v_pts = new double [n_v_pts];
	v_wts = new double [n_v_pts];
	zeta_pts = new double [n_zeta_pts];
	zeta_wts = new double [n_zeta_pts];
	for (int ir=0; ir<n_resonance; ir++)
	{
		resonances.resonance_mass[ir] *= MeVToGeV;
		resonances.resonance_decay_masses[ir][0] *= MeVToGeV;
		resonances.resonance_decay_masses[ir][1] *= MeVToGeV;
		resonances.resonance_Gamma[ir] *= MeVToGeV;
		tau_pts[ir] = new double [n_tau_pts];
		tau_wts[ir] = new double [n_tau_pts];
		s_pts[ir] = new double [n_s_pts];
		s_wts[ir] = new double [n_s_pts];
	}

//cerr << "setting gaussian integrations points..." << endl;
	//initialize all gaussian points for resonance integrals
	//syntax: int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[])
	gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0, zeta_min, zeta_max, zeta_pts, zeta_wts);
	gauss_quadrature(n_v_pts, 1, 0.0, 0.0, v_min, v_max, v_pts, v_wts);
	for (int ir = 0; ir < n_resonance; ir++)
	{
		//cerr << "working on resonance #" << ir << "..." << endl;
		double Gamma_temp = resonances.resonance_Gamma[ir];
		double m2_temp = resonances.resonance_decay_masses[ir][0];
		double m3_temp = resonances.resonance_decay_masses[ir][1];
		double M_temp = resonances.resonance_mass[ir];
		double s_min_temp = (m2_temp + m3_temp)*(m2_temp + m3_temp);
		double s_max_temp = (M_temp - particle_mass)*(M_temp - particle_mass);
		// N.B. - this is only really necessary for 3-body decays,
		//			but doesn't cause any problems for 2-body and is easier/simpler to code...
		gauss_quadrature(n_tau_pts, 5, 0.0, 0.0, 0.0, Gamma_temp, tau_pts[ir], tau_wts[ir]);
		gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min_temp, s_max_temp, s_pts[ir], s_wts[ir]);
		/*DEBUG*///cout << ir << "     " << m2_temp << "     " << m3_temp << "     " << M_temp
		/*DEBUG*///	<< "     " << particle_mass << "     " << s_min_temp << "     " << s_max_temp << endl;
		/*DEBUG*///for (int is = 0; is < n_s_pts; is++) cout << "    --> " << ir << "     " << is << "     " << s_pts[ir][is] << endl;
	}
//cerr << "finished all that stuff..." << endl;
	

   //single particle spectra for plane angle determination
   SP_pT = new double [n_SP_pT];
   SP_pT_weight = new double [n_SP_pT];
   gauss_quadrature(n_SP_pT, 1, 0.0, 0.0, SP_pT_min, SP_pT_max, SP_pT, SP_pT_weight);
   SP_pphi = new double [n_SP_pphi];
   SP_pphi_weight = new double [n_SP_pphi];
   gauss_quadrature(n_SP_pphi, 1, 0.0, 0.0, 0.0, 2*M_PI, SP_pphi, SP_pphi_weight);
   SP_p_y = 0.0e0;
   dN_dypTdpTdphi = new double* [n_SP_pT];
   cosine_iorder = new double* [n_SP_pT];
   sine_iorder = new double* [n_SP_pT];
   for(int i=0; i<n_SP_pT; i++)
   {
      dN_dypTdpTdphi[i] = new double [n_SP_pphi];
      cosine_iorder[i] = new double [n_order];
      sine_iorder[i] = new double [n_order];
   }
   dN_dydphi = new double [n_SP_pphi];
   dN_dypTdpT = new double [n_SP_pT];
   pTdN_dydphi = new double [n_SP_pphi];
   for(int i=0; i<n_SP_pphi; i++)
   {
      dN_dydphi[i] = 0.0e0;
      pTdN_dydphi[i] = 0.0e0;
      for(int j=0; j<n_SP_pT; j++) dN_dypTdpTdphi[j][i] = 0.0e0;
   }
   for(int i=0; i<n_SP_pT; i++)
   for(int j=0; j<n_order; j++)
   {
      cosine_iorder[i][j] = 0.0e0;
      sine_iorder[i][j] = 0.0e0;
   }
   for (int i=0; i<n_SP_pT; i++)
	dN_dypTdpT[i] = 0.0e0;
   plane_angle = new double [n_order];

   //pair momentum
   K_T = new double [n_localp_T];
   double dK_T = (localp_T_max - localp_T_min)/(n_localp_T - 1 + 1e-100);
   for(int i=0; i<n_localp_T; i++) K_T[i] = localp_T_min + i*dK_T;
   //K_y = p_y;
   K_y = 0.;
	ch_K_y = cosh(K_y);
	sh_K_y = sinh(K_y);
   beta_l = sh_K_y/ch_K_y;
   K_phi = new double [n_localp_phi];
   K_phi_weight = new double [n_localp_phi];
   gauss_quadrature(n_localp_phi, 1, 0.0, 0.0, localp_phi_min, localp_phi_max, K_phi, K_phi_weight);

	source_variances_array = new double [n_weighting_functions];
	for (int i=0; i<n_weighting_functions; i++) source_variances_array[i] = 0.0;

   //spatial rapidity grid
   eta_s = new double [eta_s_npts];
   eta_s_weight = new double [eta_s_npts];
   gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, eta_s_i, eta_s_f, eta_s, eta_s_weight);
   ch_eta_s = new double [eta_s_npts];
   sh_eta_s = new double [eta_s_npts];
	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		ch_eta_s[ieta] = cosh(eta_s[ieta]);
		sh_eta_s[ieta] = sinh(eta_s[ieta]);
	}

   S_func = new double* [n_localp_T];
   xs_S = new double* [n_localp_T];
   xs2_S = new double* [n_localp_T];

   R2_side = new double* [n_localp_T];

   for(int i=0; i<n_localp_T; i++)
   {
      S_func[i] = new double [n_localp_phi];
      xs_S[i] = new double [n_localp_phi];
      xs2_S[i] = new double [n_localp_phi];
      R2_side[i] = new double [n_localp_phi];
   }

//initialize all source variances and HBT radii/coeffs
for(int i=0; i<n_localp_T; i++)
{
	for(int j=0; j<n_localp_phi; j++)
	{
		S_func[i][j] = 0.;
		xs_S[i][j] = 0.;
		xs2_S[i][j] = 0.;
		R2_side[i][j] = 0.;
	}
}

//cerr << "this part was sketchy though..." << endl;
//set-up array to hold resonance momenta ( P^+ only; get P^- by taking P^-_y = -P^+_y )
/*	Pplus = new momentum_data*** [n_resonance];
	for (int ir = 0; ir < n_resonance; ir++)
	{
		Pplus[ir] = new momentum_data** [n_localp_T];
		for (int iKT = 0; iKT < n_localp_T; iKT++)
		{
			Pplus[ir][iKT] = new momentum_data* [n_localp_phi];
			for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
			{
				Pplus[ir][iKT][iKphi] = new momentum_data [eta_s_npts];
				for (int ieta = 0; ieta < eta_s_npts; ieta++)
				{
					Pplus[ir][iKT][iKphi][ieta].P0 = 0.0;
					Pplus[ir][iKT][iKphi][ieta].Po = 0.0;
					Pplus[ir][iKT][iKphi][ieta].Ps = 0.0;
					Pplus[ir][iKT][iKphi][ieta].Pl = 0.0;
				}
			}
		}
	}*/
//cerr << "...but we apparently made it through!" << endl;
		

   return;
}

void SourceVariances::Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx)
{
   //particle information
   particle_name = particle->name;
   particle_mass = particle->mass;
   particle_sign = particle->sign;
   particle_gspin = particle->gspin;
   particle_id = particle_idx;

   //erase contents of single - and two-particle spectra
   for(int i=0; i<n_SP_pphi; i++)
   {
      dN_dydphi[i] = 0.0e0;
      pTdN_dydphi[i] = 0.0e0;
      for(int j=0; j<n_SP_pT; j++) dN_dypTdpTdphi[j][i] = 0.0e0;
   }
   //erase anisotropic flows
   for(int i=0; i<n_SP_pT; i++)
   for(int j=0; j<n_order; j++)
   {
      cosine_iorder[i][j] = 0.0e0;
      sine_iorder[i][j] = 0.0e0;
   }

   //Emission function
   FO_length = FOarray_length;
   Emissionfunction_length = FO_length*eta_s_npts;
   Emissionfunction_ptr = new vector<Emissionfunction_data> (Emissionfunction_length);
   avgFOsurf_ptr = new vector<Emissionfunction_data> (FO_length*n_localp_T);

	for(int i=0; i<Emissionfunction_length; i++)
	   {
		(*Emissionfunction_ptr)[i].data = 0.0;
		(*Emissionfunction_ptr)[i].t = 0.0;
		(*Emissionfunction_ptr)[i].x = 0.0;
		(*Emissionfunction_ptr)[i].y = 0.0;
		(*Emissionfunction_ptr)[i].z = 0.0;
		(*Emissionfunction_ptr)[i].r = 0.0;
		(*Emissionfunction_ptr)[i].phi = 0.0;
		(*Emissionfunction_ptr)[i].tau = 0.0;
		(*Emissionfunction_ptr)[i].eta = 0.0;
	}

//reset only EBE source variances and EBE HBT radii/coeffs
for(int i=0; i<n_localp_T; i++)
{
	for(int j=0; j<n_localp_phi; j++)
	{
		S_func[i][j] = 0.;
		xs_S[i][j] = 0.;
		xs2_S[i][j] = 0.;

		R2_side[i][j] = 0.;
	}
}

   return;
}

SourceVariances::~SourceVariances()
{
   delete Emissionfunction_ptr;

   delete[] SP_pT;
   delete[] SP_pT_weight;
   delete[] SP_pphi;
   delete[] SP_pphi_weight;
   delete[] dN_dydphi;
   delete[] dN_dypTdpT;
   delete[] pTdN_dydphi;
   for(int i=0; i<n_SP_pT; i++)
   {
      delete[] dN_dypTdpTdphi[i];
      delete[] cosine_iorder[i];
      delete[] sine_iorder[i];
   }
   delete[] dN_dypTdpTdphi;
   delete[] cosine_iorder;
   delete[] sine_iorder;
   delete[] plane_angle;

   delete[] K_T;
   delete[] K_phi;
   delete[] K_phi_weight;
   delete[] eta_s;
   delete[] eta_s_weight;

   for(int i=0; i<n_localp_T; i++)
   {
      delete[] R2_side[i];
   }

   delete[] R2_side;

   for(int i=0; i<n_localp_T; i++)
   {
      delete[] S_func[i];
      delete[] xs_S[i];
      delete[] xs2_S[i];
   }
   delete[] S_func;
   delete[] xs_S;
   delete[] xs2_S;

   return;
}

void SourceVariances::Reset_EmissionData()
{
   Emissionfunction_length = FO_length*eta_s_npts;

   for(int i=0; i<Emissionfunction_length; i++)
   {
		(*Emissionfunction_ptr)[i].data = 0.0;
		(*Emissionfunction_ptr)[i].t = 0.0;
		(*Emissionfunction_ptr)[i].x = 0.0;
		(*Emissionfunction_ptr)[i].y = 0.0;
		(*Emissionfunction_ptr)[i].z = 0.0;
		(*Emissionfunction_ptr)[i].r = 0.0;
		(*Emissionfunction_ptr)[i].phi = 0.0;
		(*Emissionfunction_ptr)[i].tau = 0.0;
		(*Emissionfunction_ptr)[i].eta = 0.0;
   }
}

bool SourceVariances::fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

//print output to output filestream, one line at a time
void SourceVariances::Set_ofstream(ofstream& myout)
{
	global_out_stream_ptr = &myout;

	return;
}

//print output to output filestream, one line at a time
void SourceVariances::Set_path(string localpath)
{
	global_path = localpath;

	return;
}

void SourceVariances::Set_runfolder(string localrunfolder)
{
	global_runfolder = localrunfolder;

	return;
}

void SourceVariances::Set_resultsfolder_stem(string usrdef_resultsfolder_stem)
{
	global_resultsfolder_stem = usrdef_resultsfolder_stem;

	return;
}

void SourceVariances::Set_use_delta_f(bool usrdef_usedeltaf)
{
	use_delta_f = usrdef_usedeltaf;
	if (!use_delta_f)
		no_df_stem = "_no_df";
	return;
}

void SourceVariances::Set_particle_mass(double usrdef_particle_mass)
{
	particle_mass = usrdef_particle_mass;
	return;
}

void SourceVariances::Set_current_FOsurf_ptr(FO_surf* FOsurf_ptr)
{
	current_FOsurf_ptr = FOsurf_ptr;
	return;
}

//End of file

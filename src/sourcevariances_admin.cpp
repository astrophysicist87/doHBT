#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>
#include<algorithm>

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

SourceVariances::SourceVariances(particle_info* particle, particle_info* all_particles_in, int Nparticle_in,
					FO_surf* FOsurf_ptr, vector<int> chosen_resonances_in, int particle_idx, ofstream& myout)
{
	//set ofstream for output file
	global_out_stream_ptr = &myout;
	
	//particle information (both final-state particle used in HBT and all decay decay_channels)
	particle_name = particle->name;
	particle_mass = particle->mass;
	particle_sign = particle->sign;
	particle_gspin = particle->gspin;
	particle_id = particle_idx;
	target_particle_id = particle_id;
	particle_monval = particle->monval;
	S_prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;
	all_particles = all_particles_in;
	for (int icr = 0; icr < (int)chosen_resonances_in.size(); icr++)
		chosen_resonances.push_back(chosen_resonances_in[icr]);
	thermal_pions_only = false;
	Nparticle = Nparticle_in;
	read_in_all_dN_dypTdpTdphi = false;
	output_all_dN_dypTdpTdphi = true;
	currentfolderindex = -1;

	//just need this for various dummy momentum calculations
	P_eval = new double [4];
	Pp = new double [4];
	Pm = new double [4];

	//n_zeta_pts = 25;
	//n_v_pts = 25;
	//n_s_pts = 25;
	n_zeta_pts = 48;
	n_v_pts = 48;
	n_s_pts = 48;
	v_min = -1.;
	v_max = 1.;
	zeta_min = 0.;
	zeta_max = M_PI;
	
	//default: use delta_f in calculations
	use_delta_f = true;
	no_df_stem = "";
	current_FOsurf_ptr = FOsurf_ptr;
	Emissionfunction_ptr = new vector<Emissionfunction_data> (1);
//cerr << "made it inside!" << endl;
	if (INCLUDE_SOURCE_VARIANCES)
		n_weighting_functions = 15;	//AKA, number of SV's
	else
		n_weighting_functions = 1;	//just spectra
	zvec = new double [4];

//****************************************************************************************************
//OLD CODE FOR READING IN SELECTED decay_channels...
	current_resonance_mass = 0.0;
	current_resonance_mu = 0.0;
	current_resonance_Gamma = 0.0;
	current_resonance_total_br = 0.0;
	current_resonance_decay_masses = new double [2];
	current_resonance_decay_masses[0] = 0.0;
	current_resonance_decay_masses[1] = 0.0;
	previous_resonance_particle_id = -1;
	previous_decay_channel_idx = -1;				//different for each decay channel
	previous_resonance_mass = 0.0;
	previous_resonance_Gamma = 0.0;
	previous_resonance_total_br = 0.0;
	if (chosen_resonances.size() == 0)
	{
		n_decay_channels = 1;
		n_resonance = 0;
		thermal_pions_only = true;
		if (VERBOSE > 0) *global_out_stream_ptr << "Thermal pion(+) only!" << endl;
		decay_channels.resonance_particle_id = new int [n_decay_channels];
		decay_channels.resonance_idx = new int [n_decay_channels];
		decay_channels.resonance_mass = new double [n_decay_channels];
		decay_channels.nbody = new int [n_decay_channels];
		decay_channels.resonance_Gamma = new double [n_decay_channels];
		decay_channels.resonance_total_br = new double [n_decay_channels];
		decay_channels.resonance_mu = new double [n_decay_channels];
		decay_channels.resonance_gspin = new double [n_decay_channels];
		decay_channels.resonance_sign = new int [n_decay_channels];
		decay_channels.resonance_decay_masses = new double * [n_decay_channels];
		decay_channels.resonance_name = new string [n_decay_channels];
		decay_channels.include_channel = new bool [n_decay_channels];
	}
	else
	{
		//n_decay_channels is actually total number of decay channels which can generate pions
		//from chosen decay_channels
		n_decay_channels = get_number_of_decay_channels(chosen_resonances, all_particles);
		n_resonance = (int)chosen_resonances.size();
		if (VERBOSE > 0) *global_out_stream_ptr << "Computed n_decay_channels = " << n_decay_channels << endl
							<< "Computed n_resonance = " << n_resonance << endl;
		decay_channels.resonance_particle_id = new int [n_decay_channels];
		decay_channels.resonance_idx = new int [n_decay_channels];
		decay_channels.resonance_mass = new double [n_decay_channels];
		decay_channels.nbody = new int [n_decay_channels];
		decay_channels.resonance_Gamma = new double [n_decay_channels];
		decay_channels.resonance_total_br = new double [n_decay_channels];
		decay_channels.resonance_direct_br = new double [n_decay_channels];
		decay_channels.resonance_mu = new double [n_decay_channels];
		decay_channels.resonance_gspin = new double [n_decay_channels];
		decay_channels.resonance_sign = new int [n_decay_channels];
		decay_channels.resonance_decay_masses = new double * [n_decay_channels];
		decay_channels.resonance_decay_monvals = new double * [n_decay_channels];
		decay_channels.resonance_decay_Gammas = new double * [n_decay_channels];
		decay_channels.resonance_name = new string [n_decay_channels];
		decay_channels.include_channel = new bool [n_decay_channels];
		int temp_idx = 0;
		for (int icr = 0; icr < n_resonance; icr++)
		{
			particle_info particle_temp = all_particles[chosen_resonances[icr]];
			if (VERBOSE > 0) *global_out_stream_ptr << "Loading resonance: " << particle_temp.name
					<< ", chosen_resonances[" << icr << "] = " << chosen_resonances[icr] << endl;
			for (int idecay = 0; idecay < particle_temp.decays; idecay++)
			{
				if (VERBOSE > 0) *global_out_stream_ptr << "Current temp_idx = " << temp_idx << endl;
				if (temp_idx == n_decay_channels)	//i.e., all contributing decay channels have been loaded
					break;
				decay_channels.resonance_name[temp_idx] = particle_temp.name;		// set name of resonance

				//check if effective branching is too small for inclusion in source variances
				bool effective_br_is_too_small = false;
				if (particle_temp.decays_effective_branchratio[idecay] <= 1.e-12)
					effective_br_is_too_small = true;

				decay_channels.resonance_particle_id[temp_idx] = chosen_resonances[icr];	// set index of resonance in all_particles
				decay_channels.resonance_idx[temp_idx] = icr;					// set index of resonance in chosen_resonances
				decay_channels.resonance_decay_masses[temp_idx] = new double [Maxdecaypart];	// Maxdecaypart == 5
				decay_channels.resonance_decay_monvals[temp_idx] = new double [Maxdecaypart];	// Maxdecaypart == 5
				decay_channels.resonance_decay_Gammas[temp_idx] = new double [Maxdecaypart];	// Maxdecaypart == 5


				//*** SETTING RESONANCE DECAY MASSES DIFFERENTLY FOR NEW ANALYZE SF
				for (int ii = 0; ii < Maxdecaypart; ii++)
				{
					decay_channels.resonance_decay_monvals[temp_idx][ii] = particle_temp.decays_part[idecay][ii];
					if (particle_temp.decays_part[idecay][ii] == 0)
					{
						decay_channels.resonance_decay_masses[temp_idx][ii] = 0.0;
						decay_channels.resonance_decay_Gammas[temp_idx][ii] = 0.0;
					}
					else
					{
						int tempID = lookup_particle_id_from_monval(all_particles, Nparticle, particle_temp.decays_part[idecay][ii]);
						decay_channels.resonance_decay_masses[temp_idx][ii] = all_particles[tempID].mass;
						decay_channels.resonance_decay_Gammas[temp_idx][ii] = all_particles[tempID].width;
					}
				}
				decay_channels.resonance_mu[temp_idx] = particle_temp.mu;
				decay_channels.resonance_gspin[temp_idx] = particle_temp.gspin;
				decay_channels.resonance_sign[temp_idx] = particle_temp.sign;
				decay_channels.resonance_mass[temp_idx] = particle_temp.mass;
				decay_channels.nbody[temp_idx] = abs(particle_temp.decays_Npart[idecay]);
				decay_channels.resonance_Gamma[temp_idx] = particle_temp.width;
				decay_channels.resonance_total_br[temp_idx] = particle_temp.decays_effective_branchratio[idecay];
				decay_channels.resonance_direct_br[temp_idx] = particle_temp.decays_branchratio[idecay];
				
				//check if particle lifetime is too long for inclusion in source variances
				bool lifetime_is_too_long = false;
				if (CHECK_FOR_LIFETIME && decay_channels.resonance_Gamma[temp_idx] < hbarC / max_lifetime)
					lifetime_is_too_long = true;		//i.e., for lifetimes longer than 100 fm/c, skip decay channel

				//check if decay channel has too many daughters (nbody >= 4)
				bool too_many_daughters = false;
				if (CHECK_NBODY && decay_channels.nbody[temp_idx] >= 4)
					too_many_daughters = true;

				if (VERBOSE > 0) *global_out_stream_ptr << "Resonance = " << decay_channels.resonance_name[temp_idx] << ", decay channel " << idecay + 1
						<< ": mu=" << decay_channels.resonance_mu[temp_idx]
						<< ", gs=" << decay_channels.resonance_gspin[temp_idx] << ", sign=" << decay_channels.resonance_sign[temp_idx]
						<< ", M=" << decay_channels.resonance_mass[temp_idx] << ", nbody=" << decay_channels.nbody[temp_idx]
						<< ", Gamma=" << decay_channels.resonance_Gamma[temp_idx] << ", total br=" << decay_channels.resonance_total_br[temp_idx]
						<< ", direct br=" << decay_channels.resonance_direct_br[temp_idx] << endl;

				//set daughter particles masses for each decay channel
				//currently assuming no more than nbody = 3

				if (VERBOSE > 0) *global_out_stream_ptr << "Resonance = " << decay_channels.resonance_name[temp_idx] << ": ";
				for (int decay_part_idx = 0; decay_part_idx < decay_channels.nbody[temp_idx]; decay_part_idx++)
					if (VERBOSE > 0) *global_out_stream_ptr << "m[" << decay_part_idx << "] = "
						<< decay_channels.resonance_decay_masses[temp_idx][decay_part_idx] << "   "
						<< decay_channels.resonance_decay_monvals[temp_idx][decay_part_idx] << "   ";
				if (VERBOSE > 0) *global_out_stream_ptr << endl << endl;

				// if decay channel parent resonance is not too long-lived
				// and decay channel contains at least one target daughter particle,
				// include channel
				if (DO_ALL_DECAY_CHANNELS)
					decay_channels.include_channel[temp_idx] = true;
				else
					decay_channels.include_channel[temp_idx] = (!lifetime_is_too_long
											&& !too_many_daughters
											&& !effective_br_is_too_small);

				temp_idx++;
			}
		}
	}
	
	//*****************************************************************
	// Only make dN_dypTdpTdphi_moments large enough to hold all necessary resonance, not decay channels
	//*****************************************************************
	dN_dypTdpTdphi_moments = new double *** [Nparticle];
	ln_dN_dypTdpTdphi_moments = new double *** [Nparticle];
	sign_of_dN_dypTdpTdphi_moments = new double *** [Nparticle];
	for (int ir=0; ir<Nparticle; ir++)
	{
		dN_dypTdpTdphi_moments[ir] = new double ** [n_weighting_functions];
		ln_dN_dypTdpTdphi_moments[ir] = new double ** [n_weighting_functions];
		sign_of_dN_dypTdpTdphi_moments[ir] = new double ** [n_weighting_functions];
		for (int wfi=0; wfi<n_weighting_functions; wfi++)
		{
			dN_dypTdpTdphi_moments[ir][wfi] = new double * [n_interp2_pT_pts];
			ln_dN_dypTdpTdphi_moments[ir][wfi] = new double * [n_interp2_pT_pts];
			sign_of_dN_dypTdpTdphi_moments[ir][wfi] = new double * [n_interp2_pT_pts];
			for (int ipT=0; ipT<n_interp2_pT_pts; ipT++)
			{
				dN_dypTdpTdphi_moments[ir][wfi][ipT] = new double [n_interp2_pphi_pts];
				ln_dN_dypTdpTdphi_moments[ir][wfi][ipT] = new double [n_interp2_pphi_pts];
				sign_of_dN_dypTdpTdphi_moments[ir][wfi][ipT] = new double [n_interp2_pphi_pts];
				for (int ipphi=0; ipphi<n_interp2_pphi_pts; ipphi++)
					{
						dN_dypTdpTdphi_moments[ir][wfi][ipT][ipphi] = 0.0;
						ln_dN_dypTdpTdphi_moments[ir][wfi][ipT][ipphi] = 0.0;
						sign_of_dN_dypTdpTdphi_moments[ir][wfi][ipT][ipphi] = 0.0;
					}
			}
		}
	}

	s_pts = new double * [n_decay_channels];
	s_wts = new double * [n_decay_channels];
	v_pts = new double [n_v_pts];
	v_wts = new double [n_v_pts];
	zeta_pts = new double [n_zeta_pts];
	zeta_wts = new double [n_zeta_pts];
	for (int idc=0; idc<n_decay_channels; idc++)
	{
		s_pts[idc] = new double [n_s_pts];
		s_wts[idc] = new double [n_s_pts];
	}

//cerr << "setting gaussian integrations points..." << endl;
	//initialize all gaussian points for resonance integrals
	//syntax: int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[])
	gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0, zeta_min, zeta_max, zeta_pts, zeta_wts);
	gauss_quadrature(n_v_pts, 1, 0.0, 0.0, v_min, v_max, v_pts, v_wts);
	for (int idc = 0; idc < n_decay_channels; idc++)
	{
		//cerr << "working on resonance #" << idc << "..." << endl;
		double Gamma_temp = decay_channels.resonance_Gamma[idc];
		double m2_temp = decay_channels.resonance_decay_masses[idc][0];
		double m3_temp = decay_channels.resonance_decay_masses[idc][1];
		double M_temp = decay_channels.resonance_mass[idc];
		double s_min_temp = (m2_temp + m3_temp)*(m2_temp + m3_temp);
		double s_max_temp = (M_temp - particle_mass)*(M_temp - particle_mass);
		// N.B. - this is only really necessary for 3-body decays,
		//			but doesn't cause any problems for 2-body and is easier/simpler to code...
		gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min_temp, s_max_temp, s_pts[idc], s_wts[idc]);
	}
	//for (int iv = 0; iv < n_v_pts; iv++)
	//	*global_out_stream_ptr << setw(8) << setprecision(15) << "v_pts[" << iv << "] = " << v_pts[iv] << " and v_wts[" << iv << "] = " << v_wts[iv] << endl;
	//for (int izeta = 0; izeta < n_zeta_pts; izeta++)
	//	*global_out_stream_ptr << setw(8) << setprecision(15) << "zeta_pts[" << izeta << "] = " << zeta_pts[izeta] << " and zeta_wts[" << izeta << "] = " << zeta_wts[izeta] << endl;
//cerr << "finished all that stuff..." << endl;
	

   //single particle spectra for plane angle determination
   SP_pT = new double [n_SP_pT];
   SP_pT_weight = new double [n_SP_pT];
   gauss_quadrature(n_SP_pT, 1, 0.0, 0.0, SP_pT_min, SP_pT_max, SP_pT, SP_pT_weight);
   SP_pphi = new double [n_SP_pphi];
   SP_pphi_weight = new double [n_SP_pphi];
   gauss_quadrature(n_SP_pphi, 1, 0.0, 0.0, 0.0, 2.*M_PI, SP_pphi, SP_pphi_weight);
   SP_p_y = 0.0e0;

//initialize and set evenly spaced grid of px-py points in transverse plane,
//and corresponding p0 and pz points
	SPinterp2_pT = new double [n_interp2_pT_pts];
	SPinterp2_pphi = new double [n_interp2_pphi_pts];
	sin_SPinterp2_pphi = new double [n_interp2_pphi_pts];
	cos_SPinterp2_pphi = new double [n_interp2_pphi_pts];
	SPinterp2_p0 = new double* [n_interp2_pT_pts];
	SPinterp2_pz = new double* [n_interp2_pT_pts];
	double * dummywts3 = new double [n_interp2_pT_pts];
	double * dummywts4 = new double [n_interp2_pphi_pts];
	for(int ipt=0; ipt<n_interp2_pT_pts; ipt++)
	{
		SPinterp2_p0[ipt] = new double [eta_s_npts];
		SPinterp2_pz[ipt] = new double [eta_s_npts];
	}
	if (UNIFORM_SPACING)
	{
		//use uniformly spaced points in transverse momentum to make
		//interpolation simpler/faster
		for(int ipt=0; ipt<n_interp2_pT_pts; ipt++)
			SPinterp2_pT[ipt] = interp2_pT_min + (double)ipt*Del2_pT;
		for(int ipphi=0; ipphi<n_interp2_pphi_pts; ipphi++)
		{
			SPinterp2_pphi[ipphi] = interp2_pphi_min + (double)ipphi*Del2_pphi;
			sin_SPinterp2_pphi[ipphi] = sin(SPinterp2_pphi[ipphi]);
			cos_SPinterp2_pphi[ipphi] = cos(SPinterp2_pphi[ipphi]);
		}
		//for(int ipt=0; ipt<n_interp2_pT_pts; ipt++)
		//	cerr << SPinterp2_pT[ipt] << endl;
	}
	else
	{
		//try just gaussian points...
		//syntax:
		//int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[]) 
		gauss_quadrature(n_interp2_pT_pts, 5, 0.0, 0.0, 0.0, 13.0, SPinterp2_pT, dummywts3);
		//for(int ipT=0; ipT<n_interp2_pT_pts; ipT++)
		//	if (VERBOSE > 0) *global_out_stream_ptr << "SPinterp2_pT[" << ipT << "] = " << SPinterp2_pT[ipT] << endl;
		//logspace(SPinterp2_pT, interp2_pT_min, interp2_pT_max, n_interp2_pT_pts);
		//scalepoints(SPinterp2_pT, interp2_pT_min, interp2_pT_max, 0.25, n_interp2_pT_pts);
		//gauss_quadrature(n_interp2_pT_pts, 1, 0.0, 0.0, interp2_pT_min, interp2_pT_max, SPinterp2_pT, dummywts3);
		//for(int ipt=0; ipt<n_interp2_pT_pts; ipt++)
		//	cerr << SPinterp2_pT[ipt] << endl;
		gauss_quadrature(n_interp2_pphi_pts, 1, 0.0, 0.0, interp2_pphi_min, interp2_pphi_max, SPinterp2_pphi, dummywts4);
		for(int ipphi=0; ipphi<n_interp2_pphi_pts; ipphi++)
		{
			//SPinterp2_pphi[ipphi] = interp2_pphi_min + (double)ipphi*Del2_pphi;
			//if (VERBOSE > 0) *global_out_stream_ptr << setw(8) << setprecision(15) << "SPinterp2_pphi[" << ipphi << "] = " << SPinterp2_pphi[ipphi] << endl;
			sin_SPinterp2_pphi[ipphi] = sin(SPinterp2_pphi[ipphi]);
			cos_SPinterp2_pphi[ipphi] = cos(SPinterp2_pphi[ipphi]);
		}
	}

	dN_dypTdpTdphi = new double* [n_SP_pT];
	SV_dN_dypTdpTdphi = new double* [n_SP_pT];
	cosine_iorder = new double* [n_SP_pT];
	sine_iorder = new double* [n_SP_pT];
	for(int i=0; i<n_SP_pT; i++)
	{
		dN_dypTdpTdphi[i] = new double [n_SP_pphi];
		SV_dN_dypTdpTdphi[i] = new double [n_SP_pphi];
		cosine_iorder[i] = new double [n_order];
		sine_iorder[i] = new double [n_order];
	}
	dN_dydphi = new double [n_SP_pphi];
	dN_dypTdpT = new double [n_SP_pT];
	pTdN_dydphi = new double [n_SP_pphi];
	SV_dN_dydphi = new double [n_SP_pphi];
	SV_dN_dypTdpT = new double [n_SP_pT];
	SV_pTdN_dydphi = new double [n_SP_pphi];
	for(int i=0; i<n_SP_pphi; i++)
	{
		dN_dydphi[i] = 0.0e0;
		pTdN_dydphi[i] = 0.0e0;
		SV_dN_dydphi[i] = 0.0e0;
		SV_pTdN_dydphi[i] = 0.0e0;
		for(int j=0; j<n_SP_pT; j++)
		{
			dN_dypTdpTdphi[j][i] = 0.0e0;
			SV_dN_dypTdpTdphi[j][i] = 0.0e0;
		}
	}
	for(int i=0; i<n_SP_pT; i++)
	{
		dN_dypTdpT[i] = 0.0e0;
		SV_dN_dypTdpT[i] = 0.0e0;
		for(int j=0; j<n_order; j++)
		{
			cosine_iorder[i][j] = 0.0e0;
			sine_iorder[i][j] = 0.0e0;
		}
	}
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
	xo_S = new double* [n_localp_T];
	xo2_S = new double* [n_localp_T];
	xl_S = new double* [n_localp_T];
	xl2_S = new double* [n_localp_T];
	t_S = new double* [n_localp_T];
	t2_S = new double* [n_localp_T];
	xo_xs_S = new double* [n_localp_T];
	xl_xs_S = new double* [n_localp_T];
	xs_t_S = new double* [n_localp_T];
	xo_xl_S = new double* [n_localp_T];
	xo_t_S = new double* [n_localp_T];
	xl_t_S = new double* [n_localp_T];

	R2_side = new double* [n_localp_T];
	R2_side_C = new double* [n_localp_T];
	R2_side_S = new double* [n_localp_T];
	R2_out = new double* [n_localp_T];
	R2_out_C = new double* [n_localp_T];
	R2_out_S = new double* [n_localp_T];
	R2_long = new double* [n_localp_T];
	R2_long_C = new double* [n_localp_T];
	R2_long_S = new double* [n_localp_T];
	R2_outside = new double* [n_localp_T];
	R2_outside_C = new double* [n_localp_T];
	R2_outside_S = new double* [n_localp_T];
	R2_sidelong = new double* [n_localp_T];
	R2_sidelong_C = new double* [n_localp_T];
	R2_sidelong_S = new double* [n_localp_T];
	R2_outlong = new double* [n_localp_T];
	R2_outlong_C = new double* [n_localp_T];
	R2_outlong_S = new double* [n_localp_T];

	for(int i=0; i<n_localp_T; i++)
	{
		S_func[i] = new double [n_localp_phi];
		xs_S[i] = new double [n_localp_phi];
		xs2_S[i] = new double [n_localp_phi];
		xo_S[i] = new double [n_localp_phi];
		xo2_S[i] = new double [n_localp_phi];
		xl_S[i] = new double [n_localp_phi];
		xl2_S[i] = new double [n_localp_phi];
		t_S[i] = new double [n_localp_phi];
		t2_S[i] = new double [n_localp_phi];
		xo_xs_S[i] = new double [n_localp_phi];
		xl_xs_S[i] = new double [n_localp_phi];
		xs_t_S[i] = new double [n_localp_phi];
		xo_xl_S[i] = new double [n_localp_phi];
		xo_t_S[i] = new double [n_localp_phi];
		xl_t_S[i] = new double [n_localp_phi];
		
		R2_side[i] = new double [n_localp_phi];
		R2_side_C[i] = new double [n_order];
		R2_side_S[i] = new double [n_order];
		R2_out[i] = new double [n_localp_phi];
		R2_out_C[i] = new double [n_order];
		R2_out_S[i] = new double [n_order];
		R2_outside[i] = new double [n_localp_phi];
		R2_outside_C[i] = new double [n_order];
		R2_outside_S[i] = new double [n_order];
		R2_long[i] = new double [n_localp_phi];
		R2_long_C[i] = new double [n_order];
		R2_long_S[i] = new double [n_order];
		R2_sidelong[i] = new double [n_localp_phi];
		R2_sidelong_C[i] = new double [n_order];
		R2_sidelong_S[i] = new double [n_order];
		R2_outlong[i] = new double [n_localp_phi];
		R2_outlong_C[i] = new double [n_order];
		R2_outlong_S[i] = new double [n_order];
	}

	//initialize all source variances and HBT radii/coeffs
	for(int i=0; i<n_localp_T; i++)
	{
		for(int j=0; j<n_localp_phi; j++)
		{
			S_func[i][j] = 0.;
			xs_S[i][j] = 0.;
			xs2_S[i][j] = 0.;
			xo_S[i][j] = 0.;
			xo2_S[i][j] = 0.;
			xl_S[i][j] = 0.;
			xl2_S[i][j] = 0.;
			t_S[i][j] = 0.;
			t2_S[i][j] = 0.;
			xo_xs_S[i][j] = 0.;
			xl_xs_S[i][j] = 0.;
			xs_t_S[i][j] = 0.;
			xo_xl_S[i][j] = 0.;
			xo_t_S[i][j] = 0.;
			xl_t_S[i][j] = 0.;
			
			R2_side[i][j] = 0.;
			R2_out[i][j] = 0.;
			R2_long[i][j] = 0.;
			R2_outside[i][j] = 0.;
			R2_sidelong[i][j] = 0.;
			R2_outlong[i][j] = 0.;
		}
		for(int j=0; j<n_order; j++)
		{
			R2_side_C[i][j] = 0.;
			R2_side_S[i][j] = 0.;
			R2_out_C[i][j] = 0.;
			R2_out_S[i][j] = 0.;
			R2_outside_C[i][j] = 0.;
			R2_outside_S[i][j] = 0.;
			R2_long_C[i][j] = 0.;
			R2_long_S[i][j] = 0.;
			R2_sidelong_C[i][j] = 0.;
			R2_sidelong_S[i][j] = 0.;
			R2_outlong_C[i][j] = 0.;
			R2_outlong_S[i][j] = 0.;
		}
	}

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

   FO_length = FOarray_length;

//reset only EBE source variances and EBE HBT radii/coeffs
for(int i=0; i<n_localp_T; i++)
{
	for(int j=0; j<n_localp_phi; j++)
	{
		S_func[i][j] = 0.;
		xs_S[i][j] = 0.;
		xs2_S[i][j] = 0.;
		xo_S[i][j] = 0.;
		xo2_S[i][j] = 0.;
		xl_S[i][j] = 0.;
		xl2_S[i][j] = 0.;
		t_S[i][j] = 0.;
		t2_S[i][j] = 0.;
		xo_xs_S[i][j] = 0.;
		xl_xs_S[i][j] = 0.;
		xs_t_S[i][j] = 0.;
		xo_xl_S[i][j] = 0.;
		xo_t_S[i][j] = 0.;
		xl_t_S[i][j] = 0.;
		
		R2_side[i][j] = 0.;
		R2_out[i][j] = 0.;
		R2_long[i][j] = 0.;
		R2_outside[i][j] = 0.;
		R2_sidelong[i][j] = 0.;
		R2_outlong[i][j] = 0.;
	}
	for(int j=0; j<n_order; j++)
	{
		R2_side_C[i][j] = 0.;
		R2_side_S[i][j] = 0.;
		R2_out_C[i][j] = 0.;
		R2_out_S[i][j] = 0.;
		R2_outside_C[i][j] = 0.;
		R2_outside_S[i][j] = 0.;
		R2_long_C[i][j] = 0.;
		R2_long_S[i][j] = 0.;
		R2_sidelong_C[i][j] = 0.;
		R2_sidelong_S[i][j] = 0.;
		R2_outlong_C[i][j] = 0.;
		R2_outlong_S[i][j] = 0.;
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

void SourceVariances::Allocate_decay_channel_info()
{
	if (VERBOSE > 2) *global_out_stream_ptr << "Reallocating memory for decay channel information..." << endl;
	VEC_n2_v_factor = new double [n_v_pts];
	VEC_n2_zeta_factor = new double * [n_v_pts];
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
			VEC_n2_Pp[iv][izeta] = new double [4];
			VEC_n2_alpha[iv][izeta] = new double [4];
			VEC_n2_Pm[iv][izeta] = new double [4];
			VEC_n2_alpha_m[iv][izeta] = new double [4];
		}
	}
	NEW_s_pts = new double [n_s_pts];
	NEW_s_wts = new double [n_s_pts];
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
	for(int is = 0; is < n_s_pts; is++)
	{
		VEC_v_factor[is] = new double [n_v_pts];
		VEC_zeta_factor[is] = new double * [n_v_pts];
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
				VEC_Pp[is][iv][izeta] = new double [4];
				VEC_alpha[is][iv][izeta] = new double [4];
				VEC_Pm[is][iv][izeta] = new double [4];
				VEC_alpha_m[is][iv][izeta] = new double [4];
			}
		}
	}
	if (VERBOSE > 2) *global_out_stream_ptr << "Reallocated memory for decay channel information." << endl;

	return;
}

void SourceVariances::Delete_decay_channel_info()
{
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleting memory for decay channel information..." << endl;
	for(int iv = 0; iv < n_v_pts; iv++)
	{
		for(int izeta = 0; izeta < n_zeta_pts; izeta++)
		{
			delete [] VEC_n2_Pp[iv][izeta];
			delete [] VEC_n2_alpha[iv][izeta];
			delete [] VEC_n2_Pm[iv][izeta];
			delete [] VEC_n2_alpha_m[iv][izeta];
		}
		delete [] VEC_n2_MT[iv];
		delete [] VEC_n2_PPhi_tilde[iv];
		delete [] VEC_n2_PPhi_tildeFLIP[iv];
		delete [] VEC_n2_PT[iv];
		delete [] VEC_n2_Pp[iv];
		delete [] VEC_n2_alpha[iv];
		delete [] VEC_n2_Pm[iv];
		delete [] VEC_n2_alpha_m[iv];
		delete [] VEC_n2_zeta_factor[iv];
	}
	delete [] VEC_n2_v_factor;
	delete [] VEC_n2_zeta_factor;
	delete [] VEC_n2_P_Y;
	delete [] VEC_n2_MTbar ;
	delete [] VEC_n2_DeltaMT;
	delete [] VEC_n2_MTp;
	delete [] VEC_n2_MTm;
	delete [] VEC_n2_MT;
	delete [] VEC_n2_PPhi_tilde;
	delete [] VEC_n2_PPhi_tildeFLIP;
	delete [] VEC_n2_PT;
	delete [] VEC_n2_Pp;
	delete [] VEC_n2_alpha;
	delete [] VEC_n2_Pm;
	delete [] VEC_n2_alpha_m;

	for(int is = 0; is < n_s_pts; is++)
	{
		for(int iv = 0; iv < n_v_pts; iv++)
		{
			for(int izeta = 0; izeta < n_zeta_pts; izeta++)
			{
				delete [] VEC_Pp[is][iv][izeta];
				delete [] VEC_alpha[is][iv][izeta];
				delete [] VEC_Pm[is][iv][izeta];
				delete [] VEC_alpha_m[is][iv][izeta];
			}
			delete [] VEC_MT[is][iv];
			delete [] VEC_PPhi_tilde[is][iv];
			delete [] VEC_PPhi_tildeFLIP[is][iv];
			delete [] VEC_PT[is][iv];
			delete [] VEC_Pp[is][iv];
			delete [] VEC_alpha[is][iv];
			delete [] VEC_Pm[is][iv];
			delete [] VEC_alpha_m[is][iv];
			delete [] VEC_zeta_factor[is][iv];
		}
		delete [] VEC_v_factor[is];
		delete [] VEC_zeta_factor[is];
		delete [] VEC_P_Y[is];
		delete [] VEC_MTbar[is];
		delete [] VEC_DeltaMT[is];
		delete [] VEC_MTp[is];
		delete [] VEC_MTm[is];
		delete [] VEC_MT[is];
		delete [] VEC_PPhi_tilde[is];
		delete [] VEC_PPhi_tildeFLIP[is];
		delete [] VEC_PT[is];
		delete [] VEC_Pp[is];
		delete [] VEC_alpha[is];
		delete [] VEC_Pm[is];
		delete [] VEC_alpha_m[is];
	}
	delete [] NEW_s_pts;
	delete [] NEW_s_wts;
	delete [] VEC_pstar;
	delete [] VEC_Estar;
	delete [] VEC_DeltaY;
	delete [] VEC_g_s;
	delete [] VEC_s_factor;
	delete [] VEC_v_factor;
	delete [] VEC_zeta_factor;
	delete [] VEC_Yp;
	delete [] VEC_Ym;
	delete [] VEC_P_Y;
	delete [] VEC_MTbar;
	delete [] VEC_DeltaMT;
	delete [] VEC_MTp;
	delete [] VEC_MTm;
	delete [] VEC_MT;
	delete [] VEC_PPhi_tilde;
	delete [] VEC_PPhi_tildeFLIP;
	delete [] VEC_PT;
	delete [] VEC_Pp;
	delete [] VEC_alpha;
	delete [] VEC_Pm;
	delete [] VEC_alpha_m;
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleted memory for decay channel information." << endl;

	return;
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
	
	/*surf_damu = new double * [3];
	surf_pimunu = new double * [8];  //all pimunu's and bulkPi coeff --> latter probably zero
	surf_Bn_muS_muB = new double * [3];
	surf_geometry_pts = new double * [7];
	surf_particle_mu = new double * [Maxparticle];
	surf_flow = new double * [3];
	for (int i = 0; i < 3; i++)
	{
		surf_damu[i] = new double [FO_length];
		surf_Bn_muS_muB[i] = new double [FO_length];
		surf_flow[i] = new double [FO_length];
	}
	for (int i = 0; i < 7; i++)
		surf_geometry_pts[i] = new double [FO_length];
	for (int i = 0; i < 8; i++)
		surf_pimunu[i] = new double [FO_length];
	for (int i = 0; i < Maxparticle; i++)
		surf_particle_mu[i] = new double [FO_length];
	
	for (int isurf = 0; isurf < FO_length; isurf++)
	{
		FO_surf* surf = &FOsurf_ptr[isurf];
		
		surf_particle_mu[particle_id][isurf] = surf->particle_mu[particle_id];
		
		surf_damu[0][isurf] = surf->da0;
		surf_damu[1][isurf] = surf->da1;
		surf_damu[2][isurf] = surf->da2;
		
		surf_pimunu[0][isurf] = surf->pi00;
		surf_pimunu[1][isurf] = surf->pi01;
		surf_pimunu[2][isurf] = surf->pi02;
		surf_pimunu[3][isurf] = surf->pi11;
		surf_pimunu[4][isurf] = surf->pi12;
		surf_pimunu[5][isurf] = surf->pi22;
		surf_pimunu[6][isurf] = surf->pi33;
		surf_pimunu[7][isurf] = surf->bulkPi;
		
		surf_geometry_pts[0][isurf] = surf->tau;
		surf_geometry_pts[1][isurf] = surf->r;
		surf_geometry_pts[2][isurf] = surf->phi;
		surf_geometry_pts[3][isurf] = surf->xpt;
		surf_geometry_pts[4][isurf] = surf->ypt;
		surf_geometry_pts[5][isurf] = surf->sin_phi;
		surf_geometry_pts[6][isurf] = surf->cos_phi;
		
		surf_flow[0][isurf] = surf->vx;
		surf_flow[1][isurf] = surf->vy;
		surf_flow[2][isurf] = surf->gammaT;
	}*/
	
	return;
}

void SourceVariances::Get_current_decay_string(int dc_idx, string * decay_string)
{
	// N.B. - dc_idx == 0 is thermal pion(+)s in calling loop, dc_idx > 0 gives resonance decays
	//      ==> need to use dc_idx - 1 here
	*decay_string = decay_channels.resonance_name[dc_idx - 1] + " --->> ";
	int temp_monval, tempID;
	for (int decay_part_idx = 0; decay_part_idx < decay_channels.nbody[dc_idx - 1]; decay_part_idx++)
	{
		temp_monval = decay_channels.resonance_decay_monvals[dc_idx - 1][decay_part_idx];
		//if (VERBOSE > 0) *global_out_stream_ptr << "Get_current_decay_string(): temp_monval = " << temp_monval << endl;
		if (temp_monval == 0)
			continue;
		else
		{
			tempID = lookup_particle_id_from_monval(all_particles, Nparticle, temp_monval);
			*decay_string += all_particles[tempID].name;
			if (decay_part_idx < decay_channels.nbody[dc_idx - 1] - 1) *decay_string += " + ";
		}
	}
	return;
}

int SourceVariances::lookup_resonance_idx_from_particle_id(int pid)
{
	// pid - particle index in all_particles array
	// looks up location in chosen_resonances of given value particle_id
	int result = -1;

	for (int ii = 0; ii < (int)chosen_resonances.size(); ii++)
	{
		if (chosen_resonances[ii] == pid)
		{
			result = ii;
			break;
		}
	}

	// if pid is not one of the chosen_resonances, is not the target daughter (pion(+)), is not stable and has a non-zero effective branching ratio
	if (result < 0 && pid != particle_id && all_particles[pid].stable == 0 && all_particles[pid].effective_branchratio >= 1.e-12)
	{
		*global_out_stream_ptr << " *** lookup_resonance_idx_from_particle_id(): Particle_id = " << pid
					<< " (" << all_particles[pid].name <<") not found in chosen_resonances!" << endl
					<< " *** br = " << all_particles[pid].effective_branchratio << endl;
	}
	return (result);
}

double SourceVariances::lin_int(double x1, double x2, double f1, double f2, double x)
{
	double aa, bb, cc;
	double eps=1e-100;

	if (x2 == x1) 
		aa = 0.0;
	else
		aa =(f2-f1)/(x2-x1);
	bb = f1 - aa * x1;
	cc = aa*x + bb;
	
	return cc;
}

double SourceVariances::Edndp3(double ptr, double phirin, double yr, int local_pid)
{
	return Edndp3_original(ptr, phirin, yr, local_pid);
	//if (ptr > 1.0 || phirin > SPinterp2_pphi[n_interp2_pphi_pts-1] || phirin < SPinterp2_pphi[0])
	//	return Edndp3_original(ptr, phirin, yr, local_pid);
	//else
	//	return interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[local_pid][0],
	//				ptr, phirin, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING, true);
}

double SourceVariances::Edndp3_original(double ptr, double phirin, double yr, int local_pid)
{
	double phir, val;
	double f1, f2;
	int pn, npt, nphi;
	int PTCHANGE = 1.0;

	//if(phirin < 0.0){
	//	printf("ERROR: phir %15.8le < 0 !!! \n", phirin);exit(0);}
	//if(phirin > 2.0*PI){
	//	printf("ERROR: phir %15.8le > 2PI !!! \n", phirin);exit(0);}

	phir = phirin;

	pn = local_pid;

	npt = 1;
	while((ptr > SPinterp2_pT[npt]) &&
		(npt<(n_interp2_pT_pts - 1))) npt++;

	int nphi_max = n_interp2_pphi_pts-1;
	if(phir < SPinterp2_pphi[0])
	{
		f1 = lin_int(SPinterp2_pphi[nphi_max]-2.*M_PI, SPinterp2_pphi[0],
			dN_dypTdpTdphi_moments[pn][0][npt-1][nphi_max],
			dN_dypTdpTdphi_moments[pn][0][npt-1][0], phir);
		f2 = lin_int(SPinterp2_pphi[nphi_max]-2.*M_PI, SPinterp2_pphi[0],
			dN_dypTdpTdphi_moments[pn][0][npt][nphi_max],
			dN_dypTdpTdphi_moments[pn][0][npt][0], phir);
	}
	else if(phir > SPinterp2_pphi[nphi_max])
	{
		f1 = lin_int(SPinterp2_pphi[nphi_max], SPinterp2_pphi[0]+2.*M_PI,
			dN_dypTdpTdphi_moments[pn][0][npt-1][nphi_max],
			dN_dypTdpTdphi_moments[pn][0][npt-1][0], phir);
		f2 = lin_int(SPinterp2_pphi[nphi_max], SPinterp2_pphi[0]+2.*M_PI,
			dN_dypTdpTdphi_moments[pn][0][npt][nphi_max],
			dN_dypTdpTdphi_moments[pn][0][npt][0], phir);
	}
	else
	{
		nphi = 1;
		while((phir > SPinterp2_pphi[nphi])&&(nphi < nphi_max)) nphi++;
		 /* phi interpolation */
		double phi0 = SPinterp2_pphi[nphi-1];
		double phi1 = SPinterp2_pphi[nphi];
		//cerr << "nphi = " << nphi << " and npt = " << npt << endl;
		double spec0 = dN_dypTdpTdphi_moments[pn][0][npt-1][nphi-1];
		double spec1 = dN_dypTdpTdphi_moments[pn][0][npt-1][nphi];
		double spec2 = dN_dypTdpTdphi_moments[pn][0][npt][nphi-1];
		double spec3 = dN_dypTdpTdphi_moments[pn][0][npt][nphi];
		f1 = lin_int(phi0, phi1, spec0, spec1, phir);
		f2 = lin_int(phi0, phi1, spec2, spec3, phir);

	}
	if (f1 < 0 || f2 < 0)
		printf("In Edndp3(): f1 = %f, f2 = %f, ptr = %f, phir = %f, resnum = %i", f1, f2, ptr, phirin, local_pid);
	if(f1 < 0) f1 = 0.0;
	if(f2 < 0) f2 = 0.0;
	f1 = f1 + 1e-100;
	f2 = f2 + 1e-100;
	if(ptr > PTCHANGE)
	{
		f1 = log(f1);
		f2 = log(f2);
		//val = lin_int(SPinterp2_pT[npt-1], SPinterp2_pT[npt],
		//	f1, f2, ptr);
	}
	//else
	//{
	//	val = 
	//}
	val = lin_int(SPinterp2_pT[npt-1], SPinterp2_pT[npt],
		f1, f2, ptr);
         /*if(isnan(val))
         {
             printf("dEdndy3 val \n");
             printf("f1, %15.8lf f2, %15.8lf  \n", f1, f2);
             printf(" nphi  %i npt %i \n", nphi,npt);
             printf(" f1  %15.8le %15.8le  \n", f1, f2);
             printf(" phi  %15.8lf %15.8lf  \n", PHI[nphi-1], PHI[nphi]);
             printf(" pt   %15.8lf %15.8lf  \n", particle[pn].pt[npt-1],particle[pn].pt[npt]);
             printf(" phi  %15.8lf pt %15.8lf    val %15.8lf \n", phir, ptr,val);
             printf(" phi %15.8le %15.8le \n",dN_dypTdpTdphi_moments[pn][0][npt][nphi-1],
                                         dN_dypTdpTdphi_moments[pn][0][npt][nphi]);
             printf(" pt  %15.8le %15.8le \n",dN_dypTdpTdphi_moments[pn][0][npt-1][nphi-1],
                                     dN_dypTdpTdphi_moments[pn][0][npt-1][nphi]);
             exit(-1);
         }*/
	if(ptr > PTCHANGE)
		val = exp(val);

	//if (1) *global_out_stream_ptr << "In Edndp3(): {ptr, phir, local_pid, val} = " << ptr << "   " << phir << "   " << local_pid << "   " << val << endl;

	return val;
}

//End of file

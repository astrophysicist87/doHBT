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
#include "Stopwatch.h"
#include "gauss_quadrature.h"

using namespace std;

// only need to calculated interpolation grid of spacetime moments for each resonance, NOT each decay channel!
bool recycle_previous_moments = false;
bool recycle_similar_moments = false;
int reso_particle_id_of_moments_to_recycle = -1;
string reso_name_of_moments_to_recycle = "NULL";
string current_decay_channel_string = "NULL";

template < typename T >
void check_for_NaNs(string variable_name, const T variable_value, ofstream& localout)
{
	if (isnan(variable_value))
		localout << "ERROR: " << variable_name << " = " << variable_value << endl;
	return;
}

//need to define some variables for quick evaluation of S_direct
double Sdir_Y = 0.0, Sdir_R = 5., Sdir_Deltau = 1., Sdir_Deleta = 1.2;
double Sdir_eta0 = 0.0, Sdir_tau0 = 5., Sdir_etaf = 0.0, Sdir_T = 0.15;
double Sdir_prefactor, Sdir_rt, Sdir_H, Sdir_etat, Sdir_ch_Y_m_eta, Sdir_expon;
double Sdir_term1, Sdir_term2, Sdir_term3;

double SourceVariances::S_direct(double r, double eta, double tau, double MT, double PT, double cos_phi_m_Phi)
{
	//Sdir_prefactor = (2.*Sdir_Jr + 1.)/(twopi*twopi*twopi);
	Sdir_rt = r/Sdir_R;
	Sdir_term1 = 0.5*Sdir_rt*Sdir_rt;
	Sdir_term2 = 0.0*0.5*(eta-Sdir_eta0)*(eta-Sdir_eta0)/(Sdir_Deleta*Sdir_Deleta);
	Sdir_term3 = 0.5*(tau-Sdir_tau0)*(tau-Sdir_tau0)/(Sdir_Deltau*Sdir_Deltau);
	//Sdir_term2 = 0.0;
	//Sdir_term3 = 0.0;
	Sdir_H = exp( - Sdir_term1 - Sdir_term2 - Sdir_term3 ) / (M_PI * Sdir_Deltau);
	//Sdir_H = exp( - Sdir_term1 - Sdir_term2 - Sdir_term3 );
	Sdir_ch_Y_m_eta = cosh(Sdir_Y - eta);
	Sdir_expon = -(MT/Sdir_T)*Sdir_ch_Y_m_eta*cosh(Sdir_etaf*Sdir_rt)
			+ (PT/Sdir_T)*sinh(Sdir_etaf*Sdir_rt)*cos_phi_m_Phi;
	return (r * MT * Sdir_ch_Y_m_eta * Sdir_H * exp(Sdir_expon));
	//return (1.);
	//return((r > Sdir_R) ? 0. : 1.);
}

double SourceVariances::place_in_range(double phi, double min, double max)
{
	while (phi < min || phi > max)
	{
		if (phi < min) phi += twopi;
		else phi -= twopi;
	}

	return (phi);
}

// ************************************************************
// NEW AND IMPROVED version of Analyze_sourcefunction()
// ************************************************************
void SourceVariances::Analyze_sourcefunction(FO_surf* FOsurf_ptr)
{
	Stopwatch BIGsw;
	*global_out_stream_ptr << "Plane angle calculations..." << endl;
	BIGsw.tic();
	double plane_psi = 0.0;
	int iorder = USE_PLANE_PSI_ORDER;
	if (USE_PLANE_PSI_ORDER)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "Determine nth-order plane angles..." << endl;
		Determine_plane_angle(current_FOsurf_ptr, 0, true);	//uses only thermal pions...
		if (VERBOSE > 0) *global_out_stream_ptr << "Analyzing source function w.r.t. " << iorder << " th-order participant plane angle..." << endl;
		if (VERBOSE > 0) *global_out_stream_ptr << "psi = " << plane_psi << endl;
		plane_psi = plane_angle[iorder];
	}
	else
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "Analyzing source function w.r.t. psi_0 = " << plane_psi << endl;
	}
	global_plane_psi = plane_psi;
	BIGsw.toc();
	*global_out_stream_ptr << "\t ...finished in " << BIGsw.takeTime() << " seconds." << endl;

	int decay_channel_loop_cutoff = n_decay_channels;			//loop over direct pions and decay_channels

	if (read_in_all_dN_dypTdpTdphi)	//read in spectra if already calculated
	{
		//if (VERBOSE > 0) *global_out_stream_ptr << "currentfolderindex = " << currentfolderindex << endl;
		Read_in_all_dN_dypTdpTdphi(currentfolderindex);
		if (VERBOSE > 0) *global_out_stream_ptr << "************************************************************" << endl
												<< "* Read in all (thermal) space-time moments!" << endl
												<< "************************************************************" << endl << endl;

	}
	else	// calculate necessary spectra from scratch
	{
		*global_out_stream_ptr << "Setting spacetime moments grid..." << endl;
		BIGsw.tic();
		// ************************************************************
		// loop over decay_channels (idc == 0 corresponds to thermal pions)
		// ************************************************************
		for (int idc = 0; idc <= decay_channel_loop_cutoff; idc++)				//this is inefficient, but will do the job for now
		{
			// ************************************************************
			// check whether to do this decay channel
			// ************************************************************
			if (idc > 0 && thermal_pions_only)
				break;
			else if (!Do_this_decay_channel(idc))
				continue;
	
//temporary: skip omega and rho(0)
//if (idc > 1 && decay_channels.resonance_particle_id[idc-1] != 21)
//	continue;

			// ************************************************************
			// if so, set decay channel info
			// ************************************************************
			Set_current_particle_info(idc);
	
			// ************************************************************
			// decide whether to recycle old moments or calculate new moments
			// ************************************************************
			Get_spacetime_moments(FOsurf_ptr, idc);
		}	//computing all resonances' spacetime moments here first
			//THEN do phase-space integrals
	
		if (VERBOSE > 0) *global_out_stream_ptr << endl << "************************************************************"
												<< endl << "* Computed all (thermal) space-time moments!" << endl
												<< "************************************************************" << endl << endl;
		BIGsw.toc();
		*global_out_stream_ptr << "\t ...finished all (thermal) space-time moments in " << BIGsw.takeTime() << " seconds." << endl;	}
	if (output_all_dN_dypTdpTdphi)
	{
		Output_all_dN_dypTdpTdphi(currentfolderindex);
		if (VERBOSE > 0) *global_out_stream_ptr << endl << "************************************************************"
												<< endl << "* Output all (thermal) space-time moments!" << endl
												<< "************************************************************" << endl << endl;
	}

	if (SPACETIME_MOMENTS_ONLY)
		return;


	*global_out_stream_ptr << "Computing all phase-space integrals..." << endl;
	double current_dNd3p_00 = 249.013183, thermal_dNd3p_00 = 249.013183, previous_dNd3p_00 = 249.013183;
	BIGsw.tic();
	// ************************************************************
	// Compute feeddown with heaviest resonances first
	// ************************************************************
	for (int idc = 1; idc <= decay_channel_loop_cutoff; idc++)
	{
		// ************************************************************
		// check whether to do this decay channel
		// ************************************************************
		if (thermal_pions_only)
			break;
		else if (!Do_this_decay_channel(idc))
			continue;

//temporary: skip omega and rho(0)
//if (decay_channels.resonance_particle_id[idc-1] == 12 || decay_channels.resonance_particle_id[idc-1] == 10)
//	continue;

		// ************************************************************
		// if so, set decay channel info
		// ************************************************************
		Set_current_particle_info(idc);

		// ************************************************************
		// begin source variances calculations here...
		// ************************************************************
		Allocate_decay_channel_info();				// allocate needed memory
		for (int idc_DI = 0; idc_DI < current_reso_nbody; idc_DI++)
		{
			int daughter_resonance_particle_id = -1;
			if (!Do_this_daughter_particle(idc, idc_DI, &daughter_resonance_particle_id))
				continue;
			Set_current_daughter_info(idc, idc_DI);
			if (daughter_resonance_particle_id == target_particle_id)
				previous_dNd3p_00 = current_dNd3p_00;
			Do_resonance_integrals(current_resonance_particle_id, daughter_resonance_particle_id, idc);
			if (daughter_resonance_particle_id == target_particle_id)
			{
				current_dNd3p_00 = dN_dypTdpTdphi_moments[daughter_resonance_particle_id][0][0][0];
				*global_out_stream_ptr << all_particles[daughter_resonance_particle_id].name << " from " << all_particles[current_resonance_particle_id].name
										<< " in decay channel " << current_decay_channel_string << ": "
										<< current_dNd3p_00 << "   " << current_dNd3p_00-thermal_dNd3p_00 << "   " << current_dNd3p_00-previous_dNd3p_00 << endl;
			}
			//if (VERBOSE > 0) *global_out_stream_ptr << endl;
		}
		Delete_decay_channel_info();				// free up memory
	}											// END of decay channel loop

	// ************************************************************
	// now get HBT radii from source integrals and Fourier transform
	// ************************************************************
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		if (abs(K_T[iKT]) < 1.e-10)
			continue;

		double m_perp = sqrt(K_T[iKT]*K_T[iKT] + particle_mass*particle_mass);
		beta_perp = K_T[iKT]/(m_perp*cosh(K_y));

		for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			Compute_source_variances(iKT, iKphi);
			if (INCLUDE_SOURCE_VARIANCES)
			{
			Calculate_R2_side(iKT, iKphi);
			Calculate_R2_out(iKT, iKphi);
			Calculate_R2_long(iKT, iKphi);
			Calculate_R2_outside(iKT, iKphi);
			Calculate_R2_sidelong(iKT, iKphi);
			Calculate_R2_outlong(iKT, iKphi);
			}
		}

		R2_Fourier_transform(iKT, plane_psi);
	}

   return;
}

bool SourceVariances::Do_this_decay_channel(int dc_idx)
{
	string local_name = "Thermal pion(+)";
	if (dc_idx == 0)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << endl << local_name << ": doing this one." << endl;
		return true;
	}
	else
	{
		local_name = decay_channels.resonance_name[dc_idx-1];
		Get_current_decay_string(dc_idx, &current_decay_channel_string);
	}
	if (DO_ALL_DECAY_CHANNELS)
		return true;
	else if (decay_channels.include_channel[dc_idx-1])
	{
		;
	}
	else
	{
		if (VERBOSE > 0) *global_out_stream_ptr << endl << local_name << ": skipping decay " << current_decay_channel_string << "." << endl;
	}

	return (decay_channels.include_channel[dc_idx-1]);
}

// ************************************************************
// Checks whether to do daughter particle for given decay channel
// ************************************************************
bool SourceVariances::Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid)
{
	// assume dc_idx > 0
	string local_name = decay_channels.resonance_name[dc_idx - 1];

	// look up daughter particle info
	int temp_monval = decay_channels.resonance_decay_monvals[dc_idx - 1][daughter_idx];

	if (temp_monval == 0)
		return false;

	int temp_ID = lookup_particle_id_from_monval(all_particles, Nparticle, temp_monval);
	//*daughter_resonance_idx = lookup_resonance_idx_from_particle_id(temp_ID) + 1;
	*daughter_resonance_pid = temp_ID;
	// if daughter was found in chosen_resonances or is pion(+), this is correct
	particle_info temp_daughter = all_particles[temp_ID];

	if (*daughter_resonance_pid < 0 && temp_daughter.monval != particle_monval && temp_daughter.effective_branchratio >= 1.e-12)
		*global_out_stream_ptr << "Couldn't find " << temp_daughter.name << " in chosen_resonances!  Results are probably not reliable..." << endl;

	//bool daughter_does_not_contribute = ( (temp_daughter.stable == 1 || temp_daughter.effective_branchratio < 1.e-12) && temp_daughter.monval != particle_monval );
	bool daughter_does_not_contribute = ( (temp_daughter.decays_Npart[0] == 1 || temp_daughter.effective_branchratio < 1.e-12) && temp_daughter.monval != particle_monval );

	// if daughter particle gives no contribution to final pion spectra
	if (daughter_does_not_contribute)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << local_name << ": in decay " << current_decay_channel_string << ", skipping " << temp_daughter.name
												<< " (daughter_resonance_pid = " << *daughter_resonance_pid << ")." << endl;
		return false;
	}
	else
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << local_name << ": in decay " << current_decay_channel_string << ", doing " << temp_daughter.name
												<< " (daughter_resonance_pid = " << *daughter_resonance_pid << ")." << endl;
		return true;
	}
}

void SourceVariances::Set_current_particle_info(int dc_idx)
{
	if (dc_idx == 0)
	{
		muRES = particle_mu;
		signRES = particle_sign;
		gRES = particle_gspin;
		current_resonance_particle_id = target_particle_id;
		
		return;
	}
	else
	{
		// assume dc_idx > 0
		string local_name = decay_channels.resonance_name[dc_idx - 1];

		if (VERBOSE > 0) *global_out_stream_ptr << local_name << ": doing decay " << current_decay_channel_string << "." << endl
			<< "\t * " << local_name << ": setting information for this decay channel..." << endl;

		if (dc_idx > 1)
		{
			//cerr << "Setting previous decay channel information for dc_idx = " << dc_idx << endl;
			previous_resonance_particle_id = current_resonance_particle_id;		//for look-up in all_particles
			previous_decay_channel_idx = current_decay_channel_idx;			//different for each decay channel
			previous_resonance_idx = current_resonance_idx;				//different for each decay channel
			previous_resonance_mass = current_resonance_mass;
			previous_resonance_Gamma = current_resonance_Gamma;
			previous_resonance_total_br = current_resonance_total_br;
			previous_resonance_direct_br = current_resonance_direct_br;
			previous_reso_nbody = current_reso_nbody;
		}
		//cerr << "Setting current decay channel information for dc_idx = " << dc_idx << endl;
		current_decay_channel_idx = dc_idx;
		current_resonance_idx = decay_channels.resonance_idx[dc_idx-1];
		current_resonance_particle_id = decay_channels.resonance_particle_id[dc_idx-1];
		current_resonance_mass = decay_channels.resonance_mass[dc_idx-1];
		current_resonance_Gamma = decay_channels.resonance_Gamma[dc_idx-1];
		current_resonance_total_br = decay_channels.resonance_total_br[dc_idx-1];
		current_resonance_direct_br = decay_channels.resonance_direct_br[dc_idx-1];
		current_reso_nbody = decay_channels.nbody[dc_idx-1];
		
		// might want to rename these for notational consistency...
		muRES = decay_channels.resonance_mu[dc_idx-1];
		signRES = decay_channels.resonance_sign[dc_idx-1];
		gRES = decay_channels.resonance_gspin[dc_idx-1];
		
		if (dc_idx > 1)
		{
			int similar_particle_idx = -1;
			int temp_reso_idx = decay_channels.resonance_idx[dc_idx-1];
			
			if ( current_resonance_particle_id == previous_resonance_particle_id )
			{
				//previous resonance is the same as this one...
				recycle_previous_moments = true;
				recycle_similar_moments = false;
			}
			else if ( Search_for_similar_particle( temp_reso_idx, &similar_particle_idx ) )
			{
				//previous resonance is NOT the same as this one BUT this one is sufficiently similar to some preceding one...
				recycle_previous_moments = false;
				recycle_similar_moments = true;
				reso_particle_id_of_moments_to_recycle = chosen_resonances[similar_particle_idx];
			}
			else
			{
				recycle_previous_moments = false;
				recycle_similar_moments = false;
				reso_particle_id_of_moments_to_recycle = -1;	//guarantees it won't be used spuriously
			}
		}
	}
	
	return;
}

void SourceVariances::Set_current_daughter_info(int dc_idx, int daughter_idx)
{
	if (dc_idx > 1)
	{
		previous_resonance_particle_id = current_resonance_particle_id;		//for look-up in all_particles
		previous_decay_channel_idx = current_decay_channel_idx;			//different for each decay channel
		previous_resonance_idx = current_resonance_idx;
		previous_resonance_mass = current_resonance_mass;
		previous_resonance_Gamma = current_resonance_Gamma;
		previous_m2_Gamma = current_m2_Gamma;
		previous_m3_Gamma = current_m3_Gamma;
		previous_resonance_total_br = current_resonance_total_br;
		previous_resonance_direct_br = current_resonance_direct_br;
		previous_reso_nbody = current_reso_nbody;
		previous_daughter_mass = current_daughter_mass;
		previous_daughter_Gamma = current_daughter_Gamma;
	}
	current_decay_channel_idx = dc_idx;
	current_resonance_idx = decay_channels.resonance_idx[dc_idx-1];
	current_resonance_particle_id = decay_channels.resonance_particle_id[dc_idx-1];
	current_resonance_mass = decay_channels.resonance_mass[dc_idx-1];
	current_resonance_Gamma = decay_channels.resonance_Gamma[dc_idx-1];
	current_resonance_total_br = decay_channels.resonance_total_br[dc_idx-1];
	current_resonance_direct_br = decay_channels.resonance_direct_br[dc_idx-1];
	current_reso_nbody = decay_channels.nbody[dc_idx-1];
	current_daughter_mass = decay_channels.resonance_decay_masses[dc_idx-1][daughter_idx];
	current_daughter_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][daughter_idx];

	// might want to rename these for notational consistency...
	muRES = decay_channels.resonance_mu[dc_idx-1];
	signRES = decay_channels.resonance_sign[dc_idx-1];
	gRES = decay_channels.resonance_gspin[dc_idx-1];

	// set non-daughter decay masses for computing contributions to spectra of daughter
	double m2ex = 0.0, m3ex = 0.0, m4ex = 0.0;
	switch(current_reso_nbody)
	{
		case 2:
			current_resonance_decay_masses[1] = 0.0;
			current_m3_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][0];
			if (daughter_idx == 0)
			{
				current_resonance_decay_masses[0] = decay_channels.resonance_decay_masses[dc_idx-1][1];
				current_m2_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][1];
			}
			else
			{
				current_resonance_decay_masses[0] = decay_channels.resonance_decay_masses[dc_idx-1][0];
				current_m2_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][0];
			}
			break;
		case 3:
			if (daughter_idx == 0)
			{
				current_resonance_decay_masses[0] = decay_channels.resonance_decay_masses[dc_idx-1][1];
				current_resonance_decay_masses[1] = decay_channels.resonance_decay_masses[dc_idx-1][2];
				current_m2_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][1];
				current_m3_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][2];
			}
			else if (daughter_idx == 1)
			{
				current_resonance_decay_masses[0] = decay_channels.resonance_decay_masses[dc_idx-1][0];
				current_resonance_decay_masses[1] = decay_channels.resonance_decay_masses[dc_idx-1][2];
				current_m2_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][0];
				current_m3_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][2];
			}
			else
			{
				current_resonance_decay_masses[0] = decay_channels.resonance_decay_masses[dc_idx-1][0];
				current_resonance_decay_masses[1] = decay_channels.resonance_decay_masses[dc_idx-1][1];
				current_m2_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][0];
				current_m3_Gamma = decay_channels.resonance_decay_Gammas[dc_idx-1][1];
			}
			break;
		case 4:
			if (daughter_idx == 0)
			{
				m2ex = decay_channels.resonance_decay_masses[dc_idx-1][1];
				m3ex = decay_channels.resonance_decay_masses[dc_idx-1][2];
				m4ex = decay_channels.resonance_decay_masses[dc_idx-1][3];
			}
			else if (daughter_idx == 1)
			{
				m2ex = decay_channels.resonance_decay_masses[dc_idx-1][0];
				m3ex = decay_channels.resonance_decay_masses[dc_idx-1][2];
				m4ex = decay_channels.resonance_decay_masses[dc_idx-1][3];
			}
			else if (daughter_idx == 2)
			{
				m2ex = decay_channels.resonance_decay_masses[dc_idx-1][0];
				m3ex = decay_channels.resonance_decay_masses[dc_idx-1][1];
				m4ex = decay_channels.resonance_decay_masses[dc_idx-1][3];
			}
			else
			{
				m2ex = decay_channels.resonance_decay_masses[dc_idx-1][0];
				m3ex = decay_channels.resonance_decay_masses[dc_idx-1][1];
				m4ex = decay_channels.resonance_decay_masses[dc_idx-1][2];
			}
			current_resonance_decay_masses[0] = m2ex;
			current_resonance_decay_masses[1] = 0.5 * (m3ex + m4ex + current_resonance_mass - current_daughter_mass - m2ex);
			// approximation obtained from earlier resonances code
			//*global_out_stream_ptr << "Current decay " << current_decay_channel_string << ", br = " << current_resonance_direct_br
			//						<< ": {m2ex, m3ex, m4ex, m3eff} = {"
			//						<< m2ex << ", " << m3ex << ", " << m4ex << ", " << current_resonance_decay_masses[1] << "}" << endl;
			break;
		default:
			cerr << "Set_current_daughter_info(): shouldn't have ended up here, bad value of current_reso_nbody!" << endl;
			exit(1);
	}
}

bool SourceVariances::Search_for_similar_particle(int reso_idx, int * result)
{
	// for the timebeing, just search from beginning of decay_channels until similar particle is found;
	// should be more careful, since could lead to small numerical discrepancies if similar particle was
	// already recycled by some other (dissimilar) particle, but ignore this complication for now...
	*result = -1;
	
	for (int local_ir = 0; local_ir < reso_idx; local_ir++)
	{// only need to search decay_channels that have already been calculated
		if (particles_are_the_same(local_ir, reso_idx))
		{
			*result = local_ir;
			break;
		}
	}
	
	return (*result >= 0);
}

//**********************************************************************************************
bool SourceVariances::particles_are_the_same(int reso_idx1, int reso_idx2)
{
	int icr1 = chosen_resonances[reso_idx1];
	int icr2 = chosen_resonances[reso_idx2];
	//int icr1 = reso_idx1;
	//int icr2 = reso_idx2;
	if (all_particles[icr1].sign != all_particles[icr2].sign)
		return false;
	if (abs(all_particles[icr1].mass-all_particles[icr2].mass) / (all_particles[icr2].mass+1.e-30) > PARTICLE_DIFF_TOLERANCE)
		return false;
	//assume chemical potential mu is constant over entire FO surface
	double chem1 = all_particles[icr1].mu, chem2 = all_particles[icr2].mu;
	if (2.*abs(chem1 - chem2)/(chem1 + chem2 + 1.e-30) > PARTICLE_DIFF_TOLERANCE)
		return false;

	return true;
}

void SourceVariances::Recycle_spacetime_moments()
{
	for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
	for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
	for(int wfi = 0; wfi < n_weighting_functions; wfi++)
		dN_dypTdpTdphi_moments[current_resonance_particle_id][wfi][ipt][iphi] = dN_dypTdpTdphi_moments[reso_particle_id_of_moments_to_recycle][wfi][ipt][iphi];
		//N.B. - not dc_idx - 1, since spacetime moments for dc_idx = 0 are just thermal pions

	return;
}

void SourceVariances::Get_spacetime_moments(FO_surf* FOsurf_ptr, int dc_idx)
{
//**************************************************************
//Set resonance name
//**************************************************************
//debugger(__LINE__,__FILE__);
	string local_name = "Thermal pion(+)";
	if (dc_idx > 0)
		local_name = decay_channels.resonance_name[dc_idx-1];
//**************************************************************
//Decide what to do with this resonance / decay channel
//**************************************************************
	if (recycle_previous_moments && RECYCLE_ST_MOMENTS && dc_idx > 1)	// similar (but different) earlier resonance
	{
		if (VERBOSE > 0) *global_out_stream_ptr << local_name
			<< ": new parent resonance (" << decay_channels.resonance_name[current_decay_channel_idx-1] << ", dc_idx = " << current_decay_channel_idx
			<< ") same as preceding parent resonance \n\t\t--> reusing old dN_dypTdpTdphi_moments!" << endl;
		//Recycle_spacetime_moments();
	}
	else if (recycle_similar_moments && RECYCLE_ST_MOMENTS && dc_idx > 1)	// similar (but different) earlier resonance
	{
		if (VERBOSE > 0) *global_out_stream_ptr << local_name
			<< ": new parent resonance (" << decay_channels.resonance_name[current_decay_channel_idx-1] << ", dc_idx = " << current_decay_channel_idx
			<< ") sufficiently close to preceding parent resonance (" << all_particles[reso_particle_id_of_moments_to_recycle].name
			<< ", reso_particle_id = " << reso_particle_id_of_moments_to_recycle << ") \n\t\t--> reusing old dN_dypTdpTdphi_moments!" << endl;
		Recycle_spacetime_moments();
	}
	else
	{
		if (dc_idx == 0)	//if it's thermal pions
		{
			if (VERBOSE > 0) *global_out_stream_ptr << "  --> Computing dN_dypTdpTdphi_moments for thermal pion(+)!" << endl;
		}
		else if (dc_idx == 1)	//if it's the first resonance
		{
			if (VERBOSE > 0) *global_out_stream_ptr << "  --> Computing dN_dypTdpTdphi_moments for " << local_name << endl;
		}
		else			//if it's a later resonance
		{
			if (VERBOSE > 0 && !RECYCLE_ST_MOMENTS) *global_out_stream_ptr << local_name
				<< ": RECYCLE_ST_MOMENTS currently set to false \n\t\t--> calculating new dN_dypTdpTdphi_moments!" << endl;
			else if (VERBOSE > 0 && !recycle_previous_moments && !recycle_similar_moments) *global_out_stream_ptr << local_name
				<< ": new parent resonance (" << decay_channels.resonance_name[current_decay_channel_idx-1] << ", dc_idx = " << current_decay_channel_idx
				<< ") dissimilar from all preceding decay_channels \n\t\t--> calculating new dN_dypTdpTdphi_moments!" << endl;
			else
			{
				cerr << "You shouldn't have ended up here!" << endl;
				exit(1);
			}
		}
//debugger(__LINE__,__FILE__);
		Set_dN_dypTdpTdphi_moments(FOsurf_ptr, current_resonance_particle_id);
	}
//**************************************************************
//Spacetime moments now set
//**************************************************************
//debugger(__LINE__,__FILE__);
	return;
}

void SourceVariances::Set_dN_dypTdpTdphi_moments(FO_surf* FOsurf_ptr, int local_pid)
{
	time_t starttime, stoptime;
   	struct tm * timeinfo;

	double localmass = all_particles[local_pid].mass;
	string local_name = "Thermal pion(+)";
	if (local_pid != target_particle_id)
		local_name = all_particles[local_pid].name;

	for(int i=0; i<eta_s_npts; i++)
	{
		double local_eta_s = eta_s[i];
		double local_cosh = cosh(SP_p_y - local_eta_s);
		double local_sinh = sinh(SP_p_y - local_eta_s);

		for(int ipt=0; ipt<n_interp2_pT_pts; ipt++)
		{
			double mT = sqrt(localmass*localmass + SPinterp2_pT[ipt]*SPinterp2_pT[ipt]);
			SPinterp2_p0[ipt][i] = mT*local_cosh;
			SPinterp2_pz[ipt][i] = mT*local_sinh;
		}
	}
//debugger(__LINE__,__FILE__);
	//**************************************************************
	//Timing checks
	//**************************************************************
	time (&starttime);
 	timeinfo = localtime (&starttime);
 	if (VERBOSE > 0) *global_out_stream_ptr << local_name << ":" << endl << "   * Started setting spacetime moments at " << asctime(timeinfo);
	//**************************************************************
	Cal_dN_dypTdpTdphi_with_weights_polar(FOsurf_ptr, local_pid);
	//**************************************************************
	time (&stoptime);
 	timeinfo = localtime (&stoptime);
 	if (VERBOSE > 0) *global_out_stream_ptr << "   * Finished setting spacetime moments at " << asctime(timeinfo)
						<< "   * Took " << difftime(stoptime, starttime) << " seconds." << endl;
	//**************************************************************


	return;
}

void SourceVariances::Cal_dN_dypTdpTdphi_with_weights_polar(FO_surf* FOsurf_ptr, int local_pid)
{
cout << "local_pid = " << local_pid << endl;
//debugger(__LINE__,__FILE__);
	//int local_reso_idx;
	double z0, z1, z2, z3;
	double sign, degen, localmass, mu;
	if (local_pid == target_particle_id)
	{
//debugger(__LINE__,__FILE__);
		sign = particle_sign;
		degen = particle_gspin;
		localmass = particle_mass;
		mu = FOsurf_ptr[0].particle_mu[particle_id];
	}
	else
	{
//debugger(__LINE__,__FILE__);
		sign = all_particles[local_pid].sign;
		degen = all_particles[local_pid].gspin;
		localmass = all_particles[local_pid].mass;
		mu = all_particles[local_pid].mu;
	}
//debugger(__LINE__,__FILE__);
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	Tdec = (&FOsurf_ptr[0])->Tdec;
	double one_by_Tdec = 1./Tdec;
	Pdec = (&FOsurf_ptr[0])->Pdec;
	Edec = (&FOsurf_ptr[0])->Edec;
	double deltaf_prefactor = 0.;
	if (use_delta_f) deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));
//debugger(__LINE__,__FILE__);
	//declare variables used below here to see if this speeds up code:
	double tau, vx, vy, da0, da1, da2;
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
//debugger(__LINE__,__FILE__);
	for(int isurf=0; isurf<FO_length ; isurf++)
	{
		surf = &FOsurf_ptr[isurf];
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
//debugger(__LINE__,__FILE__);
		for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
		{
			pT = SPinterp2_pT[ipt];
//debugger(__LINE__,__FILE__);
			for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
			{
//debugger(__LINE__,__FILE__);
				sin_pphi = sin_SPinterp2_pphi[iphi];
				cos_pphi = cos_SPinterp2_pphi[iphi];
				px = pT*cos_pphi;
				py = pT*sin_pphi;
				sin_phi_m_pphi = sin_temp_phi * cos_pphi - cos_temp_phi * sin_pphi;
				cos_phi_m_pphi = cos_temp_phi * cos_pphi + sin_temp_phi * sin_pphi;
				z1 = temp_r * cos_phi_m_pphi;
				z2 = temp_r * sin_phi_m_pphi;
				for(int ieta=0; ieta < eta_s_npts; ieta++)
				{
//debugger(__LINE__,__FILE__);
					p0 = SPinterp2_p0[ipt][ieta];
					pz = SPinterp2_pz[ipt][ieta];
	
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
					deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);
					//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
					S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);
					//ignore points where delta f is large or emission function goes negative from pdsigma
					if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
						S_p = 0.0;

					S_p_withweight = S_p*tau*eta_s_weight[ieta];
					z0 = tau*ch_eta_s[ieta];
					z3 = tau*sh_eta_s[ieta];
					dN_dypTdpTdphi_moments[local_pid][0][ipt][iphi] += eta_even_factor*S_p_withweight;			//<1>
		//if (local_pid == 1 && isurf == 0) cout << "EF: " << p0 << "   " << px << "   " << py << "   " << pz << "   " << S_p*tau << endl;
					if (INCLUDE_SOURCE_VARIANCES)
					{
					dN_dypTdpTdphi_moments[local_pid][1][ipt][iphi] += eta_even_factor*S_p_withweight*z2;			//<x_s>
					dN_dypTdpTdphi_moments[local_pid][2][ipt][iphi] += eta_even_factor*S_p_withweight*z2*z2;			//<x^2_s>
					dN_dypTdpTdphi_moments[local_pid][3][ipt][iphi] += eta_even_factor*S_p_withweight*z1;			//<x_o>
					dN_dypTdpTdphi_moments[local_pid][4][ipt][iphi] += eta_even_factor*S_p_withweight*z1*z1;			//<x^2_o>
					dN_dypTdpTdphi_moments[local_pid][5][ipt][iphi] += eta_odd_factor*S_p_withweight*z3;			//<x_l>
					dN_dypTdpTdphi_moments[local_pid][6][ipt][iphi] += eta_even_factor*S_p_withweight*z3*z3;			//<x^2_l>
					dN_dypTdpTdphi_moments[local_pid][7][ipt][iphi] += eta_even_factor*S_p_withweight*z0;			//<t>
					dN_dypTdpTdphi_moments[local_pid][8][ipt][iphi] += eta_even_factor*S_p_withweight*z0*z0;			//<t^2>
					dN_dypTdpTdphi_moments[local_pid][9][ipt][iphi] += eta_even_factor*S_p_withweight*z2*z1;			//<x_s x_o>
					dN_dypTdpTdphi_moments[local_pid][10][ipt][iphi] += eta_odd_factor*S_p_withweight*z2*z3;			//<x_s x_l>
					dN_dypTdpTdphi_moments[local_pid][11][ipt][iphi] += eta_even_factor*S_p_withweight*z2*z0;		//<x_s t>
					dN_dypTdpTdphi_moments[local_pid][12][ipt][iphi] += eta_odd_factor*S_p_withweight*z1*z3;			//<x_o x_l>
					dN_dypTdpTdphi_moments[local_pid][13][ipt][iphi] += eta_even_factor*S_p_withweight*z1*z0;		//<x_o t>
					dN_dypTdpTdphi_moments[local_pid][14][ipt][iphi] += eta_odd_factor*S_p_withweight*z3*z0;			//<x_l t>
					}
				}
			}
		}
	}
//debugger(__LINE__,__FILE__);

	//only use this for spectra for now
	double temp;
	//set log of dN_dypTdpTdphi_moments...
	for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
	for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
	{
		temp = dN_dypTdpTdphi_moments[local_pid][0][ipt][iphi];
		ln_dN_dypTdpTdphi_moments[local_pid][0][ipt][iphi] = log(abs(temp));
		sign_of_dN_dypTdpTdphi_moments[local_pid][0][ipt][iphi] = sgn(temp);
	}
//debugger(__LINE__,__FILE__);

	return;
}

void SourceVariances::Cal_dN_dypTdpTdphi_interpolate_polar_grid(double* SP_pT, double* SP_pphi)
{
	for(int ipt = 0; ipt < n_SP_pT; ipt++)
	for(int iphi = 0; iphi < n_SP_pphi; iphi++)
	{
		double pT = SP_pT[ipt];
		double pphi = SP_pphi[iphi];
		SV_dN_dypTdpTdphi[ipt][iphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[0][0],
							pT, pphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		if (0) cout << "-2   " << pT << "  " << pphi << "  " << pT*cos(pphi) << "  " << pT*sin(pphi) << "  " << dN_dypTdpTdphi[ipt][iphi] << endl;
	}

   return;
}

void SourceVariances::Determine_plane_angle(FO_surf* FOsurf_ptr, int dc_idx, bool thermal_particles_only /*= false*/)
{
   double localmass = particle_mass;
   if (dc_idx > 0) return;	//don't really care about thermal resonance distributions at the moment,
				//so just skip this part
   //if (dc_idx > 0) localmass = decay_channels.resonance_mass[dc_idx-1];
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

	if (thermal_particles_only)	//no interpolation at all; calculate everything exactly
	{
		Cal_dN_dypTdpTdphi(p0, px, py, pz, FOsurf_ptr);
		if (VERBOSE > 0) *global_out_stream_ptr << "  --> No interpolation!  Calculating everything exactly..." << endl;
	}
	else						//using polar grid for interpolation (pT, pphi)
	{
   		Cal_dN_dypTdpTdphi_interpolate_polar_grid(SP_pT, SP_pphi);
		if (VERBOSE > 0) *global_out_stream_ptr << "  --> Interpolating on a polar grid..." << endl;
	}
   
	if (thermal_particles_only)
	{
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
	}
	else	// if not doing thermal particles only
	{
		for(int ipt=0; ipt<n_SP_pT; ipt++)
	   	{
			for(int iphi=0; iphi<n_SP_pphi; iphi++)
			{
				SV_dN_dydphi[iphi] += SV_dN_dypTdpTdphi[ipt][iphi]*SP_pT[ipt]*SP_pT_weight[ipt];
				SV_pTdN_dydphi[iphi] += SV_dN_dypTdpTdphi[ipt][iphi]*SP_pT[ipt]*SP_pT[ipt]*SP_pT_weight[ipt];
				SV_dN_dypTdpT[ipt] += SV_dN_dypTdpTdphi[ipt][iphi]*SP_pphi_weight[iphi];
			}
		}
	   	double norm = 0.0e0;
	   	for(int iphi=0; iphi<n_SP_pphi; iphi++)
			norm += SV_dN_dydphi[iphi]*SP_pphi_weight[iphi];
	   	for(int iorder=0; iorder < n_order; iorder++)
		{
			double cosine = 0.0e0;
			double sine = 0.0e0;
			for(int iphi=0; iphi<n_SP_pphi; iphi++)
			{
				cosine += SV_dN_dydphi[iphi]*cos(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
				sine += SV_dN_dydphi[iphi]*sin(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
				//get pT-differential v_n here
				for(int ipt=0; ipt<n_SP_pT; ipt++)
				{
					cosine_iorder[ipt][iorder] += SV_dN_dypTdpTdphi[ipt][iphi]*cos(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
					sine_iorder[ipt][iorder] += SV_dN_dypTdpTdphi[ipt][iphi]*sin(iorder*SP_pphi[iphi])*SP_pphi_weight[iphi];
				}
			}
			for(int ipt=0; ipt<n_SP_pT; ipt++)
			{
				cosine_iorder[ipt][iorder] /= SV_dN_dypTdpT[ipt];
				sine_iorder[ipt][iorder] /= SV_dN_dypTdpT[ipt];
			}
			cosine = cosine/norm;
			sine = sine/norm;
			if( sqrt(sine*sine + cosine*cosine) < 1e-8)
				plane_angle[iorder] = 0.0e0;
			else
				plane_angle[iorder] = atan2(sine, cosine)/double(iorder);
		}
		
		for(int ipt=0; ipt<n_SP_pT; ipt++)
			SV_dN_dypTdpT[ipt] /= (2.*M_PI);
		//cout << "Currently getting <p_T> stuff..." << endl;
		
		mean_pT = 0.;
		for(int iphi=0; iphi<n_SP_pphi; iphi++)
			mean_pT += SV_pTdN_dydphi[iphi]*SP_pphi_weight[iphi];
		mean_pT /= norm;
		plane_angle[0] = norm;
	}
	
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

void SourceVariances::Load_decay_channel_info(int dc_idx, double K_T_local, double K_phi_local)
{
	Mres = current_resonance_mass;
	Gamma = current_resonance_Gamma;
	double one_by_Gamma_Mres = hbarC/(Gamma*Mres);
	mass = current_daughter_mass;
	br = current_resonance_direct_br;	//doesn't depend on target daughter particle, just parent resonance and decay channel
	m2 = current_resonance_decay_masses[0];
	m3 = current_resonance_decay_masses[1];
	pT = K_T_local;
	current_K_phi = K_phi_local;
	cos_cKphi = cos(K_phi_local);
	sin_cKphi = sin(K_phi_local);
	n_body = current_reso_nbody;
	if (n_body == 2)
	{
		// some particles may decay to particles with more total mass than originally
		// --> broaden with resonance widths
		while ((mass + m2) > Mres)
		{
			Mres += 0.25 * current_resonance_Gamma;
		    mass -= 0.5 * current_daughter_Gamma;
		    m2 -= 0.5 * current_m2_Gamma;
		}

		mT = sqrt(mass*mass + pT*pT);

		//set up vectors of points to speed-up integrals...
		double s_loc = m2*m2;
		VEC_n2_spt = s_loc;
		double pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc) )/(2.0*Mres);
		VEC_n2_pstar = pstar_loc;
		check_for_NaNs("pstar_loc", pstar_loc, *global_out_stream_ptr);
		double g_s_loc = g(s_loc);	//for n_body == 2, doesn't actually use s_loc since result is just a factor * delta(...); just returns factor
		VEC_n2_s_factor = br/(4.*M_PI*VEC_n2_pstar);	//==g_s_loc
		double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
		VEC_n2_Estar = Estar_loc;
		double psBmT = pstar_loc / mT;
		VEC_n2_psBmT = psBmT;
		double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
		VEC_n2_DeltaY = DeltaY_loc;
		p_y = 0.0;
		VEC_n2_Yp = p_y + DeltaY_loc;
		VEC_n2_Ym = p_y - DeltaY_loc;
		check_for_NaNs("DeltaY_loc", DeltaY_loc, *global_out_stream_ptr);
/*if (VERBOSE > 0) *global_out_stream_ptr << "Working on resonance # = " << current_resonance_idx << ":" << endl
			<< setw(8) << setprecision(15)
			<< "  --> s_loc = " << s_loc << endl
			<< "  --> pstar_loc = " << pstar_loc << endl
			<< "  --> g_s_loc = " << g_s_loc << endl
			<< "  --> Estar_loc = " << Estar_loc << endl
			<< "  --> psBmT = " << psBmT << endl
			<< "  --> DeltaY_loc = " << DeltaY_loc << endl;*/
		for(int iv = 0; iv < n_v_pts; iv++)
		{
			//cerr << "In v loop# = " << iv << endl;
			double v_loc = v_pts[iv];
			double P_Y_loc = p_y + v_loc*DeltaY_loc;
			VEC_n2_P_Y[iv] = P_Y_loc;
			double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
			double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
			//if (x2 < 1e-15)
			//{
				//*global_out_stream_ptr << "Load_decay_channel_info(" << dc_idx << ", " << K_T_local << ", " << K_phi_local
				//	<< "): x2 = " << x2 << ", v = " << v_pts[iv] << endl;
				//x2 = 0.0;
			//}
			VEC_n2_v_factor[iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);
			double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
			VEC_n2_MTbar[iv] = MTbar_loc;
			double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
			VEC_n2_DeltaMT[iv] = DeltaMT_loc;
			VEC_n2_MTp[iv] = MTbar_loc + DeltaMT_loc;
			VEC_n2_MTm[iv] = MTbar_loc - DeltaMT_loc;
			check_for_NaNs("MTbar_loc", MTbar_loc, *global_out_stream_ptr);
			check_for_NaNs("DeltaMT_loc", DeltaMT_loc, *global_out_stream_ptr);
			//if (Estar_loc*Estar_loc - x2 < 1e-15)
			//{
				//*global_out_stream_ptr << "Load_decay_channel_info(" << dc_idx << ", " << K_T_local << ", " << K_phi_local
				//	<< "): Estar_loc*Estar_loc - x2 = " << Estar_loc*Estar_loc - x2 << ", v = " << v_pts[iv] << endl;
				//DeltaMT_loc = 0.0;
			//}

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
				check_for_NaNs("PT_loc", PT_loc, *global_out_stream_ptr);
				check_for_NaNs("PPhi_tilde_loc", PPhi_tilde_loc, *global_out_stream_ptr);
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
					check_for_NaNs("VEC_n2_alpha[iv][izeta][ii]", VEC_n2_alpha[iv][izeta][ii], *global_out_stream_ptr);
				}
			}
		}
	}
	else
	{
		mT = sqrt(mass*mass + pT*pT);
		double s_min_temp = (m2 + m3)*(m2 + m3);
		double s_max_temp = (Mres - mass)*(Mres - mass);
		gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min_temp, s_max_temp, NEW_s_pts, NEW_s_wts);
		Qfunc = get_Q();
/*if (VERBOSE > 0) *global_out_stream_ptr << "Working on resonance # = " << current_resonance_idx << ":" << endl
			<< setw(8) << setprecision(15)
			<< "  --> Qfunc = " << Qfunc << endl
			<< "  --> s_min_temp = " << s_min_temp << endl
			<< "  --> s_max_temp = " << s_max_temp << endl
			<< "  --> m2 = " << m2 << endl
			<< "  --> m3 = " << m3 << endl
			<< "  --> mass = " << mass << endl
			<< "  --> Mres = " << Mres << endl;*/
		//set up vectors of points to speed-up integrals...
		//cerr << "Entering loops now..." << endl;
		for (int is = 0; is < n_s_pts; is++)
		{
			//cerr << "In s loop# = " << is << endl;
			//double s_loc = s_pts[dc_idx-1][is];
			double s_loc = NEW_s_pts[is];
			double g_s_loc = g(s_loc);
			VEC_g_s[is] = g_s_loc;
			//VEC_s_factor[is] = s_wts[dc_idx-1][is]*g_s_loc;
			VEC_s_factor[is] = NEW_s_wts[is]*g_s_loc;
			double pstar_loc = sqrt(((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc))/(2.0*Mres);
			VEC_pstar[is] = pstar_loc;
			double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
			VEC_Estar[is] = Estar_loc;
			double psBmT = pstar_loc / mT;
			double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
/*if (VERBOSE > 0) *global_out_stream_ptr << "Working on resonance # = " << current_resonance_idx << ":" << endl
			<< setw(8) << setprecision(15)
			<< "  --> s_loc = " << s_loc << endl
			<< "  --> pstar_loc = " << pstar_loc << endl
			<< "  --> g_s_loc = " << g_s_loc << endl
			<< "  --> Estar_loc = " << Estar_loc << endl
			<< "  --> psBmT = " << psBmT << endl
			<< "  --> DeltaY_loc = " << DeltaY_loc << endl;*/
			VEC_DeltaY[is] = DeltaY_loc;
			p_y = 0.0;
			VEC_Yp[is] = p_y + DeltaY_loc;
			VEC_Ym[is] = p_y - DeltaY_loc;
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
						check_for_NaNs("VEC_alpha[iv][izeta][ii]", VEC_alpha[is][iv][izeta][ii], *global_out_stream_ptr);
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
		double mu = surf->particle_mu[particle_id];

		double vx = surf->vx;
		double vy = surf->vy;
		double Tdec = surf->Tdec;
		double one_by_Tdec = 1./Tdec;
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

					double expon = 0.0, f0 = 0.0;
					//now get distribution function, emission function, etc.
					if (TRUNCATE_COOPER_FRYE)
					{
						expon = (gammaT*(p0*1. - px*vx - py*vy) - mu)*one_by_Tdec;
						if (expon > 20.) continue;
						f0 = 1./(exp(expon)+sign);	//thermal equilibrium distributions
					}
					else
						f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions

					//p^mu d^3sigma_mu: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
					double pdsigma = p0*da0 + px*da1 + py*da2;

					//viscous corrections
					double Wfactor = p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
					double deltaf = 0.;
					if (use_delta_f)
						deltaf = (1. - sign*f0)*Wfactor*deltaf_prefactor;

					//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
					double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

					//ignore points where delta f is large or emission function goes negative from pdsigma
					if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
						S_p = 0.0e0;


					//double S_p = prefactor*pdsigma*f0*(1.+deltaf);
					double symmetry_factor = 1.0;
					if (ASSUME_ETA_SYMMETRIC) symmetry_factor = 2.0;
					double S_p_withweight = S_p*tau*eta_s_weight[ieta]*symmetry_factor; //symmetry_factor accounts for the assumed reflection symmetry along eta direction
					dN_dypTdpTdphi[ipt][iphi] += S_p_withweight;
				}	//end of ieta loop
			}		//end of iphi loop
		}			//end of ipt loop
	}				//end of isurf loop
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
	return 0.0;
}

void SourceVariances::Compute_source_variances(int iKT, int iKphi)
{
	double phi_K = K_phi[iKphi];
	double KT = K_T[iKT];
	if (USE_INTERP_ALT)
	{//rapidity-dependence not calculated here, but spectra are rap.-indep anyway
		S_func[iKT][iKphi] = Edndp3(KT, phi_K, 0.0, target_particle_id);
	}
	else
	{
		S_func[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][0],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
	}
	if (INCLUDE_SOURCE_VARIANCES)
	{
		xs_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][1],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xs2_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][2],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xo_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][3],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xo2_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][4],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xl_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][5],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xl2_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][6],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);		
		t_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][7],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		t2_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][8],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);		
		xo_xs_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][9],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xl_xs_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][10],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xs_t_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][11],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xo_xl_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][12],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xo_t_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][13],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
		xl_t_S[iKT][iKphi] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[target_particle_id][14],
							KT, phi_K, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING);
	}
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
		R2_side_C[iKT][Morder] = temp_sum_side_cos/(2.*M_PI);
		R2_side_S[iKT][Morder] = temp_sum_side_sin/(2.*M_PI);
		R2_out_C[iKT][Morder] = temp_sum_out_cos/(2.*M_PI);
		R2_out_S[iKT][Morder] = temp_sum_out_sin/(2.*M_PI);
		R2_outside_C[iKT][Morder] = temp_sum_outside_cos/(2.*M_PI);
		R2_outside_S[iKT][Morder] = temp_sum_outside_sin/(2.*M_PI);
		R2_long_C[iKT][Morder] = temp_sum_long_cos/(2.*M_PI);
		R2_long_S[iKT][Morder] = temp_sum_long_sin/(2.*M_PI);
		R2_sidelong_C[iKT][Morder] = temp_sum_sidelong_cos/(2.*M_PI);
		R2_sidelong_S[iKT][Morder] = temp_sum_sidelong_sin/(2.*M_PI);
		R2_outlong_C[iKT][Morder] = temp_sum_outlong_cos/(2.*M_PI);
		R2_outlong_S[iKT][Morder] = temp_sum_outlong_sin/(2.*M_PI);
	}
	return;
}

//End of file

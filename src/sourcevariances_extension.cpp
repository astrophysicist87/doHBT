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

void SourceVariances::Analyze_sourcefunction_V2(FO_surf* FOsurf_ptr)
{
	//this version uses adaptive Simpson's method routine to do resonance integrals
	//checking to see if it can be made any faster/more accurate than the Gaussian points integration of V1 above...
	time_t rawtime;
  	struct tm * timeinfo;

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

	//int decay_channel_loop_cutoff = 0;				//loop over direct pions only
	//int decay_channel_loop_cutoff = n_decay_channels;			//loop over direct pions and decay_channels
	int decay_channel_loop_cutoff = 1;					//other
	
	for (int idc = 0; idc <= decay_channel_loop_cutoff; idc++)
	{
		if (idc != 0 && idc != 1)
		{
			if (VERBOSE > 0) *global_out_stream_ptr << "idc = " << idc << ": skipping." << endl;
			continue;
		}
		if (INTERPOLATION_FORMAT == 0) return;
		for(int iKT = 0; iKT < n_localp_T; iKT++)
		{
			if (abs(K_T[iKT]) < 1.e-10) continue;
			if (iKT > 0) continue;
			for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
			{
				if (iKphi > 0) continue;
				time (&rawtime);
				timeinfo = localtime (&rawtime);
				if (VERBOSE > 0) *global_out_stream_ptr << "Starting KT = " << K_T[iKT] << " and Kphi = " << K_phi[iKphi]
					<< " for resonance #" << idc << " at " << asctime(timeinfo);
				if (idc > 0)
					Do_resonance_integrals_V2(iKT, iKphi, idc);
				Update_source_variances(iKT, iKphi, idc);

			}
		}
	}
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
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

//some declarations for the functions in this file
	const int jmax = 4;
	const double eps = 1.e-2;
	
	double one_by_Gamma_Mres, K_phi_local;
	double s_factor, v_factor, zeta_factor;
	double s_loc, g_s_loc, pstar_loc, Estar_loc, psBmT, DeltaY_loc, p_y = 0.0;
	double vsum = 0.0, v_loc, P_Y_loc, mT_ch_P_Y_p_y, x2, MTbar_loc, DeltaMT_loc;
	double zeta_loc, MT_loc, PT_loc;
	double temp_cos_PPhi_tilde_loc, temp_sin_PPhi_tilde_loc, PPhi_tilde_loc, PPhi_tilde_SHIFT, PPhi_tilde_SHIFT_FLIP;

	double alpha_t, alpha_o, alpha_s, alpha_o_m, alpha_s_m, alpha_l;
	double PKT, PKphi, PKY;

void SourceVariances::s_loop(double s_loc, double * ssum_vec)
{
	pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc) )/(2.0*Mres);
	g_s_loc = g(s_loc);	//for n_body == 2, doesn't actually use s_loc since result is just a factor * delta(...); just returns factor
	s_factor = g_s_loc;
	Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
	psBmT = pstar_loc / mT;
	DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));

	if (VERBOSE > 0) *global_out_stream_ptr << "Working on resonance # = " << current_decay_channel_idx << ":" << endl
		<< setw(8) << setprecision(15)
		<< "  --> s_loc = " << s_loc << endl
		<< "  --> pstar_loc = " << pstar_loc << endl
		<< "  --> g_s_loc = " << g_s_loc << endl
		<< "  --> Estar_loc = " << Estar_loc << endl
		<< "  --> psBmT = " << psBmT << endl
		<< "  --> DeltaY_loc = " << DeltaY_loc << endl;
	
	double * vsum_vec = new double [n_weighting_functions];
	set_to_zero(vsum_vec, n_weighting_functions);
	//adaptive_simpson_integration(&SourceVariances::v_loop, -1.0+eps, 1.0-eps, vsum_vec);
	//adaptive_simpson_integration(&SourceVariances::v_loop, -0.995556969790498, 0.995556969790498, vsum_vec);
	adaptive_simpson_integration(&SourceVariances::v_loop, -1.0+0.0001, 1.0-0.0001, vsum_vec);

	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		ssum_vec[wfi] += Mres * s_factor * vsum_vec[wfi];
	if (VERBOSE > 0) *global_out_stream_ptr << "   ssum_vec[0] = " << ssum_vec[0] << endl;
	
	return;
}

void SourceVariances::v_loop(double v_loc, double * vsum_vec)
{
	p_y = 0.0;
	P_Y_loc = p_y + v_loc*DeltaY_loc;
	mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
	x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
	MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
	DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
	v_factor = DeltaY_loc/sqrt(x2);

	if (VERBOSE > 0) *global_out_stream_ptr << "Working on resonance # = " << current_decay_channel_idx << ":" << endl
		<< setw(8) << setprecision(15)
		<< "  --> v_loc = " << v_loc << endl
		<< "  --> P_Y_loc = " << P_Y_loc << endl
		<< "  --> mT = " << mT << endl
		<< "  --> mT_ch_P_Y_p_y = " << mT_ch_P_Y_p_y << endl
		<< "  --> x2 = " << x2 << endl
		<< "  --> MTbar_loc = " << MTbar_loc << endl
		<< "  --> Estar_loc*Estar_loc - x2 = " << Estar_loc*Estar_loc - x2 << endl
		<< "  --> DeltaMT_loc = " << DeltaMT_loc << endl;
	
	double * zetasum_vec = new double [n_weighting_functions];
	set_to_zero(zetasum_vec, n_weighting_functions);
	//adaptive_simpson_integration(&SourceVariances::zeta_loop, 0.0, M_PI, zetasum_vec);
	adaptive_simpson_integration(&SourceVariances::zeta_loop, 0.0+0.0001, M_PI-0.0001, zetasum_vec);
	//adaptive_simpson_integration(&SourceVariances::zeta_loop, 0.0+0.00697909553292475, M_PI-0.00697909553292475, zetasum_vec);

	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		vsum_vec[wfi] += v_factor * zetasum_vec[wfi];
	if (VERBOSE > 0) *global_out_stream_ptr << "      vsum_vec[0] = " << vsum_vec[0] << endl;
	
	return;
}

void SourceVariances::zeta_loop(double zeta_loc, double * zetasum_vec)
{
	MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
	zeta_factor = MT_loc;
	PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
	temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
	//assume that PPhi_tilde is +ve in next step...
	temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
	PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), interp_pphi_min, interp_pphi_max);
	PPhi_tilde_SHIFT = place_in_range( K_phi_local + PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
	PPhi_tilde_SHIFT_FLIP = place_in_range( K_phi_local - PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);

	if (VERBOSE > 0) *global_out_stream_ptr << "Working on resonance # = " << current_decay_channel_idx << ":" << endl
		<< setw(8) << setprecision(15)
		<< "  --> zeta_loc = " << zeta_loc << endl
		<< "  --> MT_loc = " << MT_loc << endl
		<< "  --> PT_loc = " << PT_loc << endl
		<< "  --> PPhi_tilde_loc = " << PPhi_tilde_loc << endl;

	alpha_t = one_by_Gamma_Mres * MT_loc * cosh(P_Y_loc);
	alpha_o = one_by_Gamma_Mres * PT_loc * cos(PPhi_tilde_SHIFT);
	alpha_s = one_by_Gamma_Mres * PT_loc * sin(PPhi_tilde_SHIFT);
	alpha_o_m = one_by_Gamma_Mres * PT_loc * cos(PPhi_tilde_SHIFT_FLIP);
	alpha_s_m = one_by_Gamma_Mres * PT_loc * sin(PPhi_tilde_SHIFT_FLIP);
	alpha_l = one_by_Gamma_Mres * MT_loc * sinh(P_Y_loc);

	PKT = PT_loc;
	PKY = P_Y_loc;
	PKphi = PPhi_tilde_SHIFT;

	double * Csum_vec = new double [n_weighting_functions];
	double * rap_indep_y_of_r = new double [n_weighting_functions];
	double * y_of_r = new double [n_weighting_functions];
	set_to_zero(Csum_vec, n_weighting_functions);
	set_to_zero(rap_indep_y_of_r, n_weighting_functions);
	set_to_zero(y_of_r, n_weighting_functions);
	
	for (int tempidx = 1; tempidx <= 2; tempidx++)
	{
		if (tempidx != 1)
		{
			//Phi only changes sign, does NOT get shifted by pi!
			PKphi = PPhi_tilde_SHIFT_FLIP;		//also takes Pp --> Pm
			alpha_o = alpha_o_m;
			alpha_s = alpha_s_m;
		}
		//***this option is exact
		compute_rap_indep_spacetime_moments(current_FOsurf_ptr, current_decay_channel_idx, PKT, PKphi, rap_indep_y_of_r);
		//boost-invariant Cooper-Frye only yields boost-invariant spectra at Y=0
		//	==> shift source-variances appropriately to obtain rapidity dependence
		get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
		//now compute appropriate linear combinations (maybe shift these into preceding function eventually?)
		combine_sourcevariances(Csum_vec, y_of_r, alpha_t, alpha_o, alpha_s, alpha_l);
	}
	
	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		zetasum_vec[wfi] += zeta_factor * Csum_vec[wfi];
	if (VERBOSE > 0) *global_out_stream_ptr << "         zetasum_vec[0] = " << zetasum_vec[0]
					<< ", zeta_factor = " << zeta_factor << ", Csum_vec[0] = " << Csum_vec[0] << endl;
	
	return;
}


void SourceVariances::Do_resonance_integrals_V2(int iKT, int iKphi, int dc_idx)
{
	//if (VERBOSE > 2) *global_out_stream_ptr << "   Made it to do_all_integrals(): n_body = " << n_body << endl;
	time_t rawtime;
  	struct tm * timeinfo;
	//double * final_vec = new double [n_weighting_functions];
	//set_to_zero(final_vec, n_weighting_functions);
	double * ssum_vec = new double [n_weighting_functions];
	set_to_zero(ssum_vec, n_weighting_functions);

	current_decay_channel_idx = dc_idx;
	//cerr << "Entering Load_decay_channel_info for dc_idx = " << dc_idx << endl;
	if (dc_idx == 0)
	{
		muRES = particle_mu;
		signRES = particle_sign;
		gRES = particle_gspin;
		
		/*if (DEBUG)
		{
			cerr << "Working on resonance # = " << current_decay_channel_idx << ":" << endl
				<< "  --> muRES = " << muRES << endl
				<< "  --> signRES = " << signRES << endl
				<< "  --> gRES = " << gRES << endl;
		}*/
		
		return;
	}
	else
	{
		//cerr << "Entering loop in Load_decay_channel_info for dc_idx = " << dc_idx << endl;
		current_resonance_mass = decay_channels.resonance_mass[dc_idx-1];
		current_resonance_Gamma = decay_channels.resonance_Gamma[dc_idx-1];
		current_resonance_total_br = decay_channels.resonance_total_br[dc_idx-1];
		current_resonance_decay_masses[0] = decay_channels.resonance_decay_masses[dc_idx-1][0];
		current_resonance_decay_masses[1] = decay_channels.resonance_decay_masses[dc_idx-1][1];

		muRES = decay_channels.resonance_mu[dc_idx-1];
		signRES = decay_channels.resonance_sign[dc_idx-1];
		gRES = decay_channels.resonance_gspin[dc_idx-1];

		Mres = current_resonance_mass;
		mass = particle_mass;
		Gamma = current_resonance_Gamma;
		one_by_Gamma_Mres = hbarC/(Gamma*Mres);
		br = current_resonance_total_br;
		m2 = current_resonance_decay_masses[0];
		m3 = current_resonance_decay_masses[1];
		/*if (DEBUG)
		{
			cerr << "Working on resonance # = " << current_decay_channel_idx << ":" << endl
				<< "  --> muRES = " << muRES << endl
				<< "  --> signRES = " << signRES << endl
				<< "  --> gRES = " << gRES << endl
				<< "  --> Mres = " << Mres << endl
				<< "  --> mass = " << mass << endl
				<< "  --> Gamma = " << Gamma << endl
				<< "  --> br = " << br << endl
				<< "  --> m2 = " << m2 << endl
				<< "  --> m3 = " << m3 << endl << endl;
		}*/
		if (abs(m3) <= 1.e-6)
			n_body = 2;
		else
			n_body = 3;
		pT = K_T[iKT];
		mT = sqrt(mass*mass + pT*pT);
		K_phi_local = K_phi[iKphi];
		current_K_phi = K_phi_local;
		cos_cKphi = cos(K_phi_local);
		sin_cKphi = sin(K_phi_local);
		Qfunc = get_Q(dc_idx);
		
		if (n_body == 2)
		{
			//then g(s) is delta-function, skip s-integration entirely
    		s_loop(m2*m2, ssum_vec);
		}
		else if (n_body == 3)
		{
			cerr << "Skipping n_body == 3 stuff for simplicity..." << endl;
			return;
			double smin = (m2+m3)*(m2+m3);
			double smax = (Mres-mass)*(Mres-mass);
			adaptive_simpson_integration(&SourceVariances::s_loop, smin, smax, ssum_vec);
		}
	}

	for (int iweight = 0; iweight < n_weighting_functions; iweight++)
		integrated_spacetime_moments[dc_idx][iweight][iKT][iKphi] = ssum_vec[iweight];

	return;
}





void SourceVariances::combine_sourcevariances(double * output, double * input, double alpha_t, double alpha_o, double alpha_s, double alpha_l)
{
	//assumes all SVs being calculated
	//set_to_zero(output, 15);
	
	//[{1}_r]_{r-->\pi}
	output[0] += input[0];
	//[{xs}_r]_{r-->\pi}
	output[1] += input[1] + alpha_s*input[0];
	//[{xs2}_r]_{r-->\pi}
	output[2] += input[2] + 2.*alpha_s*input[1] + 2.*alpha_s*alpha_s*input[0];
	//[{xo}_r]_{r-->\pi}
	output[3] += input[3] + alpha_o*input[0];
	//[{xo2}_r]_{r-->\pi}
	output[4] += input[4] + 2.*alpha_o*input[3] + 2.*alpha_o*alpha_o*input[0];
	//[{xl}_r]_{r-->\pi}
	output[5] += input[5] + alpha_l*input[0];
	//[{xl2}_r]_{r-->\pi}
	output[6] += input[6] + 2.*alpha_l*input[5] + 2.*alpha_l*alpha_l*input[0];
	//[{t}_r]_{r-->\pi}
	output[7] += input[7] + alpha_t*input[0];
	//[{t2}_r]_{r-->\pi}
	output[8] += input[8] + 2.*alpha_t*input[7] + 2.*alpha_t*alpha_t*input[0];
	//[{xs_xo}_r]_{r-->\pi}
	output[9] += input[9] + alpha_s*input[3] + alpha_o*input[1] + 2.*alpha_s*alpha_o*input[0];
	//[{xs_xl}_r]_{r-->\pi}
	output[10] += input[10] + alpha_s*input[5] + alpha_l*input[1] + 2.*alpha_s*alpha_l*input[0];
	//[{xs_t}_r]_{r-->\pi}
	output[11] += input[11] + alpha_s*input[7] + alpha_t*input[1] + 2.*alpha_s*alpha_t*input[0];
	//[{xo_xl}_r]_{r-->\pi}
	output[12] += input[12] + alpha_o*input[5] + alpha_l*input[3] + 2.*alpha_o*alpha_l*input[0];
	//[{xo_t}_r]_{r-->\pi}
	output[13] += input[13] + alpha_o*input[7] + alpha_t*input[3] + 2.*alpha_o*alpha_t*input[0];
	//[{xl_t}_r]_{r-->\pi}
	output[14] += input[14] + alpha_l*input[7] + alpha_t*input[5] + 2.*alpha_l*alpha_t*input[0];
}

void SourceVariances::compute_rap_indep_spacetime_moments(FO_surf* FOsurf_ptr, int dc_idx, double KTres, double Kphires, double * rapidity_independent_y_of_r)
{
	//ALL calculations here done for Y == 0!!!
	// --> shift to non-zero rapidity afterwards
	double local_resonance_rapidity = 0.0;
	if (VERBOSE > 0) *global_out_stream_ptr << "            KTres = " << KTres << " and Kphires = " << Kphires << endl;
	
	double tpt, xopt, xspt, zpt;
	double sign, degen, localmass;
	if (dc_idx == 0)
	{
		cerr << "Shouldn't have gotten here!  Don't do resonance integrals for thermal pions!" << endl;
		exit(1);
	}
	else
	{
		sign = decay_channels.resonance_sign[dc_idx - 1];
		degen = decay_channels.resonance_gspin[dc_idx - 1];
		localmass = decay_channels.resonance_mass[dc_idx - 1];
	}
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	//these are constants along the FO surface,
	//so don't waste time updating them for each cell
	Tdec = (&FOsurf_ptr[0])->Tdec;
	double one_by_Tdec = 1./Tdec;
	Pdec = (&FOsurf_ptr[0])->Pdec;
	Edec = (&FOsurf_ptr[0])->Edec;
	double deltaf_prefactor = 0.;
	if (use_delta_f) deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	//declare variables used below here to avoid unnecessary redeclarations:
	double mu, tau, vx, vy, da0, da1, da2;
	double pi00, pi01, pi02, pi11, pi12, pi22, pi33;
	double temp_r, temp_phi, xpt, ypt, sin_temp_phi, cos_temp_phi, gammaT, expon;
	
	double sin_Kphires, cos_Kphires, sin_phi_m_Kphires, cos_phi_m_Kphires;
	double px, py, p0, pz, f0, deltaf, S_p, S_p_withweight;
	double MTres = sqrt(KTres*KTres + localmass*localmass);
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
		temp_phi = place_in_range(surf->phi, 0.0, 2.0*M_PI);
		xpt = surf->xpt;
		ypt = surf->ypt;
		sin_temp_phi = surf->sin_phi;
		cos_temp_phi = surf->cos_phi;
		gammaT = surf->gammaT;
		sin_Kphires = sin(Kphires);
		cos_Kphires = cos(Kphires);
		px = KTres*cos_Kphires;
		py = KTres*sin_Kphires;
		sin_phi_m_Kphires = sin_temp_phi * cos_Kphires - cos_temp_phi * sin_Kphires;
		cos_phi_m_Kphires = cos_temp_phi * cos_Kphires + sin_temp_phi * sin_Kphires;
		xopt = temp_r * cos_phi_m_Kphires;
		xspt = temp_r * sin_phi_m_Kphires;
		for(int ieta=0; ieta < eta_s_npts; ieta++)
		{
			//p0 = SPinterp_p0[ipt][ieta];
			//pz = SPinterp_pz[ipt][ieta];
			p0 = MTres * cosh(local_resonance_rapidity - eta_s[ieta]);
			pz = MTres * sinh(local_resonance_rapidity - eta_s[ieta]);
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
				//use analytic definition of S for code-checking
				S_p = prefactor * S_direct(temp_r, eta_s[ieta], tau, sqrt(KTres*KTres + localmass*localmass), KTres, cos_phi_m_Kphires);
			}
			else
			{
				//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
				S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);
				//ignore points where delta f is large or emission function goes negative from pdsigma
				if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
					S_p = 0.0;
			}
			if (flagneg == 1 && S_p < tol)
			{
				S_p = 0.0e0;
			}

			S_p_withweight = S_p*tau*eta_s_weight[ieta];
			rapidity_independent_y_of_r[0] += eta_even_factor*S_p_withweight;					//<1>
			rapidity_independent_y_of_r[1] += eta_even_factor*S_p_withweight*xspt;				//<x_s>
			rapidity_independent_y_of_r[2] += eta_even_factor*S_p_withweight*xspt*xspt;			//<x^2_s>
			rapidity_independent_y_of_r[3] += eta_even_factor*S_p_withweight*xopt;				//<x_o>
			rapidity_independent_y_of_r[4] += eta_even_factor*S_p_withweight*xopt*xopt;			//<x^2_o>
			rapidity_independent_y_of_r[5] += eta_odd_factor*S_p_withweight*zpt;				//<x_l>
			rapidity_independent_y_of_r[6] += eta_even_factor*S_p_withweight*zpt*zpt;			//<x^2_l>
			rapidity_independent_y_of_r[7] += eta_even_factor*S_p_withweight*tpt;				//<t>
			rapidity_independent_y_of_r[8] += eta_even_factor*S_p_withweight*tpt*tpt;			//<t^2>
			rapidity_independent_y_of_r[9] += eta_even_factor*S_p_withweight*xspt*xopt;			//<x_s x_o>
			rapidity_independent_y_of_r[10] += eta_odd_factor*S_p_withweight*xspt*zpt;			//<x_s x_l>
			rapidity_independent_y_of_r[11] += eta_even_factor*S_p_withweight*xspt*tpt;			//<x_s t>
			rapidity_independent_y_of_r[12] += eta_odd_factor*S_p_withweight*xopt*zpt;			//<x_o x_l>
			rapidity_independent_y_of_r[13] += eta_even_factor*S_p_withweight*xopt*tpt;			//<x_o t>
			rapidity_independent_y_of_r[14] += eta_odd_factor*S_p_withweight*zpt*tpt;			//<x_l t>
		}
	}
	//if (VERBOSE > 0) *global_out_stream_ptr << "            rapidity_independent_y_of_r[0] = " << rapidity_independent_y_of_r[0] << endl;

	return;
}

//***********************************************************************************************
//Some miscellaneous numerical routines that should be useful here
//***********************************************************************************************

//set all elements of a vector of given length to zero
/*void SourceVariances::set_to_zero(double * array, int arraylength)
{
	for (int arrayidx=0; arrayidx<arraylength; arrayidx++) array[arrayidx] = 0.0;
	
	return;
}*/

//implement adaptive simpson's integration routine for resonance integrals
//uses "divide-and-conquer" strategy
//more info available at: http://www-rohan.sdsu.edu/~jmahaffy/courses/s10/math541/lectures/pdf/week08/lecture.pdf
void SourceVariances::adaptive_simpson_integration(void (SourceVariances::*f) (double, double *), double a, double b, double * results)
{
	int devcount = 0;
	double relativedev = 0.0;
	double diff, d;
	//double osr = -1.e10, ostr = -1.e10;
	double del, x;
	double * tempa = new double [n_weighting_functions];
	double * tempb = new double [n_weighting_functions];
	double * tempx = new double [n_weighting_functions];
	double * str = new double [n_weighting_functions];
	double * sumr = new double [n_weighting_functions];
	double * osr = new double [n_weighting_functions];
	double * ostr = new double [n_weighting_functions];
	//int * refine_flags = new int [n_weighting_functions];
	diff = b - a;
	d = 0.5*diff;
	set_to_zero(results, n_weighting_functions);
	set_to_zero(str, n_weighting_functions);

	//str = d * (f(a) + f(b));
	(*this.*f)(a, tempa);
	(*this.*f)(b, tempb);
	for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		str[wfi] = d * (tempa[wfi] + tempb[wfi]);

	int it = 1;
	int finalorder = 0;
	for (int j = 2; j <= jmax; j++)
	{
		del = diff / (double)it;
		finalorder += it;
		x = a + 0.5 * del;
		//sumr = 0.0;
		set_to_zero(sumr, n_weighting_functions);
		for (int k = 1; k <= it; k++)
		{
			//sumr += f(x);
			(*this.*f)(x, tempx);
			for (int wfi = 0; wfi < n_weighting_functions; wfi++)
				sumr[wfi] += tempx[wfi];
			x += del;
			if (0) cout << "   --> inside adaptive_simpson_integration(): f(" << x << ") = " << tempx[0] << endl;
		}
		//str = 0.5 * (str + del * sumr);
		//sr = (4. * str - ostr) / 3.;
		for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		{
			str[wfi] = 0.5 * (str[wfi] + del * sumr[wfi]);
			results[wfi] = (4. * str[wfi] - ostr[wfi]) / 3.;
		}
		//for now, check for convergence based on spectra only (wfi == 0)
		if (EXIT_EARLY && abs(results[0] - osr[0]) < eps * abs(osr[0]))
		{
			if (0) cerr << j << endl;
			return;
		}
		/*devcount = 0;
		for (int wfi = 0; wfi < n_weighting_functions; wfi++)
			if (abs(results[wfi] - osr[wfi]) > eps * abs(osr[wfi]))
				devcount++;
		if (EXIT_EARLY && devcount == 0)
		{
			if (0) cerr << j << endl;
			return;
		}*/
		relativedev = abs(results[0]-osr[0]) / abs(osr[0]);
		for (int wfi = 0; wfi < n_weighting_functions; wfi++)
		{
			osr[wfi] = results[wfi];
			ostr[wfi] = str[wfi];
		}
		it *= 2;
	}

	if (0) cerr << "   --> exiting adaptive_simpson_integration(): " << relativedev << "   " << finalorder << endl;

	return;
}

//End of file

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

double SourceVariances::get_Q(int dc_idx)
{
	double smin = (m2+m3)*(m2+m3);
	double smax = (Mres-mass)*(Mres-mass);
	double sum = 0.;
	
	for (int is = 0; is < n_s_pts; is++)
	{
		double sp = s_pts[dc_idx-1][is];
		double f1 = (Mres+mass)*(Mres+mass) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		sum += s_wts[dc_idx-1][is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double SourceVariances::g(double s)
{
	double gs_pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s)*((Mres-mass)*(Mres-mass) - s) )/(2.0*Mres);
	double g_res = br/(4.*M_PI*gs_pstar_loc);
	if (n_body == 3)
	{
		double pre_f = (Mres*br)/(2.*M_PI*s);
		double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
		double den = Qfunc;
		g_res = pre_f * num / den;
	}
	//haven't treated (rare) case of n_body == 4 just yet...

	return g_res;
}

void SourceVariances::get_rapidity_dependence(double * rap_indep_vector, double * rap_dep_vector, double rap_val)
{
	//assumes currently calculating all 15 SVs
	//rap_indep_vector - SV vector evaluated at momentum rapidity y_p = 0
	//rap_dep_vector - SV vector evaluated at momentum rapidity y_p = rap_val
	//rap_val - value of rapidity to evaluate at
	
	double ch_rap_val = cosh(rap_val);
	double sh_rap_val = sinh(rap_val);
	
	//for SVs which don't change, just copy them over
	//[{1}_r]_{r-->\pi}
	rap_dep_vector[0] = rap_indep_vector[0];								//1
	//[{xs}_r]_{r-->\pi}
	rap_dep_vector[1] = rap_indep_vector[1];								//xs
	//[{xs2}_r]_{r-->\pi}
	rap_dep_vector[2] = rap_indep_vector[2];								//xs2
	//[{xo}_r]_{r-->\pi}
	rap_dep_vector[3] = rap_indep_vector[3];								//xo
	//[{xo2}_r]_{r-->\pi}
	rap_dep_vector[4] = rap_indep_vector[4];								//xo2
	//[{xl}_r]_{r-->\pi}
	rap_dep_vector[5] = ch_rap_val * rap_indep_vector[5]							//xl
				+ sh_rap_val * rap_indep_vector[7];						//t
	//[{xl2}_r]_{r-->\pi}
	rap_dep_vector[6] = ch_rap_val * ch_rap_val * rap_indep_vector[6]					//xl2
				+ 2. * ch_rap_val * sh_rap_val * rap_indep_vector[14]				//xlt
				+ sh_rap_val * sh_rap_val * rap_indep_vector[8];				//t2
	//[{t}_r]_{r-->\pi}
	rap_dep_vector[7] = ch_rap_val * rap_indep_vector[7]							//t
				+ sh_rap_val * rap_indep_vector[5];						//xl
	//[{t2}_r]_{r-->\pi}
	rap_dep_vector[8] = ch_rap_val * ch_rap_val * rap_indep_vector[8]					//t2
				+ 2. * ch_rap_val * sh_rap_val * rap_indep_vector[14]				//xlt
				+ sh_rap_val * sh_rap_val * rap_indep_vector[6];				//xl2
	//[{xs_xo}_r]_{r-->\pi}
	rap_dep_vector[9] = rap_indep_vector[9];								//xsxo
	//[{xs_xl}_r]_{r-->\pi}
	rap_dep_vector[10] = ch_rap_val * rap_indep_vector[10]							//xsxl
				+ sh_rap_val * rap_indep_vector[11];						//xst
	//[{xs_t}_r]_{r-->\pi}
	rap_dep_vector[11] = ch_rap_val * rap_indep_vector[11]							//xst
				+ sh_rap_val * rap_indep_vector[10];						//xsxl
	//[{xo_xl}_r]_{r-->\pi}
	rap_dep_vector[12] = ch_rap_val * rap_indep_vector[12]							//xoxl					
				+ sh_rap_val * rap_indep_vector[13];						//xot
	//[{xo_t}_r]_{r-->\pi}
	rap_dep_vector[13] = ch_rap_val * rap_indep_vector[13]							//xot
				+ sh_rap_val * rap_indep_vector[12];						//xoxl
	//[{xl_t}_r]_{r-->\pi}
	rap_dep_vector[14] = (ch_rap_val * ch_rap_val + sh_rap_val * sh_rap_val) * rap_indep_vector[14]		//xlt
				+ sh_rap_val * ch_rap_val * (rap_indep_vector[6] + rap_indep_vector[8]);	//t2
	
	return;
}


void SourceVariances::Do_resonance_integrals(int iKT, int iKphi, int dc_idx)
{
	//if (VERBOSE > 2) *global_out_stream_ptr << "   Made it to Do_resonance_integrals(): n_body = " << n_body << endl;
	time_t rawtime;
  	struct tm * timeinfo;
	int reso_idx = decay_channels.resonance_idx[dc_idx - 1] + 1;
	//double local_eta_s_wt = eta_s_weight[ieta];
	double * ssum_vec = new double [n_weighting_functions];
	double * vsum_vec = new double [n_weighting_functions];
	double * zetasum_vec = new double [n_weighting_functions];
	double * Csum_vec = new double [n_weighting_functions];
	double * rap_indep_y_of_r = new double [n_weighting_functions];
	double * y_of_r = new double [n_weighting_functions];
	set_to_zero(ssum_vec, n_weighting_functions);
	set_to_zero(vsum_vec, n_weighting_functions);
	set_to_zero(zetasum_vec, n_weighting_functions);
	set_to_zero(Csum_vec, n_weighting_functions);
	set_to_zero(rap_indep_y_of_r, n_weighting_functions);
	set_to_zero(y_of_r, n_weighting_functions);
	Qfunc = get_Q(dc_idx);

	if (n_body == 2)
	{
		//then g(s) is delta-function, skip s-integration entirely
		//double s_loc = m2*m2;
    		set_to_zero(vsum_vec, n_weighting_functions);
		for (int iv = 0; iv < n_v_pts; iv++)
		{
			double zetasum = 0.0;
			time (&rawtime);
			timeinfo = localtime (&rawtime);
			//cerr << "Starting v-loop #" << iv << " at " << asctime(timeinfo);
			set_to_zero(zetasum_vec, n_weighting_functions);
			for (int izeta = 0; izeta < n_zeta_pts; izeta++)
			{
				double alpha_t = VEC_n2_alpha[iv][izeta][0];
				double alpha_o = VEC_n2_alpha[iv][izeta][1];
				double alpha_s = VEC_n2_alpha[iv][izeta][2];
				double alpha_l = VEC_n2_alpha[iv][izeta][3];
				//cerr << "  Starting zeta-loop #" << izeta << endl;
				set_to_zero(Csum_vec, n_weighting_functions);
				double PKT = VEC_n2_PT[iv][izeta];
				double PKY = VEC_n2_P_Y[iv];
				double PKphi = VEC_n2_PPhi_tilde[iv][izeta];
				for (int tempidx = 1; tempidx <= 2; tempidx++)
				{
					if (tempidx != 1)
					{
						//Phi only changes sign, does NOT get shifted by pi!
						PKphi = VEC_n2_PPhi_tildeFLIP[iv][izeta];		//also takes Pp --> Pm
						alpha_s = VEC_n2_alpha_m[iv][izeta][2];
						alpha_o = VEC_n2_alpha_m[iv][izeta][1];
					}
					//instead of calculating each weight_function and averaging over FO surf a bazillion times,
					//just interpolate table of single particle spectra...
					//do interpolations
//if (0) cerr << iKT << "   " << iKphi << "   " << iv << "   " << izeta << "   " << tempidx << endl;
					for (int iweight = 0; iweight < n_weighting_functions; iweight++)
						rap_indep_y_of_r[iweight] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][iweight],
									PKT, PKphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING, true);
					get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
//if (dc_idx==42) cout << "-1   " << dc_idx << "   " << setw(8) << setprecision(15)
//			<< PKT << "   " << PKY << "   " << PKphi << "   " << y_of_r[0] << "   " << y_of_r[7]
//			<< "   " << y_of_r[8] << "   " << y_of_r[5] << "   " << y_of_r[6] << "   " << y_of_r[14] << endl;
//if (0) cout << "-2   " << dc_idx << "   " << setw(8) << setprecision(15)
//			<< PKT << "   " << PKY << "   " << PKphi << "   " << rap_indep_y_of_r[0] << "   " << rap_indep_y_of_r[7]
//			<< "   " << rap_indep_y_of_r[8] << "   " << rap_indep_y_of_r[5] << "   " << rap_indep_y_of_r[6] << "   " << rap_indep_y_of_r[14] << endl;
					//now compute appropriate linear combinations (maybe shift these into preceding function eventually?)
					//[{1}_r]_{r-->\pi}
					Csum_vec[0] += y_of_r[0];
					//[{xs}_r]_{r-->\pi}
					Csum_vec[1] += y_of_r[1] + alpha_s*y_of_r[0];
					//[{xs2}_r]_{r-->\pi}
					Csum_vec[2] += y_of_r[2] + 2.*alpha_s*y_of_r[1] + 2.*alpha_s*alpha_s*y_of_r[0];
					//[{xo}_r]_{r-->\pi}
					Csum_vec[3] += y_of_r[3] + alpha_o*y_of_r[0];
					//[{xo2}_r]_{r-->\pi}
					Csum_vec[4] += y_of_r[4] + 2.*alpha_o*y_of_r[3] + 2.*alpha_o*alpha_o*y_of_r[0];
					//[{xl}_r]_{r-->\pi}
					Csum_vec[5] += y_of_r[5] + alpha_l*y_of_r[0];
					//[{xl2}_r]_{r-->\pi}
					Csum_vec[6] += y_of_r[6] + 2.*alpha_l*y_of_r[5] + 2.*alpha_l*alpha_l*y_of_r[0];
//if (iKT == 0 && iKphi == 0)
//	if (0) cerr << v_pts[iv] << "   " << zeta_pts[izeta] << "   " << tempidx << "   " << K_T[iKT] << "   " << K_phi[iKphi] << "   "
//			<< PKT << "   " << PKphi << "   " << PKY << "   " << Csum_vec[0] << "   " << y_of_r[0] << endl;
					//[{t}_r]_{r-->\pi}
					Csum_vec[7] += y_of_r[7] + alpha_t*y_of_r[0];
					//[{t2}_r]_{r-->\pi}
					Csum_vec[8] += y_of_r[8] + 2.*alpha_t*y_of_r[7] + 2.*alpha_t*alpha_t*y_of_r[0];
					//[{xs_xo}_r]_{r-->\pi}
					Csum_vec[9] += y_of_r[9] + alpha_s*y_of_r[3] + alpha_o*y_of_r[1] + 2.*alpha_s*alpha_o*y_of_r[0];
					//[{xs_xl}_r]_{r-->\pi}
					Csum_vec[10] += y_of_r[10] + alpha_s*y_of_r[5] + alpha_l*y_of_r[1] + 2.*alpha_s*alpha_l*y_of_r[0];
					//[{xs_t}_r]_{r-->\pi}
					Csum_vec[11] += y_of_r[11] + alpha_s*y_of_r[7] + alpha_t*y_of_r[1] + 2.*alpha_s*alpha_t*y_of_r[0];
					//[{xo_xl}_r]_{r-->\pi}
					Csum_vec[12] += y_of_r[12] + alpha_o*y_of_r[5] + alpha_l*y_of_r[3] + 2.*alpha_o*alpha_l*y_of_r[0];
					//[{xo_t}_r]_{r-->\pi}
					Csum_vec[13] += y_of_r[13] + alpha_o*y_of_r[7] + alpha_t*y_of_r[3] + 2.*alpha_o*alpha_t*y_of_r[0];
					//[{xl_t}_r]_{r-->\pi}
					Csum_vec[14] += y_of_r[14] + alpha_l*y_of_r[7] + alpha_t*y_of_r[5] + 2.*alpha_l*alpha_t*y_of_r[0];
				}
				for (int iweight = 0; iweight < n_weighting_functions; iweight++)
					zetasum_vec[iweight] += VEC_n2_zeta_factor[iv][izeta]*Csum_vec[iweight];
				//cerr << "BIG DEBUG (dc_idx = " << dc_idx << "): " << Csum_vec[0] << endl;
			}
			for (int iweight = 0; iweight < n_weighting_functions; iweight++)
				vsum_vec[iweight] += VEC_n2_v_factor[iv]*zetasum_vec[iweight];
		}
		for (int iweight = 0; iweight < n_weighting_functions; iweight++)
			ssum_vec[iweight] += Mres*VEC_n2_s_factor*vsum_vec[iweight];
	}
	else if (n_body == 3)
	{
		for (int is = 0; is < n_s_pts; is++)
		{
			time (&rawtime);
			timeinfo = localtime (&rawtime);
			//cerr << "Starting s-loop #" << is << " at " << asctime(timeinfo);
			double vsum = 0.0;
    			set_to_zero(vsum_vec, n_weighting_functions);
			for (int iv = 0; iv < n_v_pts; iv++)
			{
				//double zetasum = 0.0;
				//cerr << "\tStarting v-loop #" << iv << endl;
				set_to_zero(zetasum_vec, n_weighting_functions);
				for (int izeta = 0; izeta < n_zeta_pts; izeta++)
				{
					double alpha_t = VEC_alpha[is][iv][izeta][0];
					double alpha_o = VEC_alpha[is][iv][izeta][1];
					double alpha_s = VEC_alpha[is][iv][izeta][2];
					double alpha_l = VEC_alpha[is][iv][izeta][3];
					//cerr << "\t\tStarting zeta-loop #" << izeta << endl;
					set_to_zero(Csum_vec, n_weighting_functions);
					double PKT = VEC_PT[is][iv][izeta];
					double PKY = VEC_P_Y[is][iv];
					double PKphi = VEC_PPhi_tilde[is][iv][izeta];
					for (int tempidx = 1; tempidx <= 2; tempidx++)
					{
						if (tempidx != 1)
						{
							PKphi = VEC_PPhi_tildeFLIP[is][iv][izeta];		//also takes Pp --> Pm
							//alpha_s *= -1.;
							alpha_s = VEC_alpha_m[is][iv][izeta][2];
							alpha_o = VEC_alpha_m[is][iv][izeta][1];
						}
						//instead of calculating each weight_function and averaging over FO surf a bazillion times,
						//just interpolate table of single particle spectra...
						for (int iweight = 0; iweight < n_weighting_functions; iweight++)
							rap_indep_y_of_r[iweight] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx][iweight],
										PKT, PKphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING, true);
						get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
//if (dc_idx==42) cout << "-1   " << dc_idx << "   " << setw(8) << setprecision(15)
//			<< PKT << "   " << PKY << "   " << PKphi << "   " << y_of_r[0] << "   " << y_of_r[7]
//			<< "   " << y_of_r[8] << "   " << y_of_r[5] << "   " << y_of_r[6] << "   " << y_of_r[14] << endl;
						//now compute appropriate linear combinations
						//[{1}_r]_{r-->\pi}
						Csum_vec[0] += y_of_r[0];
						//[{xs}_r]_{r-->\pi}
						Csum_vec[1] += y_of_r[1] + alpha_s*y_of_r[0];
						//[{xs2}_r]_{r-->\pi}
						Csum_vec[2] += y_of_r[2] + 2.*alpha_s*y_of_r[1] + 2.*alpha_s*alpha_s*y_of_r[0];
						//[{xo}_r]_{r-->\pi}
						Csum_vec[3] += y_of_r[3] + alpha_o*y_of_r[0];
						//[{xo2}_r]_{r-->\pi}
						Csum_vec[4] += y_of_r[4] + 2.*alpha_o*y_of_r[3] + 2.*alpha_o*alpha_o*y_of_r[0];
						//[{xl}_r]_{r-->\pi}
						Csum_vec[5] += y_of_r[5] + alpha_l*y_of_r[0];
						//[{xl2}_r]_{r-->\pi}
						Csum_vec[6] += y_of_r[6] + 2.*alpha_l*y_of_r[5] + 2.*alpha_l*alpha_l*y_of_r[0];
						//[{t}_r]_{r-->\pi}
						Csum_vec[7] += y_of_r[7] + alpha_t*y_of_r[0];
						//[{t2}_r]_{r-->\pi}
						Csum_vec[8] += y_of_r[8] + 2.*alpha_t*y_of_r[7] + 2.*alpha_t*alpha_t*y_of_r[0];
						//[{xs_xo}_r]_{r-->\pi}
						Csum_vec[9] += y_of_r[9] + alpha_s*y_of_r[3] + alpha_o*y_of_r[1] + 2.*alpha_s*alpha_o*y_of_r[0];
						//[{xs_xl}_r]_{r-->\pi}
						Csum_vec[10] += y_of_r[10] + alpha_s*y_of_r[5] + alpha_l*y_of_r[1] + 2.*alpha_s*alpha_l*y_of_r[0];
						//[{xs_t}_r]_{r-->\pi}
						Csum_vec[11] += y_of_r[11] + alpha_s*y_of_r[7] + alpha_t*y_of_r[1] + 2.*alpha_s*alpha_t*y_of_r[0];
						//[{xo_xl}_r]_{r-->\pi}
						Csum_vec[12] += y_of_r[12] + alpha_o*y_of_r[5] + alpha_l*y_of_r[3] + 2.*alpha_o*alpha_l*y_of_r[0];
						//[{xo_t}_r]_{r-->\pi}
						Csum_vec[13] += y_of_r[13] + alpha_o*y_of_r[7] + alpha_t*y_of_r[3] + 2.*alpha_o*alpha_t*y_of_r[0];
						//[{xl_t}_r]_{r-->\pi}
						Csum_vec[14] += y_of_r[14] + alpha_l*y_of_r[7] + alpha_t*y_of_r[5] + 2.*alpha_l*alpha_t*y_of_r[0];
					}
					for (int iweight = 0; iweight < n_weighting_functions; iweight++)
						zetasum_vec[iweight] += VEC_zeta_factor[is][iv][izeta]*Csum_vec[iweight];
				}
				for (int iweight = 0; iweight < n_weighting_functions; iweight++)
					vsum_vec[iweight] += VEC_v_factor[is][iv]*zetasum_vec[iweight];
			}
			for (int iweight = 0; iweight < n_weighting_functions; iweight++)
				ssum_vec[iweight] += Mres*VEC_s_factor[is]*vsum_vec[iweight];
		}
	}

	for (int iweight = 0; iweight < n_weighting_functions; iweight++)
		integrated_spacetime_moments[dc_idx][iweight][iKT][iKphi] = ssum_vec[iweight];

	//clean up
	delete [] ssum_vec;
	delete [] vsum_vec;
	delete [] zetasum_vec;
	delete [] Csum_vec;
	delete [] rap_indep_y_of_r;
	delete [] y_of_r;

	return;
}

void SourceVariances::Do_resonance_integrals_NEW(int parent_resonance_index, int daughter_particle_index, int decay_channel)
{
	time_t rawtime;
  	struct tm * timeinfo;
	//int reso_idx = decay_channels.resonance_idx[dc_idx - 1] + 1;
	//double local_eta_s_wt = eta_s_weight[ieta];
	double * ssum_vec = new double [n_weighting_functions];
	double * vsum_vec = new double [n_weighting_functions];
	double * zetasum_vec = new double [n_weighting_functions];
	double * Csum_vec = new double [n_weighting_functions];
	double * rap_indep_y_of_r = new double [n_weighting_functions];
	double * y_of_r = new double [n_weighting_functions];
	set_to_zero(ssum_vec, n_weighting_functions);
	set_to_zero(vsum_vec, n_weighting_functions);
	set_to_zero(zetasum_vec, n_weighting_functions);
	set_to_zero(Csum_vec, n_weighting_functions);
	set_to_zero(rap_indep_y_of_r, n_weighting_functions);
	set_to_zero(y_of_r, n_weighting_functions);
	Qfunc = get_Q(decay_channel);

	for (int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
	for (int ipphi = 0; ipphi < n_interp2_pphi_pts; ipphi++)
	{
		double local_pT = SPinterp2_pT[ipt];
		double local_pphi = SPinterp2_pphi[ipphi];
		set_to_zero(ssum_vec, n_weighting_functions);
		set_to_zero(vsum_vec, n_weighting_functions);
		set_to_zero(zetasum_vec, n_weighting_functions);
		set_to_zero(Csum_vec, n_weighting_functions);
		set_to_zero(rap_indep_y_of_r, n_weighting_functions);
		set_to_zero(y_of_r, n_weighting_functions);
		Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information
		if (n_body == 2)
		{
			//then g(s) is delta-function, skip s-integration entirely
			//double s_loc = m2*m2;
			set_to_zero(vsum_vec, n_weighting_functions);
				for (int iv = 0; iv < n_v_pts; iv++)
				{
					double zetasum = 0.0;
					time (&rawtime);
					timeinfo = localtime (&rawtime);
					set_to_zero(zetasum_vec, n_weighting_functions);
					for (int izeta = 0; izeta < n_zeta_pts; izeta++)
					{
						double alpha_t = VEC_n2_alpha[iv][izeta][0];
						double alpha_o = VEC_n2_alpha[iv][izeta][1];
						double alpha_s = VEC_n2_alpha[iv][izeta][2];
						double alpha_l = VEC_n2_alpha[iv][izeta][3];
						set_to_zero(Csum_vec, n_weighting_functions);
						double PKT = VEC_n2_PT[iv][izeta];
						double PKY = VEC_n2_P_Y[iv];
						double PKphi = VEC_n2_PPhi_tilde[iv][izeta];
						for (int tempidx = 1; tempidx <= 2; tempidx++)
						{
							if (tempidx != 1)
							{
								//Phi only changes sign, does NOT get shifted by pi!
								PKphi = VEC_n2_PPhi_tildeFLIP[iv][izeta];		//also takes Pp --> Pm
								alpha_s = VEC_n2_alpha_m[iv][izeta][2];
								alpha_o = VEC_n2_alpha_m[iv][izeta][1];
							}
							//instead of calculating each weight_function and averaging over FO surf a bazillion times,
							//just interpolate table of single particle spectra...
							//do interpolations
							for (int iweight = 0; iweight < n_weighting_functions; iweight++)
								rap_indep_y_of_r[iweight] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[parent_resonance_index][iweight],
											PKT, PKphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING, true);
							get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
							//now compute appropriate linear combinations (maybe shift these into preceding function eventually?)
							//[{1}_r]_{r-->\pi}
							Csum_vec[0] += y_of_r[0];
							//[{xs}_r]_{r-->\pi}
							Csum_vec[1] += y_of_r[1] + alpha_s*y_of_r[0];
							//[{xs2}_r]_{r-->\pi}
							Csum_vec[2] += y_of_r[2] + 2.*alpha_s*y_of_r[1] + 2.*alpha_s*alpha_s*y_of_r[0];
							//[{xo}_r]_{r-->\pi}
							Csum_vec[3] += y_of_r[3] + alpha_o*y_of_r[0];
							//[{xo2}_r]_{r-->\pi}
							Csum_vec[4] += y_of_r[4] + 2.*alpha_o*y_of_r[3] + 2.*alpha_o*alpha_o*y_of_r[0];
							//[{xl}_r]_{r-->\pi}
							Csum_vec[5] += y_of_r[5] + alpha_l*y_of_r[0];
							//[{xl2}_r]_{r-->\pi}
							Csum_vec[6] += y_of_r[6] + 2.*alpha_l*y_of_r[5] + 2.*alpha_l*alpha_l*y_of_r[0];
							//[{t}_r]_{r-->\pi}
							Csum_vec[7] += y_of_r[7] + alpha_t*y_of_r[0];
							//[{t2}_r]_{r-->\pi}
							Csum_vec[8] += y_of_r[8] + 2.*alpha_t*y_of_r[7] + 2.*alpha_t*alpha_t*y_of_r[0];
							//[{xs_xo}_r]_{r-->\pi}
							Csum_vec[9] += y_of_r[9] + alpha_s*y_of_r[3] + alpha_o*y_of_r[1] + 2.*alpha_s*alpha_o*y_of_r[0];
							//[{xs_xl}_r]_{r-->\pi}
							Csum_vec[10] += y_of_r[10] + alpha_s*y_of_r[5] + alpha_l*y_of_r[1] + 2.*alpha_s*alpha_l*y_of_r[0];
							//[{xs_t}_r]_{r-->\pi}
							Csum_vec[11] += y_of_r[11] + alpha_s*y_of_r[7] + alpha_t*y_of_r[1] + 2.*alpha_s*alpha_t*y_of_r[0];
							//[{xo_xl}_r]_{r-->\pi}
							Csum_vec[12] += y_of_r[12] + alpha_o*y_of_r[5] + alpha_l*y_of_r[3] + 2.*alpha_o*alpha_l*y_of_r[0];
							//[{xo_t}_r]_{r-->\pi}
							Csum_vec[13] += y_of_r[13] + alpha_o*y_of_r[7] + alpha_t*y_of_r[3] + 2.*alpha_o*alpha_t*y_of_r[0];
							//[{xl_t}_r]_{r-->\pi}
							Csum_vec[14] += y_of_r[14] + alpha_l*y_of_r[7] + alpha_t*y_of_r[5] + 2.*alpha_l*alpha_t*y_of_r[0];
						}// end of tempidx sum
						for (int iweight = 0; iweight < n_weighting_functions; iweight++)
							zetasum_vec[iweight] += VEC_n2_zeta_factor[iv][izeta]*Csum_vec[iweight];
					}// end of zeta sum
					for (int iweight = 0; iweight < n_weighting_functions; iweight++)
						vsum_vec[iweight] += VEC_n2_v_factor[iv]*zetasum_vec[iweight];
				}// end of v sum
				for (int iweight = 0; iweight < n_weighting_functions; iweight++)
					ssum_vec[iweight] += Mres*VEC_n2_s_factor*vsum_vec[iweight];
		}// end of nbody == 2
		else if (n_body == 3)
		{
			for (int is = 0; is < n_s_pts; is++)
			{
				double vsum = 0.0;
 		  			set_to_zero(vsum_vec, n_weighting_functions);
				for (int iv = 0; iv < n_v_pts; iv++)
				{
					set_to_zero(zetasum_vec, n_weighting_functions);
					for (int izeta = 0; izeta < n_zeta_pts; izeta++)
					{
						double alpha_t = VEC_alpha[is][iv][izeta][0];
						double alpha_o = VEC_alpha[is][iv][izeta][1];
						double alpha_s = VEC_alpha[is][iv][izeta][2];
						double alpha_l = VEC_alpha[is][iv][izeta][3];
						set_to_zero(Csum_vec, n_weighting_functions);
						double PKT = VEC_PT[is][iv][izeta];
						double PKY = VEC_P_Y[is][iv];
						double PKphi = VEC_PPhi_tilde[is][iv][izeta];
						for (int tempidx = 1; tempidx <= 2; tempidx++)
						{
							if (tempidx != 1)
							{
								PKphi = VEC_PPhi_tildeFLIP[is][iv][izeta];		//also takes Pp --> Pm
								//alpha_s *= -1.;
								alpha_s = VEC_alpha_m[is][iv][izeta][2];
								alpha_o = VEC_alpha_m[is][iv][izeta][1];
							}
							//instead of calculating each weight_function and averaging over FO surf a bazillion times,
							//just interpolate table of single particle spectra...
							for (int iweight = 0; iweight < n_weighting_functions; iweight++)
								rap_indep_y_of_r[iweight] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[parent_resonance_index][iweight],
											PKT, PKphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING, true);
							get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
							//now compute appropriate linear combinations
							//[{1}_r]_{r-->\pi}
							Csum_vec[0] += y_of_r[0];
							//[{xs}_r]_{r-->\pi}
							Csum_vec[1] += y_of_r[1] + alpha_s*y_of_r[0];
							//[{xs2}_r]_{r-->\pi}
							Csum_vec[2] += y_of_r[2] + 2.*alpha_s*y_of_r[1] + 2.*alpha_s*alpha_s*y_of_r[0];
							//[{xo}_r]_{r-->\pi}
							Csum_vec[3] += y_of_r[3] + alpha_o*y_of_r[0];
							//[{xo2}_r]_{r-->\pi}
							Csum_vec[4] += y_of_r[4] + 2.*alpha_o*y_of_r[3] + 2.*alpha_o*alpha_o*y_of_r[0];
							//[{xl}_r]_{r-->\pi}
							Csum_vec[5] += y_of_r[5] + alpha_l*y_of_r[0];
							//[{xl2}_r]_{r-->\pi}
							Csum_vec[6] += y_of_r[6] + 2.*alpha_l*y_of_r[5] + 2.*alpha_l*alpha_l*y_of_r[0];
							//[{t}_r]_{r-->\pi}
							Csum_vec[7] += y_of_r[7] + alpha_t*y_of_r[0];
							//[{t2}_r]_{r-->\pi}
							Csum_vec[8] += y_of_r[8] + 2.*alpha_t*y_of_r[7] + 2.*alpha_t*alpha_t*y_of_r[0];
							//[{xs_xo}_r]_{r-->\pi}
							Csum_vec[9] += y_of_r[9] + alpha_s*y_of_r[3] + alpha_o*y_of_r[1] + 2.*alpha_s*alpha_o*y_of_r[0];
							//[{xs_xl}_r]_{r-->\pi}
							Csum_vec[10] += y_of_r[10] + alpha_s*y_of_r[5] + alpha_l*y_of_r[1] + 2.*alpha_s*alpha_l*y_of_r[0];
							//[{xs_t}_r]_{r-->\pi}
							Csum_vec[11] += y_of_r[11] + alpha_s*y_of_r[7] + alpha_t*y_of_r[1] + 2.*alpha_s*alpha_t*y_of_r[0];
							//[{xo_xl}_r]_{r-->\pi}
							Csum_vec[12] += y_of_r[12] + alpha_o*y_of_r[5] + alpha_l*y_of_r[3] + 2.*alpha_o*alpha_l*y_of_r[0];
							//[{xo_t}_r]_{r-->\pi}
							Csum_vec[13] += y_of_r[13] + alpha_o*y_of_r[7] + alpha_t*y_of_r[3] + 2.*alpha_o*alpha_t*y_of_r[0];
							//[{xl_t}_r]_{r-->\pi}
							Csum_vec[14] += y_of_r[14] + alpha_l*y_of_r[7] + alpha_t*y_of_r[5] + 2.*alpha_l*alpha_t*y_of_r[0];
						}// end of tempidx sum
						for (int iweight = 0; iweight < n_weighting_functions; iweight++)
							zetasum_vec[iweight] += VEC_zeta_factor[is][iv][izeta]*Csum_vec[iweight];
					}// end of zeta sum
					for (int iweight = 0; iweight < n_weighting_functions; iweight++)
						vsum_vec[iweight] += VEC_v_factor[is][iv]*zetasum_vec[iweight];
				}// end of v sum
				for (int iweight = 0; iweight < n_weighting_functions; iweight++)
					ssum_vec[iweight] += Mres*VEC_s_factor[is]*vsum_vec[iweight];
			}// end of s sum
		}// end of nbody == 3
		else if (n_body == 4)
		{
			;
		}// end of nbody == 4

		for (int iweight = 0; iweight < n_weighting_functions; iweight++)
			dN_dypTdpTdphi_moments[daughter_particle_index][iweight][ipt][ipphi] += ssum_vec[iweight];
	}//end of pt, pphi loops

	//clean up
	delete [] ssum_vec;
	delete [] vsum_vec;
	delete [] zetasum_vec;
	delete [] Csum_vec;
	delete [] rap_indep_y_of_r;
	delete [] y_of_r;

	return;
}

void SourceVariances::set_to_zero(double * array, int arraylength)
{
	for (int arrayidx=0; arrayidx<arraylength; arrayidx++) array[arrayidx] = 0.0;
	
	return;
}

//End of file

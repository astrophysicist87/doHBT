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

double SourceVariances::get_Q()
{
	double smin = (m2+m3)*(m2+m3);
	double smax = (Mres-mass)*(Mres-mass);
	double sum = 0.;
	
	for (int is = 0; is < n_s_pts; ++is)
	{
		//double sp = s_pts[dc_idx-1][is];
		double sp = NEW_s_pts[is];
		double f1 = (Mres+mass)*(Mres+mass) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		//sum += s_wts[dc_idx-1][is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
		sum += NEW_s_wts[is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double SourceVariances::g(double s)
{
	double gs_pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s)*((Mres-mass)*(Mres-mass) - s) )/(2.0*Mres);
	double g_res = br/(4.*M_PI*gs_pstar_loc);
	if (n_body == 3 || n_body == 4)		//both set up to work the same way
	{
		double pre_f = (Mres*br)/(2.*M_PI*s);
		double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
		double den = Qfunc;
		g_res = pre_f * num / den;
	}

	return g_res;
}

void SourceVariances::get_rapidity_dependence(double * rap_indep_vector, double * rap_dep_vector, double rap_val)
{
	//assumes currently calculating all 15 SVs
	//rap_indep_vector - SV vector evaluated at momentum rapidity y_p = 0
	//rap_dep_vector - SV vector evaluated at momentum rapidity y_p = rap_val
	//rap_val - value of rapidity to evaluate at
	
	//for SVs which don't change, just copy them over
	//[{1}_r]_{r-->\pi}
	rap_dep_vector[0] = rap_indep_vector[0];									//1
	if (INCLUDE_SOURCE_VARIANCES)
	{
		double ch_rap_val = cosh(rap_val);
		double sh_rap_val = sinh(rap_val);
	
		//[{xs}_r]_{r-->\pi}
		rap_dep_vector[1] = rap_indep_vector[1];								//xs
		//[{xs2}_r]_{r-->\pi}
		rap_dep_vector[2] = rap_indep_vector[2];								//xs2
		//[{xo}_r]_{r-->\pi}
		rap_dep_vector[3] = rap_indep_vector[3];								//xo
		//[{xo2}_r]_{r-->\pi}
		rap_dep_vector[4] = rap_indep_vector[4];								//xo2
		//[{xl}_r]_{r-->\pi}
		rap_dep_vector[5] = ch_rap_val * rap_indep_vector[5]					//xl
					+ sh_rap_val * rap_indep_vector[7];							//t
		//[{xl2}_r]_{r-->\pi}
		rap_dep_vector[6] = ch_rap_val * ch_rap_val * rap_indep_vector[6]		//xl2
					+ 2. * ch_rap_val * sh_rap_val * rap_indep_vector[14]		//xlt
					+ sh_rap_val * sh_rap_val * rap_indep_vector[8];			//t2
		//[{t}_r]_{r-->\pi}
		rap_dep_vector[7] = ch_rap_val * rap_indep_vector[7]					//t
					+ sh_rap_val * rap_indep_vector[5];							//xl
		//[{t2}_r]_{r-->\pi}
		rap_dep_vector[8] = ch_rap_val * ch_rap_val * rap_indep_vector[8]		//t2
					+ 2. * ch_rap_val * sh_rap_val * rap_indep_vector[14]		//xlt
					+ sh_rap_val * sh_rap_val * rap_indep_vector[6];			//xl2
		//[{xs_xo}_r]_{r-->\pi}
		rap_dep_vector[9] = rap_indep_vector[9];								//xsxo	
		//[{xs_xl}_r]_{r-->\pi}
		rap_dep_vector[10] = ch_rap_val * rap_indep_vector[10]					//xsxl
					+ sh_rap_val * rap_indep_vector[11];						//xst
		//[{xs_t}_r]_{r-->\pi}
		rap_dep_vector[11] = ch_rap_val * rap_indep_vector[11]					//xst
					+ sh_rap_val * rap_indep_vector[10];						//xsxl
		//[{xo_xl}_r]_{r-->\pi}
		rap_dep_vector[12] = ch_rap_val * rap_indep_vector[12]					//xoxl					
					+ sh_rap_val * rap_indep_vector[13];						//xot
		//[{xo_t}_r]_{r-->\pi}
		rap_dep_vector[13] = ch_rap_val * rap_indep_vector[13]					//xot
					+ sh_rap_val * rap_indep_vector[12];						//xoxl
		//[{xl_t}_r]_{r-->\pi}
		rap_dep_vector[14] = (ch_rap_val * ch_rap_val + sh_rap_val * sh_rap_val) * rap_indep_vector[14]		//xlt
					+ sh_rap_val * ch_rap_val * (rap_indep_vector[6] + rap_indep_vector[8]);	//t2
	}
	
	return;
}

void SourceVariances::combine_sourcevariances(double * output, double * input, double * alpha_vec)
{
	//0, 1, 2, 3 --> t, o, s, l
	//[{1}_r]_{r-->\pi}
	output[0] += input[0];
	if (INCLUDE_SOURCE_VARIANCES)
	{
		//[{xs}_r]_{r-->\pi}
		output[1] += input[1] + alpha_vec[2]*input[0];
		//[{xs2}_r]_{r-->\pi}
		output[2] += input[2] + 2.*alpha_vec[2]*input[1] + 2.*alpha_vec[2]*alpha_vec[2]*input[0];
		//[{xo}_r]_{r-->\pi}
		output[3] += input[3] + alpha_vec[1]*input[0];
		//[{xo2}_r]_{r-->\pi}
		output[4] += input[4] + 2.*alpha_vec[1]*input[3] + 2.*alpha_vec[1]*alpha_vec[1]*input[0];
		//[{xl}_r]_{r-->\pi}
		output[5] += input[5] + alpha_vec[3]*input[0];
		//[{xl2}_r]_{r-->\pi}
		output[6] += input[6] + 2.*alpha_vec[3]*input[5] + 2.*alpha_vec[3]*alpha_vec[3]*input[0];
		//[{t}_r]_{r-->\pi}
		output[7] += input[7] + alpha_vec[0]*input[0];
		//[{t2}_r]_{r-->\pi}
		output[8] += input[8] + 2.*alpha_vec[0]*input[7] + 2.*alpha_vec[0]*alpha_vec[0]*input[0];
		//[{xs_xo}_r]_{r-->\pi}
		output[9] += input[9] + alpha_vec[2]*input[3] + alpha_vec[1]*input[1] + 2.*alpha_vec[2]*alpha_vec[1]*input[0];
		//[{xs_xl}_r]_{r-->\pi}
		output[10] += input[10] + alpha_vec[2]*input[5] + alpha_vec[3]*input[1] + 2.*alpha_vec[2]*alpha_vec[3]*input[0];
		//[{xs_t}_r]_{r-->\pi}
		output[11] += input[11] + alpha_vec[2]*input[7] + alpha_vec[0]*input[1] + 2.*alpha_vec[2]*alpha_vec[0]*input[0];
		//[{xo_xl}_r]_{r-->\pi}
		output[12] += input[12] + alpha_vec[1]*input[5] + alpha_vec[3]*input[3] + 2.*alpha_vec[1]*alpha_vec[3]*input[0];
		//[{xo_t}_r]_{r-->\pi}
		output[13] += input[13] + alpha_vec[1]*input[7] + alpha_vec[0]*input[3] + 2.*alpha_vec[1]*alpha_vec[0]*input[0];
		//[{xl_t}_r]_{r-->\pi}
		output[14] += input[14] + alpha_vec[3]*input[7] + alpha_vec[0]*input[5] + 2.*alpha_vec[3]*alpha_vec[0]*input[0];
	}
	
	return;
}

void SourceVariances::Do_resonance_integrals(int parent_resonance_particle_id, int daughter_particle_id, int decay_channel)
{
	time_t rawtime;
  	struct tm * timeinfo;
	double * ssum_vec = new double [n_weighting_functions];
	double * vsum_vec = new double [n_weighting_functions];
	double * zetasum_vec = new double [n_weighting_functions];
	double * Csum_vec = new double [n_weighting_functions];
	double * rap_indep_y_of_r = new double [n_weighting_functions];
	double * y_of_r = new double [n_weighting_functions];
	double * alphavec = new double [4];
	set_to_zero(ssum_vec, n_weighting_functions);
	set_to_zero(vsum_vec, n_weighting_functions);
	set_to_zero(zetasum_vec, n_weighting_functions);
	set_to_zero(Csum_vec, n_weighting_functions);
	set_to_zero(rap_indep_y_of_r, n_weighting_functions);
	set_to_zero(y_of_r, n_weighting_functions);

	res_sign_info = sign_of_dN_dypTdpTdphi_moments[parent_resonance_particle_id];
	res_log_info = ln_dN_dypTdpTdphi_moments[parent_resonance_particle_id];
	res_moments_info = dN_dypTdpTdphi_moments[parent_resonance_particle_id];

	//if (parent_resonance_particle_id == 75)
	//	current_level_of_output = 1;
	//else
	//	current_level_of_output = 0;
	//
	//if (current_level_of_output > 0) cout << all_particles[parent_resonance_particle_id].name << " to " << all_particles[daughter_particle_id].name << " (decay #): ";

	n_body = current_reso_nbody;

	if (n_body == 2)
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
if (ipt == 0 && ipphi == 0)
	current_level_of_output = 1;
else
	current_level_of_output = 0;
			set_to_zero(ssum_vec, n_weighting_functions);
			set_to_zero(vsum_vec, n_weighting_functions);
			set_to_zero(zetasum_vec, n_weighting_functions);
			set_to_zero(Csum_vec, n_weighting_functions);
			set_to_zero(rap_indep_y_of_r, n_weighting_functions);
			set_to_zero(y_of_r, n_weighting_functions);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

			//then g(s) is delta-function, skip s-integration entirely
			//double s_loc = m2*m2;
			for (int iv = 0; iv < n_v_pts; ++iv)
			{
				double zetasum = 0.0;
				time (&rawtime);
				timeinfo = localtime (&rawtime);
				set_to_zero(zetasum_vec, n_weighting_functions);
				for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					alphavec = VEC_n2_alpha[iv][izeta];
					set_to_zero(Csum_vec, n_weighting_functions);
					double PKT = VEC_n2_PT[iv][izeta];
					double PKY = VEC_n2_P_Y[iv];
					double PKphi = VEC_n2_PPhi_tilde[iv][izeta];
					for (int tempidx = 1; tempidx <= 2; ++tempidx)
					{
						if (tempidx != 1)
						{
							//Phi only changes sign, does NOT get shifted by pi!
							PKphi = VEC_n2_PPhi_tildeFLIP[iv][izeta];		//also takes Pp --> Pm
							alphavec = VEC_n2_alpha_m[iv][izeta];
						}
						Edndp3(PKT, PKphi, rap_indep_y_of_r);
						get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
						combine_sourcevariances(Csum_vec, y_of_r, alphavec);
					}																					// end of tempidx sum
					for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
						zetasum_vec[iweight] += VEC_n2_zeta_factor[iv][izeta]*Csum_vec[iweight];
				}																						// end of zeta sum
				for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
					vsum_vec[iweight] += VEC_n2_v_factor[iv]*zetasum_vec[iweight];
			}																							// end of v sum
			for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
				ssum_vec[iweight] += Mres*VEC_n2_s_factor*vsum_vec[iweight];

			//update all gridpoints for daughter moments
			for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
			{
				dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] += ssum_vec[iweight];
				double temp = dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi];
				ln_dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] = log(abs(temp));
				sign_of_dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] = sgn(temp);
			}
	
			if (isnan(dN_dypTdpTdphi_moments[daughter_particle_id][0][ipt][ipphi]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "dN_dypTdpTdphi_moments[daughter_particle_id][0][" << ipt << "][" << ipphi << "] = "
										<< setw(8) << setprecision(15) << dN_dypTdpTdphi_moments[daughter_particle_id][0][ipt][ipphi] << endl
										<< "  --> pt = " << local_pT << std::endl
										<< "  --> pphi = " << local_pphi << std::endl
										<< "daughter_particle_id = " << daughter_particle_id << endl
										<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
										<< "  --> Qfunc = " << Qfunc << endl
										<< "  --> n_body = " << n_body << endl
										<< "  --> gRES = " << gRES << endl
										<< "  --> Mres = " << Mres << endl
										<< "  --> mass = " << mass << endl
										<< "  --> Gamma = " << Gamma << endl
										<< "  --> br = " << br << endl
										<< "  --> m2 = " << m2 << endl
										<< "  --> m3 = " << m3 << endl << endl;
				exit(1);
			}
		}																								// end of pT, pphi loops
	}																									// end of nbody == 2
	else
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
			set_to_zero(ssum_vec, n_weighting_functions);
			set_to_zero(vsum_vec, n_weighting_functions);
			set_to_zero(zetasum_vec, n_weighting_functions);
			set_to_zero(Csum_vec, n_weighting_functions);
			set_to_zero(rap_indep_y_of_r, n_weighting_functions);
			set_to_zero(y_of_r, n_weighting_functions);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

			for (int is = 0; is < n_s_pts; ++is)
			{
				double vsum = 0.0;
 		  		set_to_zero(vsum_vec, n_weighting_functions);
				for (int iv = 0; iv < n_v_pts; ++iv)
				{
					set_to_zero(zetasum_vec, n_weighting_functions);
					for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
					{
						alphavec = VEC_alpha[is][iv][izeta];
						set_to_zero(Csum_vec, n_weighting_functions);
						double PKT = VEC_PT[is][iv][izeta];
						double PKY = VEC_P_Y[is][iv];
						double PKphi = VEC_PPhi_tilde[is][iv][izeta];
						for (int tempidx = 1; tempidx <= 2; ++tempidx)
						{
							if (tempidx != 1)
							{
								PKphi = VEC_PPhi_tildeFLIP[is][iv][izeta];		//also takes Pp --> Pm
								alphavec = VEC_alpha_m[is][iv][izeta];
							}
							Edndp3(PKT, PKphi, rap_indep_y_of_r);
							get_rapidity_dependence(rap_indep_y_of_r, y_of_r, PKY);
							//now compute appropriate linear combinations
							combine_sourcevariances(Csum_vec, y_of_r, alphavec);
						}																					// end of tempidx sum
						for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
							zetasum_vec[iweight] += VEC_zeta_factor[is][iv][izeta]*Csum_vec[iweight];
					}																						// end of zeta sum
					for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
						vsum_vec[iweight] += VEC_v_factor[is][iv]*zetasum_vec[iweight];
				}																							// end of v sum
				for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
					ssum_vec[iweight] += Mres*VEC_s_factor[is]*vsum_vec[iweight];
			}																								// end of s sum
			//update all gridpoints for daughter moments
			for (int iweight = 0; iweight < n_weighting_functions; ++iweight)
			{
				dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] += ssum_vec[iweight];
				double temp = dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi];
				ln_dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] = log(abs(temp));
				sign_of_dN_dypTdpTdphi_moments[daughter_particle_id][iweight][ipt][ipphi] = sgn(temp);
			}
	
			if (isnan(dN_dypTdpTdphi_moments[daughter_particle_id][0][ipt][ipphi]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "dN_dypTdpTdphi_moments[daughter_particle_id][0][" << ipt << "][" << ipphi << "] = "
										<< setw(8) << setprecision(15) << dN_dypTdpTdphi_moments[daughter_particle_id][0][ipt][ipphi] << endl
										<< "  --> pt = " << local_pT << std::endl
										<< "  --> pphi = " << local_pphi << std::endl
										<< "daughter_particle_id = " << daughter_particle_id << endl
										<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
										<< "  --> Qfunc = " << Qfunc << endl
										<< "  --> n_body = " << n_body << endl
										<< "  --> gRES = " << gRES << endl
										<< "  --> Mres = " << Mres << endl
										<< "  --> mass = " << mass << endl
										<< "  --> Gamma = " << Gamma << endl
										<< "  --> br = " << br << endl
										<< "  --> m2 = " << m2 << endl
										<< "  --> m3 = " << m3 << endl << endl;
				exit(1);
			}
		}																									// end of pT, pphi loops
	}																										// end of nbody == 3



	//clean up
	delete [] ssum_vec;
	delete [] vsum_vec;
	delete [] zetasum_vec;
	delete [] Csum_vec;
	delete [] rap_indep_y_of_r;
	delete [] y_of_r;

	return;
}

inline void SourceVariances::set_to_zero(double * array, size_t arraylength)
{
	for (size_t arrayidx=0; arrayidx<arraylength; ++arrayidx) array[arrayidx] = 0.0;
}

//End of file

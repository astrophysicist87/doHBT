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
	
	for (int is = 0; is < n_s_pts; is++)
	{
		double sp = s_pts[current_resonance_idx-1][is];
		double f1 = (Mres+mass)*(Mres+mass) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		sum += s_wts[current_resonance_idx-1][is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double SourceVariances::g(double s)
{
	double g_res = br/(4.*M_PI);
	if (n_body == 3)
	{
		double pre_f = (Mres*br)/(2.*M_PI*s);
		double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
		double den = Qfunc;
		double g_res = pre_f * num / den;
	}
	//haven't treated (rare) case of n_body == 4 just yet...

	return g_res;
}


void SourceVariances::Do_resonance_integrals(int iKT, int iKphi, int reso_idx)
{
	if (VERBOSE > 2) *global_out_stream_ptr << "   Made it to do_all_integrals(): n_body = " << n_body << endl;
	time_t rawtime;
  	struct tm * timeinfo;
	double ssum = 0.;
	//double local_eta_s_wt = eta_s_weight[ieta];
	double * ssum_vec = new double [n_weighting_functions];
	double * vsum_vec = new double [n_weighting_functions];
	double * zetasum_vec = new double [n_weighting_functions];
	double * Csum_vec = new double [n_weighting_functions];
	double * y_of_r = new double [n_weighting_functions];
	set_to_zero(ssum_vec, n_weighting_functions);
	set_to_zero(vsum_vec, n_weighting_functions);
	set_to_zero(zetasum_vec, n_weighting_functions);
	set_to_zero(Csum_vec, n_weighting_functions);
	set_to_zero(y_of_r, n_weighting_functions);
	Qfunc = get_Q();
	//double one_by_Gamma_M = 1./(Gamma*Mres);

	if (n_body == 2)
	{
		//then g(s) is delta-function, skip s-integration entirely
		//double s_loc = m2*m2;
		ssum = 0.0;
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
				double alpha_s = VEC_n2_alpha[iv][izeta][2];
				double alpha_o = VEC_n2_alpha[iv][izeta][1];
				double alpha_l = VEC_n2_alpha[iv][izeta][3];
				//cerr << "  Starting zeta-loop #" << izeta << endl;
				set_to_zero(Csum_vec, n_weighting_functions);
				double PKT = VEC_n2_PT[iv][izeta];
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
						y_of_r[iweight] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx-1][iweight],
									PKT, PKphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING, true);
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
					Csum_vec[4] += y_of_r[4] + 2.*alpha_o*y_of_r[1] + 2.*alpha_o*alpha_o*y_of_r[0];
					//[{xl}_r]_{r-->\pi}
					Csum_vec[5] += y_of_r[5] + alpha_l*y_of_r[0];
					//[{xl2}_r]_{r-->\pi}
					Csum_vec[6] += y_of_r[6] + 2.*alpha_l*y_of_r[1] + 2.*alpha_l*alpha_l*y_of_r[0];
					//[{t}_r]_{r-->\pi}
					Csum_vec[7] += y_of_r[7] + alpha_t*y_of_r[0];
					//[{t2}_r]_{r-->\pi}
					Csum_vec[8] += y_of_r[8] + 2.*alpha_t*y_of_r[1] + 2.*alpha_t*alpha_t*y_of_r[0];
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
				//cerr << "BIG DEBUG (reso_idx = " << reso_idx << "): " << Csum_vec[0] << endl;
			}
			for (int iweight = 0; iweight < n_weighting_functions; iweight++)
				vsum_vec[iweight] += VEC_n2_v_factor[iv]*zetasum_vec[iweight];
		}
		for (int iweight = 0; iweight < n_weighting_functions; iweight++)
			ssum_vec[iweight] += Mres*VEC_n2_s_factor*vsum_vec[iweight];
	}
	else if (n_body == 3)
	{
		ssum = 0.0;
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
					double alpha_s = VEC_alpha[is][iv][izeta][2];
					double alpha_o = VEC_alpha[is][iv][izeta][1];
					double alpha_l = VEC_alpha[is][iv][izeta][3];
					//cerr << "\t\tStarting zeta-loop #" << izeta << endl;
					set_to_zero(Csum_vec, n_weighting_functions);
					double PKT = VEC_PT[is][iv][izeta];
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
							y_of_r[iweight] = interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx-1][iweight],
										PKT, PKphi, n_interp2_pT_pts, n_interp2_pphi_pts, INTERPOLATION_KIND, UNIFORM_SPACING, true);
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
						Csum_vec[4] += y_of_r[4] + 2.*alpha_o*y_of_r[1] + 2.*alpha_o*alpha_o*y_of_r[0];
						//[{xl}_r]_{r-->\pi}
						Csum_vec[5] += y_of_r[5] + alpha_l*y_of_r[0];
						//[{xl2}_r]_{r-->\pi}
						Csum_vec[6] += y_of_r[6] + 2.*alpha_l*y_of_r[1] + 2.*alpha_l*alpha_l*y_of_r[0];
						//[{t}_r]_{r-->\pi}
						Csum_vec[7] += y_of_r[7] + alpha_t*y_of_r[0];
						//[{t2}_r]_{r-->\pi}
						Csum_vec[8] += y_of_r[8] + 2.*alpha_t*y_of_r[1] + 2.*alpha_t*alpha_t*y_of_r[0];
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

	//return ssum;
	for (int iweight = 0; iweight < n_weighting_functions; iweight++)
		integrated_spacetime_moments[reso_idx-1][iweight][iKT][iKphi] = ssum_vec[iweight];
	return;
}

void SourceVariances::set_to_zero(double * array, int arraylength)
{
	for (int arrayidx=0; arrayidx<arraylength; arrayidx++) array[arrayidx] = 0.0;
	
	return;
}

//End of file

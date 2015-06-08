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

void SourceVariances::set_pstar(double s)
{
	pstar = sqrt(((Mres+mass)*(Mres+mass) - s)*((Mres-mass)*(Mres-mass) - s)/(2.0*Mres));
	return;
}

void SourceVariances::set_Estar()
{
	Estar = sqrt(mass*mass + pstar*pstar);
	return;
}

void SourceVariances::set_DeltaY()
{
	double psBmT = pstar / mT;
	DeltaY = log(psBmT + sqrt(1.+psBmT*psBmT));
	return;
}

void SourceVariances::set_Ypm()
{
	Yp = p_y + DeltaY;
	Ym = p_y - DeltaY;
	return;
}

void SourceVariances::set_MTbar()
{
	double mT_ch_P_Y_p_y = mT*cosh(P_Y-p_y);
	double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
	double num = Estar*Mres*mT_ch_P_Y_p_y;
	MTbar = num/x2;
	return;
}

void SourceVariances::set_DeltaMT()
{
	double mT_ch_P_Y_p_y = mT*cosh(P_Y-p_y);
	double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
	double disc = Estar*Estar - x2;
	DeltaMT = Mres*pT*sqrt(disc)/x2;
	return;
}

void SourceVariances::set_MTpm()
{
	MTp = MTbar + DeltaMT;
	MTm = MTbar - DeltaMT;
	return;
}

void SourceVariances::set_MT(double zeta)
{
	MT = MTbar + zeta*DeltaMT;
	return;
}

void SourceVariances::set_P_Y(double v)
{
	P_Y = p_y + v*DeltaY;
	return;
}

void SourceVariances::set_PPhi_vars()
{	//assumes currently that pT != 0
	double cos_PPhi_tilde = (mT*MT*cosh(P_Y-p_y) - Estar*Mres)/(pT*PT+1.e-15);
	double PPhi_tilde = acos(cos_PPhi_tilde);
	PPhip = PPhi_tilde;
	PPhim = -PPhi_tilde;
	return;
}

void SourceVariances::set_Ppm()
{
	Pp[0] = MT*cosh(P_Y);
	Pp[1] = PT*cos(PPhip);
	Pp[2] = PT*sin(PPhip);
	Pp[3] = MT*sinh(P_Y);
	Pm[0] = MT*cosh(P_Y);
	Pm[1] = PT*cos(PPhim);
	Pm[2] = PT*sin(PPhim);
	Pm[3] = MT*sinh(P_Y);
	return;
}

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
	//assume n_body == 3...
	double pre_f = (Mres*br)/(2.*M_PI*s);
	double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
	//double den = get_Q();
	double den = Qfunc;
	double g_res = pre_f * num / den;

	return g_res;
}

double SourceVariances::s_integ(int wfi)
{
	double s_sum = 0.;
	Qfunc = get_Q();
	*global_out_stream_ptr << "\t\t\t + Made it to s_integ(): n_body = " << n_body << endl;
	
	if (n_body == 2)
	{
		//then g(s) is delta-function, skip s-integration entirely
		//double s_loc = m2*m2;
		set_pstar(m2*m2);
		set_Estar();
		set_DeltaY();
		set_Ypm();
		s_sum = br*v_integ(wfi)/(4.*M_PI*pstar);
	}
	else if (n_body == 3)
	{
		for (int is = 0; is < n_s_pts; is++)
		{
			double s_loc = s_pts[current_resonance_idx - 1][is];
			set_pstar(s_loc);
			set_Estar();
			set_DeltaY();
			set_Ypm();
	
			//g(s) only defined here for n_body == 3
			double s_factor = s_wts[current_resonance_idx - 1][is]*g(s_loc);
			s_sum += s_factor*v_integ(wfi);
		}
	}

	return s_sum;
}

double SourceVariances::v_integ(int wfi)
{
	double v_sum = 0.;
	
	for (int iv = 0; iv < n_v_pts; iv++)
	{
		double v_loc = v_pts[iv];
		set_P_Y(v_loc);
		set_MTbar();
		set_DeltaMT();
		set_MTpm();

		double v_factor = v_wts[iv]*DeltaY/(mT*mT*cosh(v_loc*DeltaY)*cosh(v_loc*DeltaY) - pT*pT);
		v_sum += v_factor*zeta_integ(wfi);
	}
	
	return v_sum;
}

double SourceVariances::zeta_integ(int wfi)
{
	double zeta_sum = 0.;
	
	for (int izeta = 0; izeta < n_zeta_pts; izeta++)
	{
		double zeta_loc = zeta_pts[izeta];
		set_MT(zeta_loc);
		set_PPhi_vars();
		set_Ppm();

		double zeta_factor = zeta_wts[izeta]*(MTbar + DeltaMT*cos(zeta_loc));
		//tau_integ no longer defined
		//zeta_sum += zeta_factor*tau_integ(wfi);
	}

	return zeta_sum;
}


/*double SourceVariances::do_all_integrals(int wfi, int ieta)
{
	*global_out_stream_ptr << "\t\t\t + Made it to do_all_integrals(): n_body = " << n_body << endl;
	time_t rawtime;
  	struct tm * timeinfo;
	double ssum = 0.;
	Qfunc = get_Q();
	double one_by_Gamma_M = 1./(Gamma*Mres);

	if (n_body == 2)
	{
		//then g(s) is delta-function, skip s-integration entirely
		//double s_loc = m2*m2;
		set_pstar(m2*m2);
		set_Estar();
		set_DeltaY();
		set_Ypm();
		ssum = br*v_integ(wfi)/(4.*M_PI*pstar);
	}
	else if (n_body == 3)
	{
		ssum = 0.0;
		for (int is = 0; is < n_s_pts; is++)
		{
			time (&rawtime);
			timeinfo = localtime (&rawtime);
			cerr << "Starting s-loop #" << is << " at " << asctime(timeinfo);
			double vsum = 0.0;
			for (int iv = 0; iv < n_v_pts; iv++)
			{
				double zetasum = 0.0;
				for (int izeta = 0; izeta < n_s_pts; izeta++)
				{
					double PK0, PK1, PK2, PK3, S_PK;
					double Csum = 0.0;
					PK0 = VEC_Pp[is][ieta][iv][izeta][0];
					PK1 = VEC_Pp[is][ieta][iv][izeta][1];
					PK2 = VEC_Pp[is][ieta][iv][izeta][2];
					PK3 = VEC_Pp[is][ieta][iv][izeta][3];
					//cerr << PK0 << "\t" << PK1 << "\t" << PK2 << "\t" << PK3 << endl;
					for (int tempidx = 1; tempidx <= 2; tempidx++)
					{
						if (tempidx != 1)
						PK2 *= -1.;		//takes Pp --> Pm
						for (int isurf = 0; isurf < FO_length; isurf++)
						{
							FO_surf* surf = &current_FOsurf_ptr[isurf];
							double S_PK = Emissionfunction(PK0, PK1, PK2, PK3, surf);
							//double S_PK = 1.;
							double surftau = surf->tau;
							double surfxpt = surf->xpt;
							double surfypt = surf->ypt;
							//double surftau = 1.;
							//double surfxpt = 1.;
							//double surfypt = 1.;
							
						//compute arguments of weight_function
							//use boost-invariance: p_y == local_eta_s
							zvec[0] = surftau*ch_p_y + PK0*one_by_Gamma_M;
							zvec[1] = surfxpt*cos_cKphi + surfypt*sin_cKphi + PK1*one_by_Gamma_M;
							zvec[2] = surfypt*cos_cKphi - surfxpt*sin_cKphi + PK2*one_by_Gamma_M;
							zvec[3] = surftau*sh_p_y + PK3*one_by_Gamma_M;
							//currently wrong!!!  just interested in approximate timing and optimization right now...
							Csum += S_PK*weight_function(zvec, wfi);
						}
					}
					zetasum += VEC_zeta_factor[is][iv][izeta]*Csum;
				}
				vsum += VEC_v_factor[is][iv]*zetasum;
			}
			ssum += VEC_s_factor[is]*vsum;
		}
	}

	return ssum;
}*/


void SourceVariances::do_all_integrals(int iKT, int iKphi, int reso_idx)
{
	*global_out_stream_ptr << "   Made it to do_all_integrals(): n_body = " << n_body << endl;
	time_t rawtime;
  	struct tm * timeinfo;
	double ssum = 0.;
	//double local_eta_s_wt = eta_s_weight[ieta];
	double * ssum_vec = new double [n_weighting_functions];
	double * vsum_vec = new double [n_weighting_functions];
	double * zetasum_vec = new double [n_weighting_functions];
	double * Csum_vec = new double [n_weighting_functions];
	set_to_zero(ssum_vec, n_weighting_functions);
	set_to_zero(vsum_vec, n_weighting_functions);
	set_to_zero(zetasum_vec, n_weighting_functions);
	set_to_zero(Csum_vec, n_weighting_functions);
	Qfunc = get_Q();
	double one_by_Gamma_M = 1./(Gamma*Mres);

	if (n_body == 2)
	{
		//then g(s) is delta-function, skip s-integration entirely
		//double s_loc = m2*m2;
		/*set_pstar(m2*m2);
		set_Estar();
		set_DeltaY();
		set_Ypm();
		ssum = br*v_integ(0)/(4.*M_PI*pstar);	//0 --> wfi*/
		ssum = 0.0;
    		set_to_zero(vsum_vec, n_weighting_functions);
		for (int iv = 0; iv < n_v_pts; iv++)
		{
			//double zetasum = 0.0;
			time (&rawtime);
			timeinfo = localtime (&rawtime);
			cerr << "Starting v-loop #" << iv << " at " << asctime(timeinfo);
			set_to_zero(zetasum_vec, n_weighting_functions);
			for (int izeta = 0; izeta < n_zeta_pts; izeta++)
			{
				double PK0, PK1, PK2, PK3, S_PK;
				cerr << "  Starting zeta-loop #" << izeta << endl;
				set_to_zero(Csum_vec, n_weighting_functions);
				PK0 = VEC_n2_Pp[iv][izeta][0];
				PK1 = VEC_n2_Pp[iv][izeta][1];
				PK2 = VEC_n2_Pp[iv][izeta][2];
				double PKT = VEC_n2_PT[iv][izeta];
				double PKphi = VEC_n2_PPhi_tilde[iv][izeta];
				PK3 = VEC_n2_Pp[iv][izeta][3];
				cerr << PK0 << "\t" << PK1 << "\t" << PK2 << "\t" << PK3 << endl;
				for (int tempidx = 1; tempidx <= 2; tempidx++)
				{
					if (tempidx != 1)
					{
						PK2 *= -1.;						//takes Pp --> Pm
						PKphi = VEC_n2_PPhi_tildeFLIP[iv][izeta];		//also takes Pp --> Pm
					}
					//instead of calculating each weight_function and averaging over FO surf a bazillion times,
					//just interpolate table of single particle spectra...
					for (int iweight = 0; iweight < n_weighting_functions; iweight++)
						Csum_vec[iweight] += interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx-1][iweight], PKT, PKphi, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING, true);
				}
				for (int iweight = 0; iweight < n_weighting_functions; iweight++)
					zetasum_vec[iweight] += VEC_n2_zeta_factor[iv][izeta]*Csum_vec[iweight];
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
			cerr << "Starting s-loop #" << is << " at " << asctime(timeinfo);
			double vsum = 0.0;
    			set_to_zero(vsum_vec, n_weighting_functions);
			for (int iv = 0; iv < n_v_pts; iv++)
			{
				//double zetasum = 0.0;
				cerr << "\tStarting v-loop #" << iv << endl;
				set_to_zero(zetasum_vec, n_weighting_functions);
				for (int izeta = 0; izeta < n_zeta_pts; izeta++)
				{
					double PK0, PK1, PK2, PK3, S_PK;
					cerr << "\t\tStarting zeta-loop #" << izeta << endl;
					set_to_zero(Csum_vec, n_weighting_functions);
					PK0 = VEC_Pp[is][iv][izeta][0];
					PK1 = VEC_Pp[is][iv][izeta][1];
					PK2 = VEC_Pp[is][iv][izeta][2];
					double PKT = VEC_PT[is][iv][izeta];
					double PKphi = VEC_PPhi_tilde[is][iv][izeta];
					PK3 = VEC_Pp[is][iv][izeta][3];
					cerr << PK0 << "\t" << PK1 << "\t" << PK2 << "\t" << PK3 << endl;
					for (int tempidx = 1; tempidx <= 2; tempidx++)
					{
						if (tempidx != 1)
						{
							PK2 *= -1.;						//takes Pp --> Pm
							PKphi = VEC_PPhi_tildeFLIP[is][iv][izeta];		//also takes Pp --> Pm
						}
						//instead of calculating each weight_function and averaging over FO surf a bazillion times,
						//just interpolate table of single particle spectra...
						for (int iweight = 0; iweight < n_weighting_functions; iweight++)
							Csum_vec[iweight] += interpolate2D(SPinterp2_pT, SPinterp2_pphi, dN_dypTdpTdphi_moments[reso_idx-1][iweight], PKT, PKphi, n_interp2_pT_pts, n_interp2_pphi_pts, 1, UNIFORM_SPACING, true);
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
	//SHOULD RETURN VECTOR OF SSUM_VEC SOURCEVARIANCES, BUT JUST CHECKING TIMING RIGHT NOW...
	for (int iweight = 0; iweight < n_weighting_functions; iweight++)
		integrated_spacetime_moments[reso_idx-1][iweight][iKT][iKphi] = ssum_vec[iweight];
	return;
}

void SourceVariances::set_to_zero(double * array, int arraylength)
{
	for (int arrayidx=0; arrayidx<arraylength; arrayidx++) array[arrayidx] = 0.0;
	
	return;
}

/*void SourceVariances::set_surfarrays()
{
	surf_tau_pts = new double [FO_length];
	surf_x_pts = new double [FO_length];
	surf_y_pts = new double [FO_length];

	for (int isurf = 0; isurf < FO_length; isurf++)
	{
		FO_surf* surf = &current_FOsurf_ptr[isurf];
		surf_tau_pts[isurf] = surf->tau;
		surf_x_pts[isurf] = surf->xpt;
		surf_y_pts[isurf] = surf->ypt;
	}
	
	return;
}*/

//End of file

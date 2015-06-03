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

#define USE_PLANE_PSI_ORDER 0		// specifies whether to do HBT relative to flow-plane angle,
					// and at what order: 0 - use plane_psi = 0.0, !0 - use flow-plane angle
					// at given order
//#define use_delta_f 0			// indicates whether to use delta_f corrections to distribution function
					// 0 - false
#define VERBOSE 2			// specifies level of output - 0 is lowest (no output)

using namespace std;

void doHBT::Update_avgSource_function(int i /*= -1*/, int j /*= -1*/)
{
//N.B. - avgs. only contains sums here, have not actually been averaged yet
   if (i < 0 || j < 0)
   {
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
   }
   else
   {
		avgS_func[i][j] += S_func[i][j];
		avgxs_S[i][j] += xs_S[i][j];
		avgxo_S[i][j] += xo_S[i][j];
		avgxl_S[i][j] += xl_S[i][j];
		avgt_S[i][j] += t_S[i][j];
		avgxs_t_S[i][j] += xs_t_S[i][j];
		avgxo_t_S[i][j] += xo_t_S[i][j];
		avgxl_t_S[i][j] += xl_t_S[i][j];
		avgxo_xs_S[i][j] += xo_xs_S[i][j];
		avgxl_xs_S[i][j] += xl_xs_S[i][j];
		avgxo_xl_S[i][j] += xo_xl_S[i][j];
		avgxs2_S[i][j] += xs2_S[i][j];
		avgxo2_S[i][j] += xo2_S[i][j];
		avgxl2_S[i][j] += xl2_S[i][j];
		avgt2_S[i][j] += t2_S[i][j];
   }
   return;
}


void doHBT::Update_azavgCavgSource_function(int i /*= -1*/)
{
//N.B. - avgs. only contains sums here, have not actually been averaged yet
   if (i < 0)
   {
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		azavg_Cavg_squared_S_func[iKT] += azavg_squared_S_func[iKT];
		azavg_CavgR2_side_num[iKT] += azavg_R2_side[iKT]*azavg_squared_S_func[iKT];
		azavg_CavgR2_out_num[iKT] += azavg_R2_out[iKT]*azavg_squared_S_func[iKT];
		azavg_CavgR2_outside_num[iKT] += azavg_R2_outside[iKT]*azavg_squared_S_func[iKT];
		azavg_CavgR2_long_num[iKT] += azavg_R2_long[iKT]*azavg_squared_S_func[iKT];
		azavg_CavgR2_sidelong_num[iKT] += azavg_R2_sidelong[iKT]*azavg_squared_S_func[iKT];
		azavg_CavgR2_outlong_num[iKT] += azavg_R2_outlong[iKT]*azavg_squared_S_func[iKT];
	}
   }
   else
   {
		azavg_Cavg_squared_S_func[i] += azavg_squared_S_func[i];
		azavg_CavgR2_side_num[i] += azavg_R2_side[i]*azavg_squared_S_func[i];
		azavg_CavgR2_out_num[i] += azavg_R2_out[i]*azavg_squared_S_func[i];
		azavg_CavgR2_outside_num[i] += azavg_R2_outside[i]*azavg_squared_S_func[i];
		azavg_CavgR2_long_num[i] += azavg_R2_long[i]*azavg_squared_S_func[i];
		azavg_CavgR2_sidelong_num[i] += azavg_R2_sidelong[i]*azavg_squared_S_func[i];
		azavg_CavgR2_outlong_num[i] += azavg_R2_outlong[i]*azavg_squared_S_func[i];
   }
   return;
}


void doHBT::Update_CavgSource_function(int i /*= -1*/, int j /*= -1*/)
{
//N.B. - avgs. only contains sums here, have not actually been averaged yet
   if (i < 0 || j < 0)
   {
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		CavgS_func_squared[iKT][iKphi] += S_func[iKT][iKphi]*S_func[iKT][iKphi];
		CavgR2_side_num[iKT][iKphi] += R2_side[iKT][iKphi]*S_func[iKT][iKphi]*S_func[iKT][iKphi];
		CavgR2_out_num[iKT][iKphi] += R2_out[iKT][iKphi]*S_func[iKT][iKphi]*S_func[iKT][iKphi];
		CavgR2_outside_num[iKT][iKphi] += R2_outside[iKT][iKphi]*S_func[iKT][iKphi]*S_func[iKT][iKphi];
		CavgR2_long_num[iKT][iKphi] += R2_long[iKT][iKphi]*S_func[iKT][iKphi]*S_func[iKT][iKphi];
		CavgR2_sidelong_num[iKT][iKphi] += R2_sidelong[iKT][iKphi]*S_func[iKT][iKphi]*S_func[iKT][iKphi];
		CavgR2_outlong_num[iKT][iKphi] += R2_outlong[iKT][iKphi]*S_func[iKT][iKphi]*S_func[iKT][iKphi];
	}
   }
   else
   {
		CavgS_func_squared[i][j] += S_func[i][j]*S_func[i][j];
		CavgR2_side_num[i][j] += R2_side[i][j]*S_func[i][j]*S_func[i][j];
		CavgR2_out_num[i][j] += R2_out[i][j]*S_func[i][j]*S_func[i][j];
		CavgR2_outside_num[i][j] += R2_outside[i][j]*S_func[i][j]*S_func[i][j];
		CavgR2_long_num[i][j] += R2_long[i][j]*S_func[i][j]*S_func[i][j];
		CavgR2_sidelong_num[i][j] += R2_sidelong[i][j]*S_func[i][j]*S_func[i][j];
		CavgR2_outlong_num[i][j] += R2_outlong[i][j]*S_func[i][j]*S_func[i][j];
   }
   return;
}

void doHBT::Calculate_avgSource_function(int i /*= -1*/, int j /*= -1*/)
{
//N.B. - avgs. only contains sums, doing averaging here
   if (i < 0 || j < 0)
   {
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
   }
   else
   {
		avgS_func[i][j] /= double(n_events);
		avgxs_S[i][j] /= double(n_events);
		avgxo_S[i][j] /= double(n_events);
		avgxl_S[i][j] /= double(n_events);
		avgt_S[i][j] /= double(n_events);
		avgxs_t_S[i][j] /= double(n_events);
		avgxo_t_S[i][j] /= double(n_events);
		avgxl_t_S[i][j] /= double(n_events);
		avgxo_xs_S[i][j] /= double(n_events);
		avgxl_xs_S[i][j] /= double(n_events);
		avgxo_xl_S[i][j] /= double(n_events);
		avgxs2_S[i][j] /= double(n_events);
		avgxo2_S[i][j] /= double(n_events);
		avgxl2_S[i][j] /= double(n_events);
		avgt2_S[i][j] /= double(n_events);
   }
   return;
}

void doHBT::Calculate_CavgSource_function(int i /*= -1*/, int j /*= -1*/)
{
//N.B. - avgs. only contains sums, doing averaging here
   if (i < 0 || j < 0)
   {
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		CavgS_func_squared[iKT][iKphi] /= double(n_events);
		CavgR2_side_num[iKT][iKphi] /= double(n_events);
		CavgR2_out_num[iKT][iKphi] /= double(n_events);
		CavgR2_outside_num[iKT][iKphi] /= double(n_events);
		CavgR2_long_num[iKT][iKphi] /= double(n_events);
		CavgR2_sidelong_num[iKT][iKphi] /= double(n_events);
		CavgR2_outlong_num[iKT][iKphi] /= double(n_events);
	}
   }
   else
   {
		CavgS_func_squared[i][j] /= double(n_events);
		CavgR2_side_num[i][j] /= double(n_events);
		CavgR2_out_num[i][j] /= double(n_events);
		CavgR2_outside_num[i][j] /= double(n_events);
		CavgR2_long_num[i][j] /= double(n_events);
		CavgR2_sidelong_num[i][j] /= double(n_events);
		CavgR2_outlong_num[i][j] /= double(n_events);
   }
   return;
}

void doHBT::Calculate_azavgCavgSource_function(int i /*= -1*/)
{
//N.B. - avgs. only contains sums, doing averaging here
   if (i < 0)
   {
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		azavg_Cavg_squared_S_func[iKT] /= double(n_events);
		azavg_CavgR2_side_num[iKT] /= double(n_events);
		azavg_CavgR2_out_num[iKT] /= double(n_events);
		azavg_CavgR2_outside_num[iKT] /= double(n_events);
		azavg_CavgR2_long_num[iKT] /= double(n_events);
		azavg_CavgR2_sidelong_num[iKT] /= double(n_events);
		azavg_CavgR2_outlong_num[iKT] /= double(n_events);
	}
   }
   else
   {
		azavg_Cavg_squared_S_func[i] /= double(n_events);
		azavg_CavgR2_side_num[i] /= double(n_events);
		azavg_CavgR2_out_num[i] /= double(n_events);
		azavg_CavgR2_outside_num[i] /= double(n_events);
		azavg_CavgR2_long_num[i] /= double(n_events);
		azavg_CavgR2_sidelong_num[i] /= double(n_events);
		azavg_CavgR2_outlong_num[i] /= double(n_events);
   }
   return;
}

void doHBT::Calculate_R2_side(int iKT, int iKphi)
{
   double norm = S_func[iKT][iKphi];
   double term1 = xs2_S[iKT][iKphi];
   double term2 = xs_S[iKT][iKphi];

   R2_side[iKT][iKphi] = term1/norm - term2*term2/(norm*norm);
   //cerr << "R^2_s(KT = " << K_T[iKT] << ", Kphi = " << K_phi[iKphi] << ") = " << R2_side[iKT][iKphi] << endl;
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

void doHBT::Calculate_CavgR2_side(int iKT, int iKphi)
{
   CavgR2_side[iKT][iKphi] = CavgR2_side_num[iKT][iKphi] / CavgS_func_squared[iKT][iKphi];
   return;
}

void doHBT::Calculate_CavgR2_out(int iKT, int iKphi)
{
   CavgR2_out[iKT][iKphi] = CavgR2_out_num[iKT][iKphi] / CavgS_func_squared[iKT][iKphi];
   return;
}

void doHBT::Calculate_CavgR2_outside(int iKT, int iKphi)
{
   CavgR2_outside[iKT][iKphi] = CavgR2_outside_num[iKT][iKphi] / CavgS_func_squared[iKT][iKphi];
   return;
}

void doHBT::Calculate_CavgR2_long(int iKT, int iKphi)
{
   CavgR2_long[iKT][iKphi] = CavgR2_long_num[iKT][iKphi] / CavgS_func_squared[iKT][iKphi];
   return;
}

void doHBT::Calculate_CavgR2_outlong(int iKT, int iKphi)
{
   CavgR2_outlong[iKT][iKphi] = CavgR2_outlong_num[iKT][iKphi] / CavgS_func_squared[iKT][iKphi];
   return;
}

void doHBT::Calculate_CavgR2_sidelong(int iKT, int iKphi)
{
   CavgR2_sidelong[iKT][iKphi] = CavgR2_sidelong_num[iKT][iKphi] / CavgS_func_squared[iKT][iKphi];
   return;
}

void doHBT::Calculate_azimuthally_averaged_squared_S_func(int iKT)
{
	azavg_squared_S_func[iKT] = 0.0;
	for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		squared_S_func[iKT][iKphi] = S_func[iKT][iKphi]*S_func[iKT][iKphi];
   		azavg_squared_S_func[iKT] += squared_S_func[iKT][iKphi]*K_phi_weight[iKphi];
	}

	azavg_squared_S_func[iKT] /= (2.*M_PI);

	return;
}

void doHBT::Calculate_azimuthally_averaged_R2_side(int iKT)
{
   azavg_R2_side[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_R2_side[iKT] += R2_side[iKT][iKphi]*squared_S_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_R2_side[iKT] /= (2.*M_PI*azavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_R2_out(int iKT)
{
   azavg_R2_out[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_R2_out[iKT] += R2_out[iKT][iKphi]*squared_S_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_R2_out[iKT] /= (2.*M_PI*azavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_R2_outside(int iKT)
{
   azavg_R2_outside[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_R2_outside[iKT] += R2_outside[iKT][iKphi]*squared_S_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_R2_outside[iKT] /= (2.*M_PI*azavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_R2_long(int iKT)
{
   azavg_R2_long[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_R2_long[iKT] += R2_long[iKT][iKphi]*squared_S_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_R2_long[iKT] /= (2.*M_PI*azavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_R2_outlong(int iKT)
{
   azavg_R2_outlong[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_R2_outlong[iKT] += R2_outlong[iKT][iKphi]*squared_S_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_R2_outlong[iKT] /= (2.*M_PI*azavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_R2_sidelong(int iKT)
{
   azavg_R2_sidelong[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_R2_sidelong[iKT] += R2_sidelong[iKT][iKphi]*squared_S_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_R2_sidelong[iKT] /= (2.*M_PI*azavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_squared_avgS_func(int iKT)
{
	if (iKT < 0)
	{//if iKT is negative, do it for all KT values
		for (int i = 0; i < n_localp_T; i++)
		{
			azavg_squared_avgS_func[i] = 0.0;
			
			for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
				azavg_squared_avgS_func[i] += avgS_func[i][iKphi]*avgS_func[i][iKphi]*K_phi_weight[iKphi];
			
			azavg_squared_avgS_func[i] /= (2.*M_PI);
		}
	}
	else
	{//otherwise, just do it for the specified KT value
		azavg_squared_avgS_func[iKT] = 0.0;
		
		for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
			azavg_squared_avgS_func[iKT] += avgS_func[iKT][iKphi]*avgS_func[iKT][iKphi]*K_phi_weight[iKphi];
		
		azavg_squared_avgS_func[iKT] /= (2.*M_PI);
	}

   return;
}

void doHBT::Calculate_azimuthally_averaged_avgR2_side(int iKT)
{
   azavg_avgR2_side[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_avgR2_side[iKT] += avgR2_side[iKT][iKphi]*avgS_func[iKT][iKphi]*avgS_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_avgR2_side[iKT] /= (2.*M_PI*azavg_squared_avgS_func[iKT]);

   return;
}

void doHBT::Calculate_azimuthally_averaged_avgR2_out(int iKT)
{
   azavg_avgR2_out[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_avgR2_out[iKT] += avgR2_out[iKT][iKphi]*avgS_func[iKT][iKphi]*avgS_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_avgR2_out[iKT] /= (2.*M_PI*azavg_squared_avgS_func[iKT]);

   return;
}

void doHBT::Calculate_azimuthally_averaged_avgR2_outside(int iKT)
{
   azavg_avgR2_outside[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_avgR2_outside[iKT] += avgR2_outside[iKT][iKphi]*avgS_func[iKT][iKphi]*avgS_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_avgR2_outside[iKT] /= (2.*M_PI*azavg_squared_avgS_func[iKT]);

   return;
}

void doHBT::Calculate_azimuthally_averaged_avgR2_long(int iKT)
{
   azavg_avgR2_long[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_avgR2_long[iKT] += avgR2_long[iKT][iKphi]*avgS_func[iKT][iKphi]*avgS_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_avgR2_long[iKT] /= (2.*M_PI*azavg_squared_avgS_func[iKT]);

   return;
}

void doHBT::Calculate_azimuthally_averaged_avgR2_outlong(int iKT)
{
   azavg_avgR2_outlong[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_avgR2_outlong[iKT] += avgR2_outlong[iKT][iKphi]*avgS_func[iKT][iKphi]*avgS_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_avgR2_outlong[iKT] /= (2.*M_PI*azavg_squared_avgS_func[iKT]);

   return;
}

void doHBT::Calculate_azimuthally_averaged_avgR2_sidelong(int iKT)
{
   azavg_avgR2_sidelong[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_avgR2_sidelong[iKT] += avgR2_sidelong[iKT][iKphi]*avgS_func[iKT][iKphi]*avgS_func[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_avgR2_sidelong[iKT] /= (2.*M_PI*azavg_squared_avgS_func[iKT]);

   return;
}

void doHBT::Calculate_azimuthally_averaged_Cavg_squared_S_func(int iKT)
{
	if (iKT < 0)
	{//if iKT is negative, do it for all KT values
		for (int i = 0; i < n_localp_T; i++)
		{
			azavg_Cavg_squared_S_func[i] = 0.0;
			
			for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
				azavg_Cavg_squared_S_func[i] += CavgS_func_squared[iKT][iKphi]*K_phi_weight[iKphi];
			
			azavg_Cavg_squared_S_func[i] /= (2.*M_PI);
		}
	}
	else
	{//otherwise, just do it for the specified KT value
		azavg_Cavg_squared_S_func[iKT] = 0.0;
		
		for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
			azavg_Cavg_squared_S_func[iKT] += CavgS_func_squared[iKT][iKphi]*K_phi_weight[iKphi];
		
		azavg_Cavg_squared_S_func[iKT] /= (2.*M_PI);
	}

   return;
}


void doHBT::Calculate_azimuthally_averaged_CavgR2_side(int iKT)
{
   //azavg_CavgR2_side[iKT] = azavg_CavgR2_side_num[iKT] / azavg_Cavg_squared_S_func[iKT];
   azavg_CavgR2_side[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_CavgR2_side[iKT] += CavgR2_side_num[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_CavgR2_side[iKT] /= (2.*M_PI*azavg_Cavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_CavgR2_out(int iKT)
{
   //azavg_CavgR2_out[iKT] = azavg_CavgR2_out_num[iKT] / azavg_Cavg_squared_S_func[iKT];
   azavg_CavgR2_out[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_CavgR2_out[iKT] += CavgR2_out_num[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_CavgR2_out[iKT] /= (2.*M_PI*azavg_Cavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_CavgR2_outside(int iKT)
{
   //azavg_CavgR2_outside[iKT] = azavg_CavgR2_outside_num[iKT] / azavg_Cavg_squared_S_func[iKT];
   azavg_CavgR2_outside[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_CavgR2_outside[iKT] += CavgR2_outside_num[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_CavgR2_outside[iKT] /= (2.*M_PI*azavg_Cavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_CavgR2_long(int iKT)
{
   //azavg_CavgR2_long[iKT] = azavg_CavgR2_long_num[iKT] / azavg_Cavg_squared_S_func[iKT];
   azavg_CavgR2_long[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_CavgR2_long[iKT] += CavgR2_long_num[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_CavgR2_long[iKT] /= (2.*M_PI*azavg_Cavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_CavgR2_outlong(int iKT)
{
   //azavg_CavgR2_outlong[iKT] = azavg_CavgR2_outlong_num[iKT] / azavg_Cavg_squared_S_func[iKT];
   azavg_CavgR2_outlong[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_CavgR2_outlong[iKT] += CavgR2_outlong_num[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_CavgR2_outlong[iKT] /= (2.*M_PI*azavg_Cavg_squared_S_func[iKT]);
   return;
}

void doHBT::Calculate_azimuthally_averaged_CavgR2_sidelong(int iKT)
{
   //azavg_CavgR2_sidelong[iKT] = azavg_CavgR2_sidelong_num[iKT] / azavg_Cavg_squared_S_func[iKT];
   azavg_CavgR2_sidelong[iKT] = 0.0;

   for (int iKphi = 0; iKphi < n_localp_phi; iKphi++)
   	azavg_CavgR2_sidelong[iKT] += CavgR2_sidelong_num[iKT][iKphi]*K_phi_weight[iKphi];

   azavg_CavgR2_sidelong[iKT] /= (2.*M_PI*azavg_Cavg_squared_S_func[iKT]);
   return;
}

//End of file

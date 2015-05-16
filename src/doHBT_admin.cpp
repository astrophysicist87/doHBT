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

doHBT::doHBT()
{
   //plumberg_test_variable = 3.14159;
   //default initial event and # of events for avg'ing routines
   initial_event = 1;	//defaults for now
   n_events = 1000;	//defaults for now
   checkpoint_index = 1;
   //for (int i=1; i<=1000; ++i) eventvector.push_back(i);
   //default: use delta_f in calculations
   use_delta_f = true;
   no_df_stem = "";
   //default: append to existing files rather than overwrite them
   append_output = true;

   Emissionfunction_ptr = new vector<Emissionfunction_data> (1);

   //single particle spectra for plane angle determination
   SP_pT = new double [n_SP_pT];
   SP_pT_weight = new double [n_SP_pT];
   SP_pphi = new double [n_SP_pphi];
   SP_pphi_weight = new double [n_SP_pphi];
   gauss_quadrature(n_SP_pphi, 1, 0.0, 0.0, 0.0, 2*M_PI, SP_pphi, SP_pphi_weight);
   if (SYMMETRIC_PT_PTS)	//use equally spaced points
   {
	double h = (SP_pT_max - SP_pT_min) / (double) (n_SP_pT - 1);
	for (int ipt=0; ipt<n_SP_pT; ipt++)
	{
		SP_pT[ipt] = SP_pT_min + double(ipt)*h;
		SP_pT_weight[ipt] = h;
	}
	SP_pT_weight[0] *= 0.5;
	SP_pT_weight[n_SP_pT-1] *= 0.5;
	//for (int ipt=0; ipt<n_SP_pT; ipt++) cout << SP_pT[ipt] << "   " << SP_pT_weight[ipt] << endl;
   }
   else				//use gaussian spaced points
   {
	gauss_quadrature(n_SP_pT, 1, 0.0, 0.0, SP_pT_min, SP_pT_max, SP_pT, SP_pT_weight);
   }
   SP_p_y = 0.0e0;
   dN_dypTdpTdphi = new double* [n_SP_pT];
   cosine_iorder = new double* [n_SP_pT];
   sine_iorder = new double* [n_SP_pT];
   anisotropic_flows_pTdiff = new double* [n_SP_pT];
   anisotropic_flows_pTdiff_psin = new double* [n_SP_pT];
   EdNd3p_cfs = new double* [n_SP_pT];
   EdNd3p_phases = new double* [n_SP_pT];
   for(int i=0; i<n_SP_pT; i++)
   {
      dN_dypTdpTdphi[i] = new double [n_SP_pphi];
      cosine_iorder[i] = new double [n_order];
      sine_iorder[i] = new double [n_order];
      anisotropic_flows_pTdiff[i] = new double [n_order];
      anisotropic_flows_pTdiff_psin[i] = new double [n_order];
      EdNd3p_cfs[i] = new double [n_order];
      EdNd3p_phases[i] = new double [n_order];
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
      anisotropic_flows_pTdiff[i][j] = 0.0e0;
      anisotropic_flows_pTdiff_psin[i][j] = 0.0e0;
   }
   for (int i=0; i<n_SP_pT; i++)
	dN_dypTdpT[i] = 0.0e0;
   plane_angle = new double [n_order];
   avgplane_angle = new double [n_order];
   Cavgplane_angle = new double [n_order];
   anisotropic_flows = new double [n_order];

   //relative momentum
   q_out = new double [qnpts];
   q_side = new double [qnpts];
   q_long = new double [qnpts];
   for(int i=0; i<qnpts; i++)
   {
      q_out[i] = init_q + (double)i * delta_q;
      q_side[i] = init_q + (double)i * delta_q;
      q_long[i] = init_q + (double)i * delta_q;
   }

   Correl_1D_out = new double [qnpts];
   Correl_1D_side = new double [qnpts];
   Correl_1D_long = new double [qnpts];
   Correl_1D_out_err = new double [qnpts];
   Correl_1D_side_err = new double [qnpts];
   Correl_1D_long_err = new double [qnpts];
   for(int i=0; i<qnpts; i++)
   {
      Correl_1D_out[i] = 0.0;
      Correl_1D_side[i] = 0.0;
      Correl_1D_long[i] = 0.0;
      Correl_1D_out_err[i] = 0.0;
      Correl_1D_side_err[i] = 0.0;
      Correl_1D_long_err[i] = 0.0;
   }

   Correl_3D = new double** [qnpts];
   Correl_3D_err = new double** [qnpts];
   for(int i=0; i<qnpts; i++)
   {
      Correl_3D[i] = new double* [qnpts];
      Correl_3D_err[i] = new double* [qnpts];
      for(int j=0; j<qnpts; j++)
         Correl_3D[i][j] = new double [qnpts];
      for(int j=0; j<qnpts; j++)
         Correl_3D_err[i][j] = new double [qnpts];
   }
   for(int i=0; i<qnpts; i++)
      for(int j=0; j<qnpts; j++)
         for(int k=0; k<qnpts; k++)
         {
            Correl_3D[i][j][k] = 0.0;
            Correl_3D_err[i][j][k] = 0.0;
         }

   //pair momentum
   K_T = new double [n_localp_T];
   double dK_T = (localp_T_max - localp_T_min)/(n_localp_T - 1 + 1e-100);
   for(int i=0; i<n_localp_T; i++) K_T[i] = localp_T_min + i*dK_T;
   //K_y = p_y;
   K_y = 0.;
   beta_l = sinh(K_y)/cosh(K_y);
   K_phi = new double [n_localp_phi];
   K_phi_weight = new double [n_localp_phi];
   gauss_quadrature(n_localp_phi, 1, 0.0, 0.0, localp_phi_min, localp_phi_max, K_phi, K_phi_weight);

   //spatial rapidity grid
   eta_s = new double [eta_s_npts];
   ch_eta_s = new double [eta_s_npts];
   sh_eta_s = new double [eta_s_npts];
   eta_s_weight = new double [eta_s_npts];
   gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, eta_s_i, eta_s_f, eta_s, eta_s_weight);
   for (int ieta = 0; ieta < eta_s_npts; ieta++)
   {
	ch_eta_s[ieta] = cosh(eta_s[ieta]);
	sh_eta_s[ieta] = sinh(eta_s[ieta]);
   }

   S_func = new double* [n_localp_T];
   squared_S_func = new double* [n_localp_T];
   xs_S = new double* [n_localp_T];
   xo_S = new double* [n_localp_T];
   xl_S = new double* [n_localp_T];
   t_S = new double* [n_localp_T];
   xs_t_S = new double* [n_localp_T];
   xo_t_S = new double* [n_localp_T];
   xl_t_S = new double* [n_localp_T];
   xo_xs_S = new double* [n_localp_T];
   xl_xs_S = new double* [n_localp_T];
   xo_xl_S = new double* [n_localp_T];
   xs2_S = new double* [n_localp_T];
   xo2_S = new double* [n_localp_T];
   xl2_S = new double* [n_localp_T];
   t2_S = new double* [n_localp_T];

   azavg_squared_S_func = new double [n_localp_T];
   azavg_S_func = new double [n_localp_T];

   xs_t_cos = new double* [n_localp_T];
   xo_t_cos = new double* [n_localp_T];
   xl_t_cos = new double* [n_localp_T];
   xo_xs_cos = new double* [n_localp_T];
   xl_xs_cos = new double* [n_localp_T];
   xo_xl_cos = new double* [n_localp_T];
   xs2_cos = new double* [n_localp_T];
   xo2_cos = new double* [n_localp_T];
   xl2_cos = new double* [n_localp_T];
   t2_cos = new double* [n_localp_T];
   xs_t_sin = new double* [n_localp_T];
   xo_t_sin = new double* [n_localp_T];
   xl_t_sin = new double* [n_localp_T];
   xo_xs_sin = new double* [n_localp_T];
   xl_xs_sin = new double* [n_localp_T];
   xo_xl_sin = new double* [n_localp_T];
   xs2_sin = new double* [n_localp_T];
   xo2_sin = new double* [n_localp_T];
   xl2_sin = new double* [n_localp_T];
   t2_sin = new double* [n_localp_T];

   avgS_func = new double* [n_localp_T];
   avgxs_S = new double* [n_localp_T];
   avgxo_S = new double* [n_localp_T];
   avgxl_S = new double* [n_localp_T];
   avgt_S = new double* [n_localp_T];
   avgxs_t_S = new double* [n_localp_T];
   avgxo_t_S = new double* [n_localp_T];
   avgxl_t_S = new double* [n_localp_T];
   avgxo_xs_S = new double* [n_localp_T];
   avgxl_xs_S = new double* [n_localp_T];
   avgxo_xl_S = new double* [n_localp_T];
   avgxs2_S = new double* [n_localp_T];
   avgxo2_S = new double* [n_localp_T];
   avgxl2_S = new double* [n_localp_T];
   avgt2_S = new double* [n_localp_T];

   azavg_squared_avgS_func = new double [n_localp_T];

   CavgS_func_squared = new double* [n_localp_T];
   CavgR2_side_num = new double* [n_localp_T];
   CavgR2_out_num = new double* [n_localp_T];
   CavgR2_long_num = new double* [n_localp_T];
   CavgR2_outside_num = new double* [n_localp_T];
   CavgR2_sidelong_num = new double* [n_localp_T];
   CavgR2_outlong_num = new double* [n_localp_T];

   azavg_Cavg_squared_S_func = new double [n_localp_T];
   azavg_CavgR2_side_num = new double [n_localp_T];
   azavg_CavgR2_out_num = new double [n_localp_T];
   azavg_CavgR2_long_num = new double [n_localp_T];
   azavg_CavgR2_outside_num = new double [n_localp_T];
   azavg_CavgR2_sidelong_num = new double [n_localp_T];
   azavg_CavgR2_outlong_num = new double [n_localp_T];

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

   azavg_R2_side = new double [n_localp_T];
   azavg_R2_out = new double [n_localp_T];
   azavg_R2_long = new double [n_localp_T];
   azavg_R2_outside = new double [n_localp_T];
   azavg_R2_sidelong = new double [n_localp_T];
   azavg_R2_outlong = new double [n_localp_T];

   R2_side_err = new double* [n_localp_T];
   R2_out_err = new double* [n_localp_T];
   R2_long_err = new double* [n_localp_T];
   R2_outside_err = new double* [n_localp_T];
   R2_sidelong_err = new double* [n_localp_T];
   R2_outlong_err = new double* [n_localp_T];

   lambda_Correl = new double* [n_localp_T];
   lambda_Correl_err = new double* [n_localp_T];

   avgR2_side = new double* [n_localp_T];
   avgR2_side_C = new double* [n_localp_T];
   avgR2_side_S = new double* [n_localp_T];
   avgR2_out = new double* [n_localp_T];
   avgR2_out_C = new double* [n_localp_T];
   avgR2_out_S = new double* [n_localp_T];
   avgR2_long = new double* [n_localp_T];
   avgR2_long_C = new double* [n_localp_T];
   avgR2_long_S = new double* [n_localp_T];
   avgR2_outside = new double* [n_localp_T];
   avgR2_outside_C = new double* [n_localp_T];
   avgR2_outside_S = new double* [n_localp_T];
   avgR2_sidelong = new double* [n_localp_T];
   avgR2_sidelong_C = new double* [n_localp_T];
   avgR2_sidelong_S = new double* [n_localp_T];
   avgR2_outlong = new double* [n_localp_T];
   avgR2_outlong_C = new double* [n_localp_T];
   avgR2_outlong_S = new double* [n_localp_T];

   azavg_avgR2_side = new double [n_localp_T];
   azavg_avgR2_out = new double [n_localp_T];
   azavg_avgR2_long = new double [n_localp_T];
   azavg_avgR2_outside = new double [n_localp_T];
   azavg_avgR2_sidelong = new double [n_localp_T];
   azavg_avgR2_outlong = new double [n_localp_T];

   CavgR2_side = new double* [n_localp_T];
   CavgR2_side_C = new double* [n_localp_T];
   CavgR2_side_S = new double* [n_localp_T];
   CavgR2_out = new double* [n_localp_T];
   CavgR2_out_C = new double* [n_localp_T];
   CavgR2_out_S = new double* [n_localp_T];
   CavgR2_long = new double* [n_localp_T];
   CavgR2_long_C = new double* [n_localp_T];
   CavgR2_long_S = new double* [n_localp_T];
   CavgR2_outside = new double* [n_localp_T];
   CavgR2_outside_C = new double* [n_localp_T];
   CavgR2_outside_S = new double* [n_localp_T];
   CavgR2_sidelong = new double* [n_localp_T];
   CavgR2_sidelong_C = new double* [n_localp_T];
   CavgR2_sidelong_S = new double* [n_localp_T];
   CavgR2_outlong = new double* [n_localp_T];
   CavgR2_outlong_C = new double* [n_localp_T];
   CavgR2_outlong_S = new double* [n_localp_T];

   azavg_CavgR2_side = new double [n_localp_T];
   azavg_CavgR2_out = new double [n_localp_T];
   azavg_CavgR2_long = new double [n_localp_T];
   azavg_CavgR2_outside = new double [n_localp_T];
   azavg_CavgR2_sidelong = new double [n_localp_T];
   azavg_CavgR2_outlong = new double [n_localp_T];

   for(int i=0; i<n_localp_T; i++)
   {
      S_func[i] = new double [n_localp_phi];
      squared_S_func[i] = new double [n_localp_phi];
      xs_S[i] = new double [n_localp_phi];
      xo_S[i] = new double [n_localp_phi];
      xl_S[i] = new double [n_localp_phi];
      t_S[i] = new double [n_localp_phi];
      xs_t_S[i] = new double [n_localp_phi];
      xo_t_S[i] = new double [n_localp_phi];
      xl_t_S[i] = new double [n_localp_phi];
      xo_xs_S[i] = new double [n_localp_phi];
      xl_xs_S[i] = new double [n_localp_phi];
      xo_xl_S[i] = new double [n_localp_phi];
      xs2_S[i] = new double [n_localp_phi];
      xo2_S[i] = new double [n_localp_phi];
      xl2_S[i] = new double [n_localp_phi];
      t2_S[i] = new double [n_localp_phi];

      azavg_squared_S_func[i] = 0.0;
      azavg_S_func[i] = 0.0;

      xs_t_cos[i] = new double [n_order];
      xo_t_cos[i] = new double [n_order];
      xl_t_cos[i] = new double [n_order];
      xo_xs_cos[i] = new double [n_order];
      xl_xs_cos[i] = new double [n_order];
      xo_xl_cos[i] = new double [n_order];
      xs2_cos[i] = new double [n_order];
      xo2_cos[i] = new double [n_order];
      xl2_cos[i] = new double [n_order];
      t2_cos[i] = new double [n_order];
      xs_t_sin[i] = new double [n_order];
      xo_t_sin[i] = new double [n_order];
      xl_t_sin[i] = new double [n_order];
      xo_xs_sin[i] = new double [n_order];
      xl_xs_sin[i] = new double [n_order];
      xo_xl_sin[i] = new double [n_order];
      xs2_sin[i] = new double [n_order];
      xo2_sin[i] = new double [n_order];
      xl2_sin[i] = new double [n_order];
      t2_sin[i] = new double [n_order];

      avgS_func[i] = new double [n_localp_phi];
      avgxs_S[i] = new double [n_localp_phi];
      avgxo_S[i] = new double [n_localp_phi];
      avgxl_S[i] = new double [n_localp_phi];
      avgt_S[i] = new double [n_localp_phi];
      avgxs_t_S[i] = new double [n_localp_phi];
      avgxo_t_S[i] = new double [n_localp_phi];
      avgxl_t_S[i] = new double [n_localp_phi];
      avgxo_xs_S[i] = new double [n_localp_phi];
      avgxl_xs_S[i] = new double [n_localp_phi];
      avgxo_xl_S[i] = new double [n_localp_phi];
      avgxs2_S[i] = new double [n_localp_phi];
      avgxo2_S[i] = new double [n_localp_phi];
      avgxl2_S[i] = new double [n_localp_phi];
      avgt2_S[i] = new double [n_localp_phi];

      CavgS_func_squared[i] = new double [n_localp_phi];
      CavgR2_side_num[i] = new double [n_localp_phi];
      CavgR2_out_num[i] = new double [n_localp_phi];
      CavgR2_long_num[i] = new double [n_localp_phi];
      CavgR2_outside_num[i] = new double [n_localp_phi];
      CavgR2_sidelong_num[i] = new double [n_localp_phi];
      CavgR2_outlong_num[i] = new double [n_localp_phi];

      azavg_Cavg_squared_S_func[i] = 0.0;
      azavg_CavgR2_side_num[i] = 0.0;
      azavg_CavgR2_out_num[i] = 0.0;
      azavg_CavgR2_long_num[i] = 0.0;
      azavg_CavgR2_outside_num[i] = 0.0;
      azavg_CavgR2_sidelong_num[i] = 0.0;
      azavg_CavgR2_outlong_num[i] = 0.0;

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

      azavg_R2_side[i] = 0.0;
      azavg_R2_out[i] = 0.0;
      azavg_R2_outside[i] = 0.0;
      azavg_R2_long[i] = 0.0;
      azavg_R2_sidelong[i] = 0.0;
      azavg_R2_outlong[i] = 0.0;

      R2_side_err[i] = new double [n_localp_phi];
      R2_out_err[i] = new double [n_localp_phi];
      R2_long_err[i] = new double [n_localp_phi];
      R2_outside_err[i] = new double [n_localp_phi];
      R2_sidelong_err[i] = new double [n_localp_phi];
      R2_outlong_err[i] = new double [n_localp_phi];

      lambda_Correl[i] = new double [n_localp_phi];
      lambda_Correl_err[i] = new double [n_localp_phi];

      avgR2_side[i] = new double [n_localp_phi];
      avgR2_side_C[i] = new double [n_order];
      avgR2_side_S[i] = new double [n_order];
      avgR2_out[i] = new double [n_localp_phi];
      avgR2_out_C[i] = new double [n_order];
      avgR2_out_S[i] = new double [n_order];
      avgR2_outside[i] = new double [n_localp_phi];
      avgR2_outside_C[i] = new double [n_order];
      avgR2_outside_S[i] = new double [n_order];
      avgR2_long[i] = new double [n_localp_phi];
      avgR2_long_C[i] = new double [n_order];
      avgR2_long_S[i] = new double [n_order];
      avgR2_sidelong[i] = new double [n_localp_phi];
      avgR2_sidelong_C[i] = new double [n_order];
      avgR2_sidelong_S[i] = new double [n_order];
      avgR2_outlong[i] = new double [n_localp_phi];
      avgR2_outlong_C[i] = new double [n_order];
      avgR2_outlong_S[i] = new double [n_order];

      azavg_avgR2_side[i] = 0.0;
      azavg_avgR2_out[i] = 0.0;
      azavg_avgR2_outside[i] = 0.0;
      azavg_avgR2_long[i] = 0.0;
      azavg_avgR2_sidelong[i] = 0.0;
      azavg_avgR2_outlong[i] = 0.0;

      CavgR2_side[i] = new double [n_localp_phi];
      CavgR2_side_C[i] = new double [n_order];
      CavgR2_side_S[i] = new double [n_order];
      CavgR2_out[i] = new double [n_localp_phi];
      CavgR2_out_C[i] = new double [n_order];
      CavgR2_out_S[i] = new double [n_order];
      CavgR2_outside[i] = new double [n_localp_phi];
      CavgR2_outside_C[i] = new double [n_order];
      CavgR2_outside_S[i] = new double [n_order];
      CavgR2_long[i] = new double [n_localp_phi];
      CavgR2_long_C[i] = new double [n_order];
      CavgR2_long_S[i] = new double [n_order];
      CavgR2_sidelong[i] = new double [n_localp_phi];
      CavgR2_sidelong_C[i] = new double [n_order];
      CavgR2_sidelong_S[i] = new double [n_order];
      CavgR2_outlong[i] = new double [n_localp_phi];
      CavgR2_outlong_C[i] = new double [n_order];
      CavgR2_outlong_S[i] = new double [n_order];

      azavg_CavgR2_side[i] = 0.0;
      azavg_CavgR2_out[i] = 0.0;
      azavg_CavgR2_outside[i] = 0.0;
      azavg_CavgR2_long[i] = 0.0;
      azavg_CavgR2_sidelong[i] = 0.0;
      azavg_CavgR2_outlong[i] = 0.0;
   }

//initialize all source variances and HBT radii/coeffs
for(int i=0; i<n_localp_T; i++)
{
   //reset azimuthally averaged radii first
      azavg_squared_S_func[i] = 0.0;
      azavg_S_func[i] = 0.0;

      azavg_Cavg_squared_S_func[i] = 0.0;
      azavg_CavgR2_side_num[i] = 0.0;
      azavg_CavgR2_out_num[i] = 0.0;
      azavg_CavgR2_long_num[i] = 0.0;
      azavg_CavgR2_outside_num[i] = 0.0;
      azavg_CavgR2_sidelong_num[i] = 0.0;
      azavg_CavgR2_outlong_num[i] = 0.0;

      azavg_R2_side[i] = 0.0;
      azavg_R2_out[i] = 0.0;
      azavg_R2_outside[i] = 0.0;
      azavg_R2_long[i] = 0.0;
      azavg_R2_sidelong[i] = 0.0;
      azavg_R2_outlong[i] = 0.0;

      azavg_avgR2_side[i] = 0.0;
      azavg_avgR2_out[i] = 0.0;
      azavg_avgR2_outside[i] = 0.0;
      azavg_avgR2_long[i] = 0.0;
      azavg_avgR2_sidelong[i] = 0.0;
      azavg_avgR2_outlong[i] = 0.0;

      azavg_CavgR2_side[i] = 0.0;
      azavg_CavgR2_out[i] = 0.0;
      azavg_CavgR2_outside[i] = 0.0;
      azavg_CavgR2_long[i] = 0.0;
      azavg_CavgR2_sidelong[i] = 0.0;
      azavg_CavgR2_outlong[i] = 0.0;

	for(int j=0; j<n_localp_phi; j++)
	{
		S_func[i][j] = 0.;
		squared_S_func[i][j] = 0.;
		xs_S[i][j] = 0.;
		xo_S[i][j] = 0.;
		xl_S[i][j] = 0.;
		t_S[i][j] = 0.;
		xs_t_S[i][j] = 0.;
		xo_t_S[i][j] = 0.;
		xl_t_S[i][j] = 0.;
		xo_xs_S[i][j] = 0.;
		xl_xs_S[i][j] = 0.;
		xo_xl_S[i][j] = 0.;
		xs2_S[i][j] = 0.;
		xo2_S[i][j] = 0.;
		xl2_S[i][j] = 0.;
		t2_S[i][j] = 0.;

		avgS_func[i][j] = 0.;
		avgxs_S[i][j] = 0.;
		avgxo_S[i][j] = 0.;
		avgxl_S[i][j] = 0.;
		avgt_S[i][j] = 0.;
		avgxs_t_S[i][j] = 0.;
		avgxo_t_S[i][j] = 0.;
		avgxl_t_S[i][j] = 0.;
		avgxo_xs_S[i][j] = 0.;
		avgxl_xs_S[i][j] = 0.;
		avgxo_xl_S[i][j] = 0.;
		avgxs2_S[i][j] = 0.;
		avgxo2_S[i][j] = 0.;
		avgxl2_S[i][j] = 0.;
		avgt2_S[i][j] = 0.;

		CavgS_func_squared[i][j] = 0.;
		CavgR2_side_num[i][j] = 0.;
		CavgR2_out_num[i][j] = 0.;
		CavgR2_outside_num[i][j] = 0.;
		CavgR2_long_num[i][j] = 0.;
		CavgR2_sidelong_num[i][j] = 0.;
		CavgR2_outlong_num[i][j] = 0.;

		R2_side[i][j] = 0.;
		R2_out[i][j] = 0.;
		R2_outside[i][j] = 0.;
		R2_long[i][j] = 0.;
		R2_sidelong[i][j] = 0.;
		R2_outlong[i][j] = 0.;

		avgR2_side[i][j] = 0.;
		avgR2_out[i][j] = 0.;
		avgR2_outside[i][j] = 0.;
		avgR2_long[i][j] = 0.;
		avgR2_sidelong[i][j] = 0.;
		avgR2_outlong[i][j] = 0.;

		CavgR2_side[i][j] = 0.;
		CavgR2_out[i][j] = 0.;
		CavgR2_outside[i][j] = 0.;
		CavgR2_long[i][j] = 0.;
		CavgR2_sidelong[i][j] = 0.;
		CavgR2_outlong[i][j] = 0.;

		R2_side_err[i][j] = 0.;
		R2_out_err[i][j] = 0.;
		R2_long_err[i][j] = 0.;
		R2_outside_err[i][j] = 0.;
		R2_sidelong_err[i][j] = 0.;
		R2_outlong_err[i][j] = 0.;

		lambda_Correl[i][j] = 0.;
		lambda_Correl_err[i][j] = 0.;
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

		xs_t_cos[i][j] = 0.;
		xo_t_cos[i][j] = 0.;
		xl_t_cos[i][j] = 0.;
		xo_xs_cos[i][j] = 0.;
		xl_xs_cos[i][j] = 0.;
		xo_xl_cos[i][j] = 0.;
		xs2_cos[i][j] = 0.;
		xo2_cos[i][j] = 0.;
		xl2_cos[i][j] = 0.;
		t2_cos[i][j] = 0.;
		xs_t_sin[i][j] = 0.;
		xo_t_sin[i][j] = 0.;
		xl_t_sin[i][j] = 0.;
		xo_xs_sin[i][j] = 0.;
		xl_xs_sin[i][j] = 0.;
		xo_xl_sin[i][j] = 0.;
		xs2_sin[i][j] = 0.;
		xo2_sin[i][j] = 0.;
		xl2_sin[i][j] = 0.;
		t2_sin[i][j] = 0.;

		avgR2_side_C[i][j] = 0.;
		avgR2_side_S[i][j] = 0.;
		avgR2_out_C[i][j] = 0.;
		avgR2_out_S[i][j] = 0.;
		avgR2_outside_C[i][j] = 0.;
		avgR2_outside_S[i][j] = 0.;
		avgR2_long_C[i][j] = 0.;
		avgR2_long_S[i][j] = 0.;
		avgR2_sidelong_C[i][j] = 0.;
		avgR2_sidelong_S[i][j] = 0.;
		avgR2_outlong_C[i][j] = 0.;
		avgR2_outlong_S[i][j] = 0.;

		CavgR2_side_C[i][j] = 0.;
		CavgR2_side_S[i][j] = 0.;
		CavgR2_out_C[i][j] = 0.;
		CavgR2_out_S[i][j] = 0.;
		CavgR2_outside_C[i][j] = 0.;
		CavgR2_outside_S[i][j] = 0.;
		CavgR2_long_C[i][j] = 0.;
		CavgR2_long_S[i][j] = 0.;
		CavgR2_sidelong_C[i][j] = 0.;
		CavgR2_sidelong_S[i][j] = 0.;
		CavgR2_outlong_C[i][j] = 0.;
		CavgR2_outlong_S[i][j] = 0.;
	}
}


   return;
}

void doHBT::Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx)
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
      anisotropic_flows_pTdiff[i][j] = 0.0e0;
      anisotropic_flows_pTdiff_psin[i][j] = 0.0e0;
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
      (*Emissionfunction_ptr)[i].CDF_value = 0.0;
   }

//reset only EBE source variances and EBE HBT radii/coeffs
for(int i=0; i<n_localp_T; i++)
{
   //reset azimuthally averaged radii first
      azavg_squared_S_func[i] = 0.0;
      azavg_S_func[i] = 0.0;

      azavg_Cavg_squared_S_func[i] = 0.0;
      azavg_CavgR2_side_num[i] = 0.0;
      azavg_CavgR2_out_num[i] = 0.0;
      azavg_CavgR2_long_num[i] = 0.0;
      azavg_CavgR2_outside_num[i] = 0.0;
      azavg_CavgR2_sidelong_num[i] = 0.0;
      azavg_CavgR2_outlong_num[i] = 0.0;

      azavg_R2_side[i] = 0.0;
      azavg_R2_out[i] = 0.0;
      azavg_R2_outside[i] = 0.0;
      azavg_R2_long[i] = 0.0;
      azavg_R2_sidelong[i] = 0.0;
      azavg_R2_outlong[i] = 0.0;

      azavg_avgR2_side[i] = 0.0;
      azavg_avgR2_out[i] = 0.0;
      azavg_avgR2_outside[i] = 0.0;
      azavg_avgR2_long[i] = 0.0;
      azavg_avgR2_sidelong[i] = 0.0;
      azavg_avgR2_outlong[i] = 0.0;

      azavg_CavgR2_side[i] = 0.0;
      azavg_CavgR2_out[i] = 0.0;
      azavg_CavgR2_outside[i] = 0.0;
      azavg_CavgR2_long[i] = 0.0;
      azavg_CavgR2_sidelong[i] = 0.0;
      azavg_CavgR2_outlong[i] = 0.0;

	for(int j=0; j<n_localp_phi; j++)
	{
		S_func[i][j] = 0.;
		squared_S_func[i][j] = 0.;
		xs_S[i][j] = 0.;
		xo_S[i][j] = 0.;
		xl_S[i][j] = 0.;
		t_S[i][j] = 0.;
		xs_t_S[i][j] = 0.;
		xo_t_S[i][j] = 0.;
		xl_t_S[i][j] = 0.;
		xo_xs_S[i][j] = 0.;
		xl_xs_S[i][j] = 0.;
		xo_xl_S[i][j] = 0.;
		xs2_S[i][j] = 0.;
		xo2_S[i][j] = 0.;
		xl2_S[i][j] = 0.;
		t2_S[i][j] = 0.;

		avgS_func[i][j] = 0.;
		avgxs_S[i][j] = 0.;
		avgxo_S[i][j] = 0.;
		avgxl_S[i][j] = 0.;
		avgt_S[i][j] = 0.;
		avgxs_t_S[i][j] = 0.;
		avgxo_t_S[i][j] = 0.;
		avgxl_t_S[i][j] = 0.;
		avgxo_xs_S[i][j] = 0.;
		avgxl_xs_S[i][j] = 0.;
		avgxo_xl_S[i][j] = 0.;
		avgxs2_S[i][j] = 0.;
		avgxo2_S[i][j] = 0.;
		avgxl2_S[i][j] = 0.;
		avgt2_S[i][j] = 0.;

		CavgS_func_squared[i][j] = 0.;
		CavgR2_side_num[i][j] = 0.;
		CavgR2_out_num[i][j] = 0.;
		CavgR2_outside_num[i][j] = 0.;
		CavgR2_long_num[i][j] = 0.;
		CavgR2_sidelong_num[i][j] = 0.;
		CavgR2_outlong_num[i][j] = 0.;

		R2_side[i][j] = 0.;
		R2_out[i][j] = 0.;
		R2_outside[i][j] = 0.;
		R2_long[i][j] = 0.;
		R2_sidelong[i][j] = 0.;
		R2_outlong[i][j] = 0.;

		avgR2_side[i][j] = 0.;
		avgR2_out[i][j] = 0.;
		avgR2_outside[i][j] = 0.;
		avgR2_long[i][j] = 0.;
		avgR2_sidelong[i][j] = 0.;
		avgR2_outlong[i][j] = 0.;

		CavgR2_side[i][j] = 0.;
		CavgR2_out[i][j] = 0.;
		CavgR2_outside[i][j] = 0.;
		CavgR2_long[i][j] = 0.;
		CavgR2_sidelong[i][j] = 0.;
		CavgR2_outlong[i][j] = 0.;
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

		xs_t_cos[i][j] = 0.;
		xo_t_cos[i][j] = 0.;
		xl_t_cos[i][j] = 0.;
		xo_xs_cos[i][j] = 0.;
		xl_xs_cos[i][j] = 0.;
		xo_xl_cos[i][j] = 0.;
		xs2_cos[i][j] = 0.;
		xo2_cos[i][j] = 0.;
		xl2_cos[i][j] = 0.;
		t2_cos[i][j] = 0.;
		xs_t_sin[i][j] = 0.;
		xo_t_sin[i][j] = 0.;
		xl_t_sin[i][j] = 0.;
		xo_xs_sin[i][j] = 0.;
		xl_xs_sin[i][j] = 0.;
		xo_xl_sin[i][j] = 0.;
		xs2_sin[i][j] = 0.;
		xo2_sin[i][j] = 0.;
		xl2_sin[i][j] = 0.;
		t2_sin[i][j] = 0.;

		avgR2_side_C[i][j] = 0.;
		avgR2_side_S[i][j] = 0.;
		avgR2_out_C[i][j] = 0.;
		avgR2_out_S[i][j] = 0.;
		avgR2_outside_C[i][j] = 0.;
		avgR2_outside_S[i][j] = 0.;
		avgR2_long_C[i][j] = 0.;
		avgR2_long_S[i][j] = 0.;
		avgR2_sidelong_C[i][j] = 0.;
		avgR2_sidelong_S[i][j] = 0.;
		avgR2_outlong_C[i][j] = 0.;
		avgR2_outlong_S[i][j] = 0.;

		CavgR2_side_C[i][j] = 0.;
		CavgR2_side_S[i][j] = 0.;
		CavgR2_out_C[i][j] = 0.;
		CavgR2_out_S[i][j] = 0.;
		CavgR2_outside_C[i][j] = 0.;
		CavgR2_outside_S[i][j] = 0.;
		CavgR2_long_C[i][j] = 0.;
		CavgR2_long_S[i][j] = 0.;
		CavgR2_sidelong_C[i][j] = 0.;
		CavgR2_sidelong_S[i][j] = 0.;
		CavgR2_outlong_C[i][j] = 0.;
		CavgR2_outlong_S[i][j] = 0.;
	}
}

   return;
}

doHBT::~doHBT()
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
      delete[] anisotropic_flows_pTdiff[i];
      delete[] anisotropic_flows_pTdiff_psin[i];
   }
   delete[] dN_dypTdpTdphi;
   delete[] cosine_iorder;
   delete[] sine_iorder;
   delete[] anisotropic_flows_pTdiff;
   delete[] anisotropic_flows_pTdiff_psin;
   delete[] plane_angle;
   delete[] avgplane_angle;
   delete[] Cavgplane_angle;
   delete[] anisotropic_flows;

   delete[] K_T;
   delete[] K_phi;
   delete[] K_phi_weight;
   delete[] eta_s;
   delete[] ch_eta_s;
   delete[] sh_eta_s;
   delete[] eta_s_weight;

   for(int i=0; i<n_localp_T; i++)
   {
      delete[] R2_side[i];
      delete[] R2_side_C[i];
      delete[] R2_side_S[i];
      delete[] R2_out[i];
      delete[] R2_out_C[i];
      delete[] R2_out_S[i];
      delete[] R2_outside[i];
      delete[] R2_outside_C[i];
      delete[] R2_outside_S[i];
      delete[] R2_long[i];
      delete[] R2_long_C[i];
      delete[] R2_long_S[i];
      delete[] R2_sidelong[i];
      delete[] R2_sidelong_C[i];
      delete[] R2_sidelong_S[i];
      delete[] R2_outlong[i];
      delete[] R2_outlong_C[i];
      delete[] R2_outlong_S[i];

      delete[] R2_side_err[i];
      delete[] R2_out_err[i];
      delete[] R2_outside_err[i];
      delete[] R2_long_err[i];
      delete[] R2_sidelong_err[i];
      delete[] R2_outlong_err[i];

      delete[] lambda_Correl[i];
      delete[] lambda_Correl_err[i];

      delete[] xs_t_cos[i];
      delete[] xo_t_cos[i];
      delete[] xl_t_cos[i];
      delete[] xo_xs_cos[i];
      delete[] xl_xs_cos[i];
      delete[] xo_xl_cos[i];
      delete[] xs2_cos[i];
      delete[] xo2_cos[i];
      delete[] xl2_cos[i];
      delete[] t2_cos[i];
      delete[] xs_t_sin[i];
      delete[] xo_t_sin[i];
      delete[] xl_t_sin[i];
      delete[] xo_xs_sin[i];
      delete[] xl_xs_sin[i];
      delete[] xo_xl_sin[i];
      delete[] xs2_sin[i];
      delete[] xo2_sin[i];
      delete[] xl2_sin[i];
      delete[] t2_sin[i];

      delete[] avgR2_side[i];
      delete[] avgR2_side_C[i];
      delete[] avgR2_side_S[i];
      delete[] avgR2_out[i];
      delete[] avgR2_out_C[i];
      delete[] avgR2_out_S[i];
      delete[] avgR2_outside[i];
      delete[] avgR2_outside_C[i];
      delete[] avgR2_outside_S[i];
      delete[] avgR2_long[i];
      delete[] avgR2_long_C[i];
      delete[] avgR2_long_S[i];
      delete[] avgR2_sidelong[i];
      delete[] avgR2_sidelong_C[i];
      delete[] avgR2_sidelong_S[i];
      delete[] avgR2_outlong[i];
      delete[] avgR2_outlong_C[i];
      delete[] avgR2_outlong_S[i];

      delete[] CavgR2_side[i];
      delete[] CavgR2_side_C[i];
      delete[] CavgR2_side_S[i];
      delete[] CavgR2_out[i];
      delete[] CavgR2_out_C[i];
      delete[] CavgR2_out_S[i];
      delete[] CavgR2_outside[i];
      delete[] CavgR2_outside_C[i];
      delete[] CavgR2_outside_S[i];
      delete[] CavgR2_long[i];
      delete[] CavgR2_long_C[i];
      delete[] CavgR2_long_S[i];
      delete[] CavgR2_sidelong[i];
      delete[] CavgR2_sidelong_C[i];
      delete[] CavgR2_sidelong_S[i];
      delete[] CavgR2_outlong[i];
      delete[] CavgR2_outlong_C[i];
      delete[] CavgR2_outlong_S[i];
   }

   delete[] R2_side;
   delete[] R2_side_C;
   delete[] R2_side_S;
   delete[] R2_out;
   delete[] R2_out_C;
   delete[] R2_out_S;
   delete[] R2_outside;
   delete[] R2_outside_C;
   delete[] R2_outside_S;
   delete[] R2_long;
   delete[] R2_long_C;
   delete[] R2_long_S;
   delete[] R2_sidelong;
   delete[] R2_sidelong_C;
   delete[] R2_sidelong_S;
   delete[] R2_outlong;
   delete[] R2_outlong_C;
   delete[] R2_outlong_S;

   delete[] azavg_R2_side;
   delete[] azavg_R2_out;
   delete[] azavg_R2_outside;
   delete[] azavg_R2_long;
   delete[] azavg_R2_sidelong;
   delete[] azavg_R2_outlong;

   delete[] R2_side_err;
   delete[] R2_out_err;
   delete[] R2_outside_err;
   delete[] R2_long_err;
   delete[] R2_sidelong_err;
   delete[] R2_outlong_err;

   delete[] lambda_Correl;
   delete[] lambda_Correl_err;

   delete[] xs_t_cos;
   delete[] xo_t_cos;
   delete[] xl_t_cos;
   delete[] xo_xs_cos;
   delete[] xl_xs_cos;
   delete[] xo_xl_cos;
   delete[] xs2_cos;
   delete[] xo2_cos;
   delete[] xl2_cos;
   delete[] t2_cos;
   delete[] xs_t_sin;
   delete[] xo_t_sin;
   delete[] xl_t_sin;
   delete[] xo_xs_sin;
   delete[] xl_xs_sin;
   delete[] xo_xl_sin;
   delete[] xs2_sin;
   delete[] xo2_sin;
   delete[] xl2_sin;
   delete[] t2_sin;

   delete[] avgR2_side;
   delete[] avgR2_side_C;
   delete[] avgR2_side_S;
   delete[] avgR2_out;
   delete[] avgR2_out_C;
   delete[] avgR2_out_S;
   delete[] avgR2_outside;
   delete[] avgR2_outside_C;
   delete[] avgR2_outside_S;
   delete[] avgR2_long;
   delete[] avgR2_long_C;
   delete[] avgR2_long_S;
   delete[] avgR2_sidelong;
   delete[] avgR2_sidelong_C;
   delete[] avgR2_sidelong_S;
   delete[] avgR2_outlong;
   delete[] avgR2_outlong_C;
   delete[] avgR2_outlong_S;

   delete[] azavg_avgR2_side;
   delete[] azavg_avgR2_out;
   delete[] azavg_avgR2_outside;
   delete[] azavg_avgR2_long;
   delete[] azavg_avgR2_sidelong;
   delete[] azavg_avgR2_outlong;

   delete[] CavgR2_side;
   delete[] CavgR2_side_C;
   delete[] CavgR2_side_S;
   delete[] CavgR2_out;
   delete[] CavgR2_out_C;
   delete[] CavgR2_out_S;
   delete[] CavgR2_outside;
   delete[] CavgR2_outside_C;
   delete[] CavgR2_outside_S;
   delete[] CavgR2_long;
   delete[] CavgR2_long_C;
   delete[] CavgR2_long_S;
   delete[] CavgR2_sidelong;
   delete[] CavgR2_sidelong_C;
   delete[] CavgR2_sidelong_S;
   delete[] CavgR2_outlong;
   delete[] CavgR2_outlong_C;
   delete[] CavgR2_outlong_S;

   delete[] azavg_CavgR2_side;
   delete[] azavg_CavgR2_out;
   delete[] azavg_CavgR2_outside;
   delete[] azavg_CavgR2_long;
   delete[] azavg_CavgR2_sidelong;
   delete[] azavg_CavgR2_outlong;

   for(int i=0; i<n_localp_T; i++)
   {
      delete[] S_func[i];
      delete[] squared_S_func[i];
      delete[] xs_S[i];
      delete[] xo_S[i];
      delete[] xl_S[i];
      delete[] t_S[i];
      delete[] xs_t_S[i];
      delete[] xo_t_S[i];
      delete[] xl_t_S[i];
      delete[] xo_xs_S[i];
      delete[] xl_xs_S[i];
      delete[] xo_xl_S[i];
      delete[] xs2_S[i];
      delete[] xo2_S[i];
      delete[] xl2_S[i];
      delete[] t2_S[i];

      delete[] avgS_func[i];
      delete[] avgxs_S[i];
      delete[] avgxo_S[i];
      delete[] avgxl_S[i];
      delete[] avgt_S[i];
      delete[] avgxs_t_S[i];
      delete[] avgxo_t_S[i];
      delete[] avgxl_t_S[i];
      delete[] avgxo_xs_S[i];
      delete[] avgxl_xs_S[i];
      delete[] avgxo_xl_S[i];
      delete[] avgxs2_S[i];
      delete[] avgxo2_S[i];
      delete[] avgxl2_S[i];
      delete[] avgt2_S[i];

      delete[] CavgS_func_squared[i];
      delete[] CavgR2_side_num[i];
      delete[] CavgR2_out_num[i];
      delete[] CavgR2_outside_num[i];
      delete[] CavgR2_long_num[i];
      delete[] CavgR2_sidelong_num[i];
      delete[] CavgR2_outlong_num[i];
   }
   delete[] S_func;
   delete[] squared_S_func;
   delete[] xs_S;
   delete[] xo_S;
   delete[] xl_S;
   delete[] t_S;
   delete[] xs_t_S;
   delete[] xo_t_S;
   delete[] xl_t_S;
   delete[] xo_xs_S;
   delete[] xl_xs_S;
   delete[] xo_xl_S;
   delete[] xs2_S;
   delete[] xo2_S;
   delete[] xl2_S;
   delete[] t2_S;

   delete[] avgS_func;
   delete[] avgxs_S;
   delete[] avgxo_S;
   delete[] avgxl_S;
   delete[] avgt_S;
   delete[] avgxs_t_S;
   delete[] avgxo_t_S;
   delete[] avgxl_t_S;
   delete[] avgxo_xs_S;
   delete[] avgxl_xs_S;
   delete[] avgxo_xl_S;
   delete[] avgxs2_S;
   delete[] avgxo2_S;
   delete[] avgxl2_S;
   delete[] avgt2_S;

   delete[] CavgS_func_squared;
   delete[] CavgR2_side_num;
   delete[] CavgR2_out_num;
   delete[] CavgR2_outside_num;
   delete[] CavgR2_long_num;
   delete[] CavgR2_sidelong_num;
   delete[] CavgR2_outlong_num;

   delete[] azavg_Cavg_squared_S_func;
   delete[] azavg_CavgR2_side_num;
   delete[] azavg_CavgR2_out_num;
   delete[] azavg_CavgR2_outside_num;
   delete[] azavg_CavgR2_long_num;
   delete[] azavg_CavgR2_sidelong_num;
   delete[] azavg_CavgR2_outlong_num;

   delete[] q_out;
   delete[] q_side;
   delete[] q_long;

   delete[] Correl_1D_out;
   delete[] Correl_1D_side;
   delete[] Correl_1D_long;
   delete[] Correl_1D_out_err;
   delete[] Correl_1D_side_err;
   delete[] Correl_1D_long_err;

   for(int i=0; i<qnpts; i++)
   {
      for(int j=0; j< qnpts; j++)
          delete[] Correl_3D[i][j];
      delete[] Correl_3D[i];
      for(int j=0; j< qnpts; j++)
          delete[] Correl_3D_err[i][j];
      delete[] Correl_3D_err[i];
   }
   delete[] Correl_3D;
   delete[] Correl_3D_err;

   return;
}

void doHBT::Reset_EmissionData()
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
      (*Emissionfunction_ptr)[i].CDF_value = 0.0;
   }
}

void doHBT::Reset_source_variances_and_HBT_radii()
{
	//reset only EBE source variances and EBE HBT radii/coeffs
	for(int i=0; i<n_localp_T; i++)
	{
		//reset azimuthally averaged radii first
		
		azavg_Cavg_squared_S_func[i] = 0.0;
		azavg_CavgR2_side_num[i] = 0.0;
		azavg_CavgR2_out_num[i] = 0.0;
		azavg_CavgR2_long_num[i] = 0.0;
		azavg_CavgR2_outside_num[i] = 0.0;
		azavg_CavgR2_sidelong_num[i] = 0.0;
		azavg_CavgR2_outlong_num[i] = 0.0;
		
		azavg_R2_side[i] = 0.0;
		azavg_R2_out[i] = 0.0;
		azavg_R2_outside[i] = 0.0;
		azavg_R2_long[i] = 0.0;
		azavg_R2_sidelong[i] = 0.0;
		azavg_R2_outlong[i] = 0.0;
		
		azavg_avgR2_side[i] = 0.0;
		azavg_avgR2_out[i] = 0.0;
		azavg_avgR2_outside[i] = 0.0;
		azavg_avgR2_long[i] = 0.0;
		azavg_avgR2_sidelong[i] = 0.0;
		azavg_avgR2_outlong[i] = 0.0;
		
		azavg_CavgR2_side[i] = 0.0;
		azavg_CavgR2_out[i] = 0.0;
		azavg_CavgR2_outside[i] = 0.0;
		azavg_CavgR2_long[i] = 0.0;
		azavg_CavgR2_sidelong[i] = 0.0;
		azavg_CavgR2_outlong[i] = 0.0;

		for(int j=0; j<n_localp_phi; j++)
		{
			S_func[i][j] = 0.;
			squared_S_func[i][j] = 0.;
			xs_S[i][j] = 0.;
			xo_S[i][j] = 0.;
			xl_S[i][j] = 0.;
			t_S[i][j] = 0.;
			xs_t_S[i][j] = 0.;
			xo_t_S[i][j] = 0.;
			xl_t_S[i][j] = 0.;
			xo_xs_S[i][j] = 0.;
			xl_xs_S[i][j] = 0.;
			xo_xl_S[i][j] = 0.;
			xs2_S[i][j] = 0.;
			xo2_S[i][j] = 0.;
			xl2_S[i][j] = 0.;
			t2_S[i][j] = 0.;
	
			avgS_func[i][j] = 0.;
			avgxs_S[i][j] = 0.;
			avgxo_S[i][j] = 0.;
			avgxl_S[i][j] = 0.;
			avgt_S[i][j] = 0.;
			avgxs_t_S[i][j] = 0.;
			avgxo_t_S[i][j] = 0.;
			avgxl_t_S[i][j] = 0.;
			avgxo_xs_S[i][j] = 0.;
			avgxl_xs_S[i][j] = 0.;
			avgxo_xl_S[i][j] = 0.;
			avgxs2_S[i][j] = 0.;
			avgxo2_S[i][j] = 0.;
			avgxl2_S[i][j] = 0.;
			avgt2_S[i][j] = 0.;
	
			CavgS_func_squared[i][j] = 0.;
			CavgR2_side_num[i][j] = 0.;
			CavgR2_out_num[i][j] = 0.;
			CavgR2_outside_num[i][j] = 0.;
			CavgR2_long_num[i][j] = 0.;
			CavgR2_sidelong_num[i][j] = 0.;
			CavgR2_outlong_num[i][j] = 0.;
	
			R2_side[i][j] = 0.;
			R2_out[i][j] = 0.;
			R2_outside[i][j] = 0.;
			R2_long[i][j] = 0.;
			R2_sidelong[i][j] = 0.;
			R2_outlong[i][j] = 0.;
	
			avgR2_side[i][j] = 0.;
			avgR2_out[i][j] = 0.;
			avgR2_outside[i][j] = 0.;
			avgR2_long[i][j] = 0.;
			avgR2_sidelong[i][j] = 0.;
			avgR2_outlong[i][j] = 0.;
	
			CavgR2_side[i][j] = 0.;
			CavgR2_out[i][j] = 0.;
			CavgR2_outside[i][j] = 0.;
			CavgR2_long[i][j] = 0.;
			CavgR2_sidelong[i][j] = 0.;
			CavgR2_outlong[i][j] = 0.;
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
	
			xs_t_cos[i][j] = 0.;
			xo_t_cos[i][j] = 0.;
			xl_t_cos[i][j] = 0.;
			xo_xs_cos[i][j] = 0.;
			xl_xs_cos[i][j] = 0.;
			xo_xl_cos[i][j] = 0.;
			xs2_cos[i][j] = 0.;
			xo2_cos[i][j] = 0.;
			xl2_cos[i][j] = 0.;
			t2_cos[i][j] = 0.;
			xs_t_sin[i][j] = 0.;
			xo_t_sin[i][j] = 0.;
			xl_t_sin[i][j] = 0.;
			xo_xs_sin[i][j] = 0.;
			xl_xs_sin[i][j] = 0.;
			xo_xl_sin[i][j] = 0.;
			xs2_sin[i][j] = 0.;
			xo2_sin[i][j] = 0.;
			xl2_sin[i][j] = 0.;
			t2_sin[i][j] = 0.;
	
			avgR2_side_C[i][j] = 0.;
			avgR2_side_S[i][j] = 0.;
			avgR2_out_C[i][j] = 0.;
			avgR2_out_S[i][j] = 0.;
			avgR2_outside_C[i][j] = 0.;
			avgR2_outside_S[i][j] = 0.;
			avgR2_long_C[i][j] = 0.;
			avgR2_long_S[i][j] = 0.;
			avgR2_sidelong_C[i][j] = 0.;
			avgR2_sidelong_S[i][j] = 0.;
			avgR2_outlong_C[i][j] = 0.;
			avgR2_outlong_S[i][j] = 0.;
	
			CavgR2_side_C[i][j] = 0.;
			CavgR2_side_S[i][j] = 0.;
			CavgR2_out_C[i][j] = 0.;
			CavgR2_out_S[i][j] = 0.;
			CavgR2_outside_C[i][j] = 0.;
			CavgR2_outside_S[i][j] = 0.;
			CavgR2_long_C[i][j] = 0.;
			CavgR2_long_S[i][j] = 0.;
			CavgR2_sidelong_C[i][j] = 0.;
			CavgR2_sidelong_S[i][j] = 0.;
			CavgR2_outlong_C[i][j] = 0.;
			CavgR2_outlong_S[i][j] = 0.;
		}
	}
}


bool doHBT::fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

void doHBT::debugger()
{
	cerr << "You made it to checkpoint #" << checkpoint_index << "!" << endl;
	checkpoint_index++;
	return;
}


//End of file

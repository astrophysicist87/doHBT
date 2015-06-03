#ifndef doHBT_H
#define doHBT_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>
#include<vector>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "readindata.h"
#include "parameters.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

typedef struct
{
   double t, x, y, z;
   double r, phi;
   double tau, eta;
   double data;
   double CDF_value;
}Emissionfunction_data;

struct Correlationfunction3D_data
{
  size_t data_length;
  double *q_o;
  double *q_s;
  double *q_l;
  double *y;
  double *sigma;
};

struct Correlationfunction1D_data
{
  size_t data_length;
  double *q;
  double *y;
  double *sigma;
};

int Fittarget_correlfun1D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun1D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun1D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);
int Fittarget_correlfun3D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun3D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun3D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);
int Fittarget_correlfun3D_f_withlambda (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun3D_df_withlambda (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun3D_fdf_withlambda (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);

class doHBT
{
   private:
      //particle information 
      string particle_name;
      double particle_mass;
      int particle_id;     //particle id
      double particle_sign;   //+/- 1 for Fermi/Bose statistics for baryon/meson
      double particle_gspin;  //particle degeneracy       

      //single particle spectra for plane angle determination
      //int n_order;
      double* SP_pT;
      double* SP_pphi;
      double* SP_pT_weight;
      double* SP_pphi_weight;
      double SP_p_y;
      double** dN_dypTdpTdphi;
      double* dN_dypTdpT;
      double** cosine_iorder;
      double** sine_iorder;
      double** EdNd3p_cfs;
      double** EdNd3p_phases;
      double* dN_dydphi;
      double* pTdN_dydphi;
      double* plane_angle;
      double* avgplane_angle;
      double* Cavgplane_angle;
      double* anisotropic_flows;
      double** anisotropic_flows_pTdiff;
      double** anisotropic_flows_pTdiff_psin;
      double global_plane_psi;
      double mean_pT;
      double ** eventByEventMultiplicities;
     
      //pair momentum
      double K_y;
      double* K_T;
      double* K_phi;
      double* K_phi_weight;
      double beta_perp, beta_l;
    
      //spatial rapidity grid
      double* eta_s;
      double* ch_eta_s;
      double* sh_eta_s;
      double* eta_s_weight;

      //Emission function
      vector<Emissionfunction_data>* Emissionfunction_ptr;
      int FO_length;
      int Emissionfunction_length;
      vector<Emissionfunction_data>* avgFOsurf_ptr;

      double spectra;

      double* q_out;
      double* q_side;
      double* q_long;

      //store correlation functions
      double* Correl_1D_out;
      double* Correl_1D_out_err;
      double* Correl_1D_side;
      double* Correl_1D_side_err;
      double* Correl_1D_long;
      double* Correl_1D_long_err;
      double*** Correl_3D;
      double*** Correl_3D_err;


      //source variances
      double **S_func, **squared_S_func;
      double **xs_S, **xo_S, **xl_S, **t_S;
      double **xs_t_S, **xo_t_S, **xl_t_S;
      double **xo_xs_S, **xl_xs_S, **xo_xl_S;
      double **xs2_S, **xo2_S, **xl2_S, **t2_S;
      double **xs_t_cos, **xo_t_cos, **xl_t_cos;
      double **xo_xs_cos, **xl_xs_cos, **xo_xl_cos;
      double **xs2_cos, **xo2_cos, **xl2_cos, **t2_cos;
      double **xs_t_sin, **xo_t_sin, **xl_t_sin;
      double **xo_xs_sin, **xl_xs_sin, **xo_xl_sin;
      double **xs2_sin, **xo2_sin, **xl2_sin, **t2_sin;

      //azimuthally averaged source variances and HBT radii
      double *azavg_S_func, *azavg_squared_S_func;
      double *azavg_R2_side, *azavg_R2_out, *azavg_R2_long, *azavg_R2_outside, *azavg_R2_sidelong, *azavg_R2_outlong;

      //azimuthally averaged, EBE avgd. source variances and HBT radii
      double *azavg_avgS_func, *azavg_squared_avgS_func;
      double *azavg_avgR2_side, *azavg_avgR2_out, *azavg_avgR2_long, *azavg_avgR2_outside, *azavg_avgR2_sidelong, *azavg_avgR2_outlong;
      double *azavg_Cavg_squared_S_func;
      double *azavg_CavgR2_side_num, *azavg_CavgR2_out_num, *azavg_CavgR2_long_num;
      double *azavg_CavgR2_outside_num, *azavg_CavgR2_sidelong_num, *azavg_CavgR2_outlong_num;
      double *azavg_CavgR2_side, *azavg_CavgR2_out, *azavg_CavgR2_long, *azavg_CavgR2_outside, *azavg_CavgR2_sidelong, *azavg_CavgR2_outlong;

      //EBE avgd. source variances
      double **avgS_func;
      double **avgxs_S, **avgxo_S, **avgxl_S, **avgt_S;
      double **avgxs_t_S, **avgxo_t_S, **avgxl_t_S;
      double **avgxo_xs_S, **avgxl_xs_S, **avgxo_xl_S;
      double **avgxs2_S, **avgxo2_S, **avgxl2_S, **avgt2_S;

      double **CavgS_func_squared;
      double **CavgR2_side_num, **CavgR2_out_num, **CavgR2_long_num;
      double **CavgR2_outside_num, **CavgR2_sidelong_num, **CavgR2_outlong_num;


      //HBT radii coefficients
      double **R2_side, **R2_out, **R2_long, **R2_outside, **R2_sidelong, **R2_outlong;
      double **R2_side_err, **R2_out_err, **R2_long_err, **R2_outside_err, **R2_sidelong_err, **R2_outlong_err;
      double **R2_side_C, **R2_side_S;
      double **R2_out_C, **R2_out_S;
      double **R2_long_C, **R2_long_S;
      double **R2_outside_C, **R2_outside_S;
      double **R2_outlong_C, **R2_outlong_S;
      double **R2_sidelong_C, **R2_sidelong_S;

      double **lambda_Correl;
      double **lambda_Correl_err;
      //double **lambda_Correl_C, **lambda_Correl_S;

      //HBT radii coefficients from event-averaged emission function Sbar
      double **avgR2_side, **avgR2_out, **avgR2_long, **avgR2_outside, **avgR2_sidelong, **avgR2_outlong;
      double **avgR2_side_C, **avgR2_side_S;
      double **avgR2_out_C, **avgR2_out_S;
      double **avgR2_long_C, **avgR2_long_S;
      double **avgR2_outside_C, **avgR2_outside_S;
      double **avgR2_outlong_C, **avgR2_outlong_S;
      double **avgR2_sidelong_C, **avgR2_sidelong_S;

      //HBT radii coefficients from curvature of <C>_{ev} (i.e., source variances)
      double **CavgR2_side, **CavgR2_out, **CavgR2_long, **CavgR2_outside, **CavgR2_sidelong, **CavgR2_outlong;
      double **CavgR2_side_C, **CavgR2_side_S;
      double **CavgR2_out_C, **CavgR2_out_S;
      double **CavgR2_long_C, **CavgR2_long_S;
      double **CavgR2_outside_C, **CavgR2_outside_S;
      double **CavgR2_outlong_C, **CavgR2_outlong_S;
      double **CavgR2_sidelong_C, **CavgR2_sidelong_S;

      //miscellaneous
      ofstream * global_out_stream_ptr;
      int global_folderindex;
      string global_path;
      string global_runfolder;
      string global_resultsfolder_stem;
      string no_df_stem;
	vector<int> eventvector;
      int checkpoint_index;

   public:
      doHBT();
      ~doHBT();

      void Determine_plane_angle(FO_surf* FOsurf_ptr);
      void Get_EdNd3p_cfs(FO_surf* FOsurf_ptr, bool determine_plane_angle = true);
      void Get_azimuthally_averaged_multiplicities();
      void Get_azimuthally_averaged_EBE_HBTradii();
      void Determine_avgplane_angle();
      void Determine_Cavgplane_angle();
      void Analyze_sourcefunction(FO_surf* FOsurf_ptr);
      void quick_Analyze_sourcefunction();
      void quick_Analyze_sourcefunction_vars();
      void Analyze_AVG_sourcefunction();
      void Analyze_CAVG_sourcefunction();
      void Analyze_azimuthally_averaged_AVG_sourcefunction();
      void Analyze_azimuthally_averaged_CAVG_sourcefunction();
      void Reset_EmissionData();
      void Reset_source_variances_and_HBT_radii();
      void Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx);
      void SetEmissionData(FO_surf* FOsurf_ptr, double K_T_local, double K_phi_local, bool includezeroes);
      void Output_Emission_Function(int iKT, int iKphi, int folderindex);
      double Average_Emission_Function_on_FOsurface(FO_surf* FOsurf_ptr, int FOcell, int iKT);
      void Average_sourcefunction_on_FOsurface(FO_surf* FOsurf_ptr);
      double Average_Emission_Function_on_FOsurface(FO_surf* FOsurf_ptr, int FOcell, int iKT, int iKphi);
      void Average_sourcefunction_on_FOsurface(FO_surf* FOsurf_ptr, int iKphi);
      void Output_avgEmission_Function_on_FOsurface(int folderindex);
      void Update_avgSource_function(int iKT = -1, int iKphi = -1);
      void Update_CavgSource_function(int iKT = -1, int iKphi = -1);
      void Calculate_avgSource_function(int i = -1, int j = -1);
      void Calculate_CavgSource_function(int i = -1, int j = -1);
      //void Update_avgSource_function(int iKT = -1);		//the azimuthally independent versions
      void Update_azavgCavgSource_function(int iKT = -1);		//the azimuthally independent versions
      //void Calculate_avgSource_function(int i = -1);		//the azimuthally independent versions
      void Calculate_azavgCavgSource_function(int i = -1);		//the azimuthally independent versions
      bool fexists(const char *filename);

      void Cal_dN_dypTdpTdphi(double** SP_p0, double** SP_px, double** SP_py, double** SP_pz, FO_surf* FOsurf_ptr);
      double Emissionfunction(double p0, double px, double py, double pz, FO_surf* surf);

      void R2_Fourier_transform(int, double);
      void avgR2_Fourier_transform(int, double);
      void CavgR2_Fourier_transform(int, double);
      void Svars_Fourier_transform(int iKT, double plane_psi);
      void Get_source_variances(int, int);
      void Calculate_HBTradii_from_C_ev(int, int);
      void Get_HBTradii_from_Cbar_and_Cavg();
      void Get_HBTradii_from_azimuthally_averaged_Cavg_and_Cbar();
      void Get_HBTradii_from_Cbar_and_Cavg_random(int);
      void Get_HBTradii_from_azimuthally_averaged_Cavg_random(int);
      void Get_HBTradii_from_Cbar_and_Cavg_fixed_events(vector<int> list_of_events);
      void Simulate_subensemble_averaging();
      void Simulate_azimuthally_averaged_subensemble_averaging();

      //input and output
      void Readin_results(int);
      void Readin_HBTev_results_only(int);
      void Readin_ev_plane_psi(int);
      void Readin_AVG_results();
      void Read_in_event_multiplicity(int);
      void Output_event_multiplicity(int);
      void Output_ev_plane_psi(int);
      void Output_results(int);
      void Output_GF_results(int);
      void Output_Svars_results(int);
      void Output_azimuthally_averaged_results(int folderindex);
      void Readin_azimuthally_averaged_results(int folderindex);
      void Output_AVG_results();
      void Output_CAVG_results();
      void Output_azimuthally_averaged_Cbar_and_CAVG_results();
      void Update_subensemble_indexfile(int iM, int ibin);
      void Output_CAVG_random_results(int iM, int ibin);
      void Output_azimuthally_averaged_CAVG_random_results(int iM, int ibin);
      void Output_HBTcfsev_results_only(int);
      void Output_dN_dypTdpTdphi(int folderindex);
      void Output_dN_dypTdpT(int folderindex);
      void Output_ev_plane_psis(int folderindex);
      void Output_avgplane_psis();
      void Output_Cavgplane_psis();
      void Output_ev_anisotropic_flows(int folderindex);
      void Output_ev_anisotropic_flows_pTdiff(int folderindex);
      void Output_ev_EdNd3p_cfs(int folderindex);
      void Output_ev_mean_pT(int folderindex);

      //functions to calcuate R2ij, various avgd. R2ij, etc.
      void Calculate_R2_side(int, int);
      void Calculate_R2_out(int, int);
      void Calculate_R2_outside(int, int);
      void Calculate_R2_long(int, int);
      void Calculate_R2_sidelong(int, int);
      void Calculate_R2_outlong(int, int);

      void Calculate_avgR2_side(int, int);
      void Calculate_avgR2_out(int, int);
      void Calculate_avgR2_outside(int, int);
      void Calculate_avgR2_long(int, int);
      void Calculate_avgR2_sidelong(int, int);
      void Calculate_avgR2_outlong(int, int);

      void Calculate_CavgR2_side(int, int);
      void Calculate_CavgR2_out(int, int);
      void Calculate_CavgR2_outside(int, int);
      void Calculate_CavgR2_long(int, int);
      void Calculate_CavgR2_sidelong(int, int);
      void Calculate_CavgR2_outlong(int, int);

      void Calculate_azimuthally_averaged_squared_S_func(int iKT);
      void Calculate_azimuthally_averaged_squared_avgS_func(int iKT);
      void Calculate_azimuthally_averaged_Cavg_squared_S_func(int iKT);

      void Calculate_azimuthally_averaged_R2_side(int);
      void Calculate_azimuthally_averaged_R2_out(int);
      void Calculate_azimuthally_averaged_R2_outside(int);
      void Calculate_azimuthally_averaged_R2_long(int);
      void Calculate_azimuthally_averaged_R2_sidelong(int);
      void Calculate_azimuthally_averaged_R2_outlong(int);

      void Calculate_azimuthally_averaged_avgR2_side(int);
      void Calculate_azimuthally_averaged_avgR2_out(int);
      void Calculate_azimuthally_averaged_avgR2_outside(int);
      void Calculate_azimuthally_averaged_avgR2_long(int);
      void Calculate_azimuthally_averaged_avgR2_sidelong(int);
      void Calculate_azimuthally_averaged_avgR2_outlong(int);

      void Calculate_azimuthally_averaged_CavgR2_side(int);
      void Calculate_azimuthally_averaged_CavgR2_out(int);
      void Calculate_azimuthally_averaged_CavgR2_outside(int);
      void Calculate_azimuthally_averaged_CavgR2_long(int);
      void Calculate_azimuthally_averaged_CavgR2_sidelong(int);
      void Calculate_azimuthally_averaged_CavgR2_outlong(int);

//correlation function stuff
      void Cal_correlationfunction_1D(int, int);
      void Cal_correlationfunction_3D(int, int);
      int Read_correlationfunction_1D(int, int);
      int Read_correlationfunction_3D(int, int);
      void Output_Correlationfunction_1D(int, int, int);
      void Output_Correlationfunction_3D(int, int, int);
      void Fit_Correlationfunction1D(char, int, int);
      void Fit_Correlationfunction3D(int, int);
      void Fit_Correlationfunction3D_withlambda(int, int);
      int print_fit_state_1D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
      int print_fit_state_3D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
      int print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
      inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr);
      inline double get_fit_err (int i, gsl_matrix * covariance_ptr);
      void Get_GF_HBTradii(FO_surf* FOsurf_ptr, int folderindex);

//miscellaneous
      void Set_ofstream(ofstream& myout);
      void Set_path(string path);
      void Set_resultsfolder_stem(string usrdef_stem);
      void Set_runfolder(string runfolder);
      void Set_use_delta_f(bool usrdef_usedeltaf);
      void Set_particle_mass(double usrdef_particle_mass);
      void debugger();
      void dummy_function(int i = -1, int j = -1);

      //parameters that the user is free to define
      double plumberg_test_variable;
      bool use_delta_f;
      bool append_output;
      int n_events;
      int total_Nev;	//probably 1000 or something
      int initial_event;
      int MmaxSEA;	//maximum number of times to sub-ensemble average
};

#endif

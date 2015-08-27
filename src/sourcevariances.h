#ifndef SOURCEVARIANCES_H
#define SOURCEVARIANCES_H

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
#include <gsl/gsl_rng.h>// gsl random number generators
#include <gsl/gsl_randist.h>  // gsl random number distributions
#include <gsl/gsl_vector.h>   // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>     // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "SV_readindata.h"
#include "SV_parameters.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

typedef struct
{
   double t, x, y, z;
   double r, phi;
   double tau, eta;
   double data;
}Emissionfunction_data;

typedef struct
{
   int * resonance_particle_id;		// keeps track of current resonance's index in all_particles array
   int * resonance_idx;			// keeps track of current resonance's index in chosen_resonances vector
   int * nbody;
   int * resonance_sign;
   double * resonance_mass;
   double * resonance_mu;
   double * resonance_gspin;
   double * resonance_Gamma;
   double * resonance_total_br;
   double * resonance_direct_br;
   double ** resonance_decay_masses;
   double ** resonance_decay_monvals;
   double ** resonance_decay_Gammas;
   string * resonance_name;
   bool * include_channel;
}decay_info;

struct Correlationfunction1D_data
{
  size_t data_length;
  double * q;
  double * y;
  double * sigma;
};

int Fittarget_correlfun1D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun1D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun1D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);


class SourceVariances
{
	private:
		//particle information 
		string particle_name;
		double particle_mass;
		int particle_monval;
		int particle_id;     //particle id
		double particle_sign;   //+/- 1 for Fermi/Bose statistics for baryon/meson
		double particle_gspin;  //particle degeneracy 
		double particle_mu;
		particle_info * all_particles;
		vector<int> chosen_resonances;
		bool thermal_pions_only;
		int Nparticle;
		int target_particle_id;		//the particle whose spectra (with resonance contributions) you want to compute
		int current_level_of_output;

		int n_zeta_pts, n_v_pts, n_s_pts;
		double v_min, v_max, zeta_min, zeta_max, s_min, s_max;
		
		//integrated space-time moments vars.
		double S_r_to_part, S_xs_r_to_part, S_xo_r_to_part, S_xl_r_to_part, S_t_r_to_part;
		double S_xs2_r_to_part, S_xo2_r_to_part, S_xl2_r_to_part, S_t2_r_to_part;
		double S_xs_xo_r_to_part, S_xs_xl_r_to_part, S_xs_t_r_to_part, S_xo_xl_r_to_part, S_xo_t_r_to_part, S_xl_t_r_to_part;
		
		//array to hold previous and current resonance info
		decay_info decay_channels;
		int current_decay_channel_idx, current_resonance_particle_id, previous_resonance_idx, current_resonance_idx, current_reso_nbody;
		double current_resonance_mu, current_resonance_mass, current_resonance_Gamma, current_m2_Gamma, current_m3_Gamma;
		double current_resonance_total_br, current_resonance_direct_br, current_daughter_mass, current_daughter_Gamma;
		double * current_resonance_decay_masses, * P_eval, * alpha_mu;
		int previous_decay_channel_idx, previous_resonance_particle_id, previous_reso_nbody;
		double previous_resonance_mu, previous_resonance_mass, previous_resonance_Gamma, previous_m2_Gamma, previous_m3_Gamma;
		double previous_resonance_total_br, previous_resonance_direct_br, previous_daughter_mass, previous_daughter_Gamma;
		double * previous_resonance_decay_masses;

		//*************************************************************
		//freeze-out surface interpolation arrays (prefix: "FOI_"), etc...
		int FOI_np0pts, FOI_npTpts, FOI_npphipts, FOI_npzpts, FOI_nmupts; 
		int FOI_netaspts, FOI_nMpts; 
		double * FOI_p0, * FOI_pT, * FOI_pphi, * FOI_sin_pphi, * FOI_cos_pphi;
		double * FOI_pz, * FOI_mu, * FOI_eta_s, * FOI_M;
		double ****** FOI_source_variances;
		//*************************************************************
		
		//array to temporarily hold results of resonance SV integrations
		int n_weighting_functions;  //number of source variances to consider
		//double * source_variances_array;
		double **** integrated_spacetime_moments;
		double **** dN_dypTdpTdphi_moments;
		double **** ln_dN_dypTdpTdphi_moments;
		double **** sign_of_dN_dypTdpTdphi_moments;
		double **** all_particles_dN_dypTdpTdphi_moments;
	
		//needed for resonance calculations
		//kinematic info
		double pstar, Estar, Yp, Ym, DeltaY;
		double MTbar, DeltaMT, MTp, MTm;
		double Qfunc;
		//pair momentum info, currently assumes pT != 0
		double p_y, pT, pphi, mT, mass;
		double ch_p_y, sh_p_y;
		//resonance momentum info
		double P_Y, PT, PPhi, MT, Mres, PPhip, PPhim;
		double m2, m3, Gamma, br, m2Gamma, m3Gamma;
		double * Pp;
		double * Pm;
		double * zvec;

		//SP momentum arrays for interpolation grid
		double * SPinterp_pT;
		double * SPinterp_pphi;
		double * sin_SPinterp_pphi;
		double * cos_SPinterp_pphi;
		double ** SPinterp_p0;
		double ** SPinterp_pz;

		FO_surf* current_FOsurf_ptr;
		//FO surface info that is constant - to save time
		double Tdec, Edec, Pdec, muRES, signRES, gRES;
		double S_prefactor;
	
		//single particle spectra for plane angle determination
		//int n_order;
		double SP_p_y, global_plane_psi, mean_pT;
		double * SP_pT, * SP_pphi, * SP_pT_weight, * SP_pphi_weight, * dN_dypTdpT, * SV_dN_dypTdpT;
		double * dN_dydphi, * SV_dN_dydphi, * pTdN_dydphi, * SV_pTdN_dydphi, * plane_angle;
		double ** dN_dypTdpTdphi, ** SV_dN_dypTdpTdphi, ** cosine_iorder, ** sine_iorder;
		     
		//pair momentum
		double K_y, ch_K_y, sh_K_y;
		double current_K_phi, cos_cKphi, sin_cKphi;
		double beta_perp, beta_l;
		double * K_T, * K_phi, * K_phi_weight;
		    
		//spatial rapidity grid
		double * eta_s, * ch_eta_s, * sh_eta_s, * eta_s_weight;

		//points and weights for resonance integrals
		double * zeta_pts, * v_pts, ** s_pts, * NEW_s_pts;
		double * zeta_wts, * v_wts, ** s_wts, * NEW_s_wts;

		//some arrays to save unnecessary multiple calculations for resonances
		//use these for n_body = 2
		double VEC_n2_spt, VEC_n2_pstar, VEC_n2_Estar, VEC_n2_psBmT, VEC_n2_DeltaY, VEC_n2_Yp, VEC_n2_Ym;
		double * VEC_n2_P_Y, * VEC_n2_MTbar, * VEC_n2_DeltaMT, * VEC_n2_MTp, * VEC_n2_MTm;
		double ** VEC_n2_MT, ** VEC_n2_PPhi_tilde, ** VEC_n2_PPhi_tildeFLIP, ** VEC_n2_PT;
		double *** VEC_n2_Pp, *** VEC_n2_Pm, *** VEC_n2_alpha, *** VEC_n2_alpha_m;
		//double ** VEC_n2_PpT, ** VEC_n2_Ppphi;
		double VEC_n2_s_factor;
		double * VEC_n2_v_factor;
		double ** VEC_n2_zeta_factor;
		double VEC_n2_g_s;
		//use these for n_body = 3
		double * VEC_pstar, * VEC_Estar, * VEC_DeltaY, * VEC_Yp, * VEC_Ym, * VEC_s_factor, * VEC_g_s;
		double ** VEC_P_Y, ** VEC_MTbar, ** VEC_DeltaMT, ** VEC_MTp, ** VEC_MTm, ** VEC_v_factor;
		double *** VEC_MT, *** VEC_PPhi_tilde, *** VEC_PPhi_tildeFLIP, *** VEC_PT, *** VEC_zeta_factor;
		double **** VEC_Pp, **** VEC_alpha, **** VEC_Pm, **** VEC_alpha_m;
		
		//array to hold momenta to be integrated over in resonance calculations
		//momentum_data**** Pplus;
		
		//Emission function
		vector<Emissionfunction_data>* Emissionfunction_ptr;
		int FO_length;
		int Emissionfunction_length;
		vector<Emissionfunction_data>* avgFOsurf_ptr;
		
		double spectra;
		
		double * q_out, * q_side, * q_long;
		
		//store correlation functions
		double * Correl_1D_out;
		double * Correl_1D_out_err;
		double * Correl_1D_side;
		double * Correl_1D_side_err;
		double * Correl_1D_long;
		double * Correl_1D_long_err;

		//source variances
		double **S_func;
		double **xs_S, **xo_S, **xl_S, **t_S;
		double **xs2_S, **xo2_S, **xl2_S, **t2_S;
		double **xo_xs_S, **xl_xs_S, **xs_t_S, **xo_xl_S, **xo_t_S, **xl_t_S;
		
		//HBT radii coefficients
		double **R2_side, **R2_out, **R2_long, **R2_outside, **R2_sidelong, **R2_outlong;
		double **R2_side_C, **R2_side_S;
		double **R2_out_C, **R2_out_S;
		double **R2_long_C, **R2_long_S;
		double **R2_outside_C, **R2_outside_S;
		double **R2_outlong_C, **R2_outlong_S;
		double **R2_sidelong_C, **R2_sidelong_S;

		double *** res_sign_info, *** res_log_info, *** res_moments_info;

		
		//miscellaneous
		ofstream * global_out_stream_ptr;
		int global_folderindex;
		string global_path;
		string global_runfolder;
		string global_resultsfolder_stem;
		string no_df_stem;
		int n_resonance, n_decay_channels;
		int n_body;

		//some private methods		
		bool particles_are_the_same(int idx1, int idx2);
		bool Search_for_similar_particle(int dc_idx, int * result);

	public:
		SourceVariances(particle_info* particle, particle_info* all_particles_in, int Nparticle,
				FO_surf* FOsurf_ptr, vector<int> chosen_resonances, int particle_idx, ofstream& myout);
		~SourceVariances();

		void Determine_plane_angle(FO_surf* FOsurf_ptr, int dc_idx, bool thermal_particles_only = false);
		void Analyze_sourcefunction(FO_surf* FOsurf_ptr);
		void Analyze_sourcefunction_V1(FO_surf* FOsurf_ptr);
		void Analyze_sourcefunction_V2(FO_surf* FOsurf_ptr);
		void Analyze_sourcefunction_V3(FO_surf* FOsurf_ptr);
		void Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx);
		bool fexists(const char *filename);

		void Set_dN_dypTdpTdphi_moments(FO_surf* FOsurf_ptr, int dc_idx);
		void Cal_dN_dypTdpTdphi(double** SP_p0, double** SP_px, double** SP_py, double** SP_pz, FO_surf* FOsurf_ptr);
		void Cal_dN_dypTdpTdphi_with_weights_polar(FO_surf* FOsurf_ptr, int local_pid);
		void Cal_dN_dypTdpTdphi_with_weights_polar_V2(FO_surf* FOsurf_ptr, int local_pid);
		void Cal_dN_dypTdpTdphi_interpolate_polar_grid(double* SP_pT, double* SP_pphi);
		double Emissionfunction(double p0, double px, double py, double pz, FO_surf* surf);
		double weight_function(double PK[], int weight_function_index);
		void Do_resonance_integrals(int iKT, int iKphi, int dc_idx);
		void get_rapidity_dependence(double * rap_indep_vector, double * rap_dep_vector, double rap_val);
		void Set_current_daughter_info(int dc_idx, int daughter_idx);
		void Set_current_particle_info(int dc_idx);
		bool Do_this_decay_channel(int dc_idx);
		bool Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid);
		void Get_spacetime_moments(FO_surf* FOsurf_ptr, int dc_idx);
		void Recycle_spacetime_moments();
		void Allocate_decay_channel_info();
		void Load_decay_channel_info(int dc_idx, double K_T_local, double K_phi_local);
		void Delete_decay_channel_info();
		void combine_sourcevariances(double * output, double * input, double * alpha_vec);
		void compute_rap_indep_spacetime_moments(FO_surf* FOsurf_ptr, int dc_idx, double KTres, double Kphires, double * rapidity_independent_y_of_r);
		void Compute_source_variances(int iKT, int iKphi);

		void Get_source_variances(int, int);
		void Calculate_R2_side(int, int);
		void Calculate_R2_out(int, int);
		void Calculate_R2_long(int, int);
		void Calculate_R2_outside(int, int);
		void Calculate_R2_sidelong(int, int);
		void Calculate_R2_outlong(int, int);
		void R2_Fourier_transform(int iKT, double plane_psi);

		//miscellaneous
		void Set_ofstream(ofstream& myout);
		void Set_path(string path);
		void Set_resultsfolder_stem(string usrdef_stem);
		void Set_runfolder(string runfolder);
		void Set_use_delta_f(bool usrdef_usedeltaf);
		void Set_particle_mass(double usrdef_particle_mass);
		void Set_current_FOsurf_ptr(FO_surf* FOsurf_ptr);
		double get_Q();
		double g(double s);
		inline void set_to_zero(double * array, int arraylength);
		void adaptive_simpson_integration(void (SourceVariances::*f) (double, double *), double a, double b, double * results);
		double S_direct(double r, double eta, double tau, double MT, double PT, double cos_phi_m_Phi);
		double place_in_range(double phi, double min, double max);
		void Get_current_decay_string(int dc_idx, string * decay_string);
		int lookup_resonance_idx_from_particle_id(int particle_id);
		static inline double lin_int(double x1, double one_by_x2_m_x1, double f1, double f2, double x);
		double Edndp3(double ptr, double phir, int local_pid, int wfi);
		void Edndp3(double ptr, double phir, double * results);

		// input and output function prototypes
		void Output_SVdN_dypTdpTdphi(int folderindex);
		void Output_SVdN_dypTdpT(int folderindex);
		void Output_dN_dypTdpTdphi(int folderindex);
		void Output_dN_dypTdpT(int folderindex);
		void Output_dN_dypTdpTdphi_grid(int folderindex, int dc_idx);
		void Output_all_dN_dypTdpTdphi(int folderindex);
		void Output_total_target_dN_dypTdpTdphi(int folderindex);
		void Output_results(int folderindex);
		void Readin_results(int folderindex);
		void Read_in_all_dN_dypTdpTdphi(int folderindex);
		void Output_chosen_resonances();

		//parameters that the user is free to define
		double plumberg_test_variable;
		bool use_delta_f;
		bool append_output;
		int n_events;
		int initial_event, currentfolderindex;
		bool read_in_all_dN_dypTdpTdphi, output_all_dN_dypTdpTdphi;
};

#endif

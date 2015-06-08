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

/*typedef struct
{
   double P0, Po, Ps, Pl;
}momentum_data;*/

typedef struct
{
   double* resonance_mass;
   double* resonance_Gamma;
   double* resonance_total_br;
   double** resonance_decay_masses;
}resonance_info;


class SourceVariances
{
	private:
		//particle information 
		string particle_name;
		double particle_mass;
		int particle_id;     //particle id
		double particle_sign;   //+/- 1 for Fermi/Bose statistics for baryon/meson
		double particle_gspin;  //particle degeneracy 
		double particle_mu;
		particle_info * all_particles;

		int n_zeta_pts;
		int n_v_pts;
		int n_s_pts;
		double v_min;
		double v_max;
		double zeta_min;
		double zeta_max;
		double s_min, s_max;
		
		//array to hold resonance info
		resonance_info resonances;
		int current_resonance_idx;
		double current_resonance_mass;
		double current_resonance_Gamma;
		double current_resonance_total_br;
		double* current_resonance_decay_masses;
		double* P_eval;
		
		//array to temporarily hold results of resonance SV integrations
		int n_weighting_functions;  //number of source variances to consider
		//double * source_variances_array;
		double **** integrated_spacetime_moments;
		double **** dN_dypTdpTdphi_moments;
		double **** ln_dN_dypTdpTdphi_moments;
		double **** sign_of_dN_dypTdpTdphi_moments;
	
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
		double m2, m3, Gamma, br;
		double* Pp;
		double* Pm;
		double * zvec;

		//SP momentum arrays for interpolation grid
		double* SPinterp1_px;
		double* SPinterp1_py;
		double*** SPinterp1_p0;
		double*** SPinterp1_pz;
		double* SPinterp2_pT;
		double* SPinterp2_pphi;
		double* sin_SPinterp2_pphi;
		double* cos_SPinterp2_pphi;
		double** SPinterp2_p0;
		double** SPinterp2_pz;

		FO_surf* current_FOsurf_ptr;
		//FO surface info that is constant - to save time
		double Tdec, Edec, Pdec, muRES, signRES, gRES;
		double S_prefactor;
		double ** surf_damu, ** surf_pimunu, ** surf_Bn_muS_muB, ** surf_geometry_pts, ** surf_particle_mu, ** surf_flow;
	
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
		double* dN_dydphi;
		double* pTdN_dydphi;
		double* plane_angle;
		double global_plane_psi;
		double mean_pT;
		     
		//pair momentum
		double K_y, ch_K_y, sh_K_y;
		double* K_T;
		double* K_phi;
		double* K_phi_weight;
		//double K0, Kx, Ky, Kz;
		//inelegant way of making K_phi accessible to resonance integrals
		//without passing it everywhere...
		double current_K_phi, cos_cKphi, sin_cKphi;
		double beta_perp, beta_l;
		    
		//spatial rapidity grid
		double* eta_s;
		double* ch_eta_s;
		double* sh_eta_s;
		double* eta_s_weight;

		//points and weights for resonance integrals
		double* zeta_pts;
		double* v_pts;
		double** s_pts;
		double* zeta_wts;
		double* v_wts;
		double** s_wts;

		//some arrays to save unnecessary multiple calculations for resonances
		//use these for n_body = 2
		double VEC_n2_pstar;
		double VEC_n2_Estar;
		double VEC_n2_DeltaY;
		double VEC_n2_Yp;
		double VEC_n2_Ym;
		double * VEC_n2_P_Y;
		double * VEC_n2_MTbar;
		double * VEC_n2_DeltaMT;
		double * VEC_n2_MTp;
		double * VEC_n2_MTm;
		double ** VEC_n2_MT;
		double ** VEC_n2_PPhi_tilde, ** VEC_n2_PPhi_tildeFLIP, ** VEC_n2_PT;
		double *** VEC_n2_Pp;
		//double ** VEC_n2_PpT, ** VEC_n2_Ppphi;
		double VEC_n2_s_factor;
		double * VEC_n2_v_factor;
		double ** VEC_n2_zeta_factor;
		double VEC_n2_g_s;
		//use these for n_body = 3
		double * VEC_pstar;
		double * VEC_Estar;
		double * VEC_DeltaY;
		double * VEC_Yp;
		double * VEC_Ym;
		double ** VEC_P_Y;
		double ** VEC_MTbar;
		double ** VEC_DeltaMT;
		double ** VEC_MTp;
		double ** VEC_MTm;
		double *** VEC_MT;
		double *** VEC_PPhi_tilde, *** VEC_PPhi_tildeFLIP, *** VEC_PT;
		double **** VEC_Pp;
		//double *** VEC_PpT, *** VEC_Ppphi;
		double * VEC_s_factor;
		double ** VEC_v_factor;
		double *** VEC_zeta_factor;
		//double * VEC_tau_factor;
		double * VEC_g_s;
		
		//array to hold momenta to be integrated over in resonance calculations
		//momentum_data**** Pplus;
		
		//Emission function
		vector<Emissionfunction_data>* Emissionfunction_ptr;
		int FO_length;
		int Emissionfunction_length;
		vector<Emissionfunction_data>* avgFOsurf_ptr;
		
		double spectra;
		
		//source variances
		double **S_func;
		double **xs_S;
		double **xs2_S;
		
		//HBT radii coefficients
		double **R2_side;
		
		//miscellaneous
		ofstream * global_out_stream_ptr;
		int global_folderindex;
		string global_path;
		string global_runfolder;
		string global_resultsfolder_stem;
		string no_df_stem;
		int n_resonance;
		int n_body;

		//some private methods		
		bool particles_are_the_same(int idx1, int idx2);

	public:
		SourceVariances(particle_info* particle, particle_info* all_particles_in);
		~SourceVariances();

		void Determine_plane_angle(FO_surf* FOsurf_ptr);
		void Determine_plane_angle_check(FO_surf* FOsurf_ptr);
		void Analyze_sourcefunction(FO_surf* FOsurf_ptr);
		void Analyze_sourcefunction_check(FO_surf* FOsurf_ptr);
		void Reset_EmissionData();
		void Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx);
		void SetEmissionData(FO_surf* FOsurf_ptr, double K_T_local, double K_phi_local);
		bool fexists(const char *filename);

		void Set_dN_dypTdpTdphi_moments(FO_surf* FOsurf_ptr, int reso_idx);
		void Cal_dN_dypTdpTdphi(double** SP_p0, double** SP_px, double** SP_py, double** SP_pz, FO_surf* FOsurf_ptr);
		void Cal_dN_dypTdpTdphi_with_weights_cartesian(FO_surf* FOsurf_ptr, int reso_idx);
		void Cal_dN_dypTdpTdphi_with_weights_polar(FO_surf* FOsurf_ptr, int reso_idx);
		void Cal_dN_dypTdpTdphi_with_weights_polar_NEW(FO_surf* FOsurf_ptr, int reso_idx);
		double loop_over_FO_surface(FO_surf* FOsurf_ptr, double p0, double px, double py, double pz, double mu);
		void Cal_dN_dypTdpTdphi_interpolate_cartesian_grid(double** SP_px, double** SP_py);
		void Cal_dN_dypTdpTdphi_interpolate_polar_grid(double* SP_pT, double* SP_pphi);
		double Emissionfunction(double p0, double px, double py, double pz, FO_surf* surf);
		double weight_function(double PK[], int weight_function_index);
		//void Do_resonance_integrals(FO_surf* FOsurf_ptr, double K_T_local, double K_phi_local, int reso_idx);
		void Do_resonance_integrals(FO_surf* FOsurf_ptr, int iKT, int iKphi, int reso_idx);
		void Load_resonance_info(int reso_idx, double K_T_local, double K_phi_local);
		void Update_source_variances(int iKT, int iKphi, int reso_idx);

		void Get_source_variances(int, int);
		void Calculate_R2_side(int, int);

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
		void do_all_integrals(int iKT, int iKphi, int reso_idx);
		void set_to_zero(double * array, int arraylength);
		//void set_surfarrays();

		// input and output function prototypes
		void Output_SVdN_dypTdpTdphi(int folderindex);
		void Output_SVdN_dypTdpT(int folderindex);

		//parameters that the user is free to define
		double plumberg_test_variable;
		bool use_delta_f;
		bool append_output;
		int n_events;
		int initial_event;
};

#endif

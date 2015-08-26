#ifndef SV_PARAMETERS_H
#define SV_PARAMETERS_H

#include<string>
#include<sstream>

using namespace std;

#define SYMMETRIC_PT_PTS 		0		// chooses whether or not to use gaussian spacing or symmetric spacing for pt points
#define INTERPOLATION_KIND		0		// selects kind of interpolation
							// 0 - linear
							// 1 - cubic
							// 2 - polynomial
#define UNIFORM_SPACING			false		// specifies uniform or non-uniform grid spacing for interpolation
#define ASSUME_ETA_SYMMETRIC 		1		// 1 means integrate only over eta_s = 0..eta_s_max, and multiply by 2 or 0 to speed up calculations
							// 0 means just integrate over given range of eta_s without worrying about symmetry
#define GROUPING_PARTICLES 		0		// set to 1 to perform calculations for similar particles together
#define PARTICLE_DIFF_TOLERANCE 	0.00		// particles with mass and chemical potential (for each FO-cell) difference less than this value
							// will be considered to be identical (b/c Cooper-Frye)
#define USE_PLANE_PSI_ORDER 		0		// specifies whether to do HBT relative to flow-plane angle,
							// and at what order: 0 - use plane_psi = 0.0, !0 - use flow-plane angle at given order
#define TRUNCATE_COOPER_FRYE		false		// ignore contributions to CF integral which are extremely small --> speeds up code by factor of 3-4
#define VERBOSE 			1		// specifies level of output - 0 is lowest (no output)
#define DEBUG				false		// flag for output of debugging statements
#define EXIT_EARLY			1		// 1 decides whether to exit early from adaptive simpson's method routine
							// 0 forces routine to exit only when number of iterations j = jmax = 4
#define RECYCLE_ST_MOMENTS		true		// decides whether to recompute spacetime moments for each decay channel
							// only need this switch to see how much time this saves and make sure it's bug-free
#define SPACETIME_MOMENTS_ONLY		false		// duh
#define CHECK_FOR_LIFETIME		true		// true means skip particles which are too long-lived
#define DO_ALL_DECAY_CHANNELS		false		// duh
#define INCLUDE_SOURCE_VARIANCES	true		// false means do spectra only
#define USE_INTERP_ALT			true		// uses routine from iS.e code for comparison

const double hbarC=0.197327053;  //GeV*fm
const double twopi = 2.*M_PI;
const double MeVToGeV = 0.001;

//particle information
const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

//spatial rapidity information
const int eta_s_npts = 15;
const double eta_s_i = 0.0;
//const int eta_s_npts = 40;
//const double eta_s_i = -5.0;
const double eta_s_f = 4.0;

//relative momentum information
const int qnpts = 51;
const double delta_q = 0.02;
//const double init_q = -0.3;
const double init_q = 0.;

//single particle spectra info
const int n_SP_pT = 15;
//const int n_SP_pT = 101;
const int n_SP_pphi = 48;
const double SP_pT_min = 0.0;
const double SP_pT_max = 3.0;
//const double SP_pT_max = 5.0;

//parameters for interpolation grid
//  - polar
//const int n_interp_pT_pts = 77;
//const int n_interp_pphi_pts = 11;
const int n_interp_pT_pts = 15;
const int n_interp_pphi_pts = 48;
const double interp_pT_min = 0.01;
const double interp_pphi_min = 0.0;
const double interp_pT_max = 4.05;
const double interp_pphi_max = 2.*M_PI;
const double Del2_pT = (interp_pT_max - interp_pT_min) / (double)(n_interp_pT_pts-1);
const double Del2_pphi = (interp_pphi_max - interp_pphi_min) / (double)(n_interp_pphi_pts-1);

//correlation function info
const int corrfuncdim = 1;
const bool lambdaflag = true;

//pair momentum info
const int n_localp_T = 14;
const double localp_T_min = 0.05;
const double localp_T_max = 0.7;
const int n_localp_phi = 51;
const double localp_phi_min = 0.0;
const double localp_phi_max = 2*M_PI;

const int n_order = 1;
//const int n_order = 10;

const string path = "/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.08/NEW_TDEP_V4/NEW_TDEP_V4_results-avg-";
const string runfolder = "/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.08/NEW_TDEP_V4";

//const double tol = 1e-15;  //tolarence
const double tol = 0.0;  //tolerance
const int flagneg = 0;		//neglect all points that are negative
				//choose flagneg == 0 to agree with iS.e
				//choose flagneg == 1 for the real world

const int MCint_calls = 5000;  //# of calls for monte carlo integration

const size_t fit_max_iterations = 1000;  // stop at this point if not converged 
const double fit_tolarence = 1e-6;

//misc. resonance info
const double max_lifetime = 100.;	// fm/c

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

#endif

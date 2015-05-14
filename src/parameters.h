#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
#include<sstream>
using namespace std;

#define SYMMETRIC_PT_PTS 0		// chooses whether or not to use gaussian spacing or symmetric spacing for pt points
#define INTERPOLATION_FORMAT 2		// 0 - no interpolation, calculate everything exactly (bad idea)
					// 1 - Cartesian grid spacing in px - py
					// 2 - polar grid spacing in pT - pphi
#define UNIFORM_SPACING true		// specifies uniform or non-uniform grid spacing for interpolation
#define ASSUME_ETA_SYMMETRIC 1		// 1 means integrate only over eta_s = 0..eta_s_max, and multiply by 2 or 0 to speed up calculations
					// 0 means just integrate over given range of eta_s without worrying about symmetry

const double hbarC=0.197327053;  //GeV*fm

//particle information
const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

//spatial rapidity information
const int eta_s_npts = 20;
const double eta_s_i = 0.0;
//const int eta_s_npts = 40;
//const double eta_s_i = -4.0;
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
//  - Cartesian
const int n_interp1_px_pts = 51;
const int n_interp1_py_pts = 51;
const double interp1_px_min = 0.01;
const double interp1_py_min = -0.05;
const double interp1_px_max = 0.21;
const double interp1_py_max = 0.05;
//const double interp1_px_min = -0.1;
//const double interp1_py_min = -0.1;
//const double interp1_px_max = 0.1;
//const double interp1_py_max = 0.1;
const double Del1_x = (interp1_px_max - interp1_px_min) / (double)(n_interp1_px_pts-1);
const double Del1_y = (interp1_py_max - interp1_py_min) / (double)(n_interp1_py_pts-1);
//  - polar
const int n_interp2_pT_pts = 31;
const int n_interp2_pphi_pts = 24;
const double interp2_pT_min = 0.01;
const double interp2_pphi_min = 0.0;
const double interp2_pT_max = 3.0;
const double interp2_pphi_max = 2.*M_PI;
//const double interp2_px_min = -0.1;
//const double interp2_py_min = -0.1;
//const double interp2_px_max = 0.1;
//const double interp2_py_max = 0.1;
const double Del2_pT = (interp2_pT_max - interp2_pT_min) / (double)(n_interp2_pT_pts-1);
const double Del2_pphi = (interp2_pphi_max - interp2_pphi_min) / (double)(n_interp2_pphi_pts-1);

//correlation function info
const int corrfuncdim = 1;
const bool lambdaflag = true;

//pair momentum info
const int n_localp_T = 21;
const double localp_T_min = 0.0;
const double localp_T_max = 2.0;
const int n_localp_phi = 51;
const double localp_phi_min = 0.0;
const double localp_phi_max = 2*M_PI;

const int n_order = 4;
//const int n_order = 10;

//const string path = "/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.08/NEW_TDEP_V4/NEW_TDEP_V4_results-";
//const string runfolder = "/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.08/NEW_TDEP_V4";

//const string path = "/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.20/results-";
//const string runfolder = "/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.20";

const string path = "/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.08/NEW_TDEP_V4/NEW_TDEP_V4_results-avg-";
const string runfolder = "/home/plumberg.1/HBTwidths_viscosity_dependence/RESULTS/RESULTS_etaBYs_0.08/NEW_TDEP_V4";

const double tol = 1e-15;  //tolarence
const int flagneg = 1;     //neglect all points that are negative

const int MCint_calls = 5000;  //# of calls for monte carlo integration

const size_t fit_max_iterations = 1000;  // stop at this point if not converged 
const double fit_tolarence = 1e-6;

//const int n_events = 2;
//const int initial_event = 1;

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

//needed some extra variables for check.cpp
const double Rad = 3.5;  //fm

#endif

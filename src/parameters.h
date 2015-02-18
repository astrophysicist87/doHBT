#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
#include<sstream>
using namespace std;

const double hbarC=0.197327053;  //GeV*fm

//particle information
const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

//spatial rapidity information
const int eta_s_npts = 40;
const double eta_s_i = 0.0;
const double eta_s_f = 4.0;

//relative momentum information
const int qnpts = 51;
const double delta_q = 0.02;
//const double init_q = -0.3;
const double init_q = 0.;

//single particle spectra info
const int n_SP_pT = 15;
const int n_SP_pphi = 48;
const double SP_pT_min = 0.0;
const double SP_pT_max = 3.0;

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

const int n_events = 10;
const int initial_event = 1;

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

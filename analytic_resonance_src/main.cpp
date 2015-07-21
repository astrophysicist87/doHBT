#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <string>
#include <time.h>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <limits>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_matrix.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_linalg.h>

using namespace std;

#include "../analytic_resonance_class_src/analytic_resonance.h"
#include "../analytic_resonance_class_src/gauss.h"
//#include "CMDreader.h"

bool output_to_screen = true;
bool output_to_file = true;

size_t n;
string input_filename;
string input_filepath;

int main(int argc, char *argv[])
{
	/*int success = Read_in_parameters_from_CMD(argc, argv);
	if (success > 0)
	{
		cerr << "Problem reading in command-line parameters!" << endl;
		exit(1);
	}
	if (output_to_screen)
	{
		cout << "Initial parameters are:" << endl
			<< "eps_2_bar = " << eps_2_bar << endl
			<< "v_2_bar = " << v_2_bar << endl
			<< "psi_2_bar = " << psi_2_bar << endl
			<< "eta_f = " << eta_f << endl
			<< "Rad = " << Rad << endl
			<< "T0 = " << T0 << endl;
	}*/
	
/**************************************************************************/
/*************************MAIN PART OF THE PROGRAM*************************/
/**************************************************************************/

	gsl_set_error_handler_off();
	
	AnalyticResonance ar;
	
	//set up stuff for integrations
	ar.Create_integrations_points_and_weights();
	ar.Set_xi_integration_points();

	//get spectra
	//Compute_direct_pion_spectra_OLD();
	/*for (int ires = 0; ires <= 0; ires++)
	for (int iKT = 0; iKT < n_K_T; iKT++)
	{
		cout << "ires = " << ires << ", K_T = " << K_T[iKT] << ":   "  << Y_integrated_direct_resonance_spectra(K_T[iKT], ires) << endl;
		Direct_contributions_to_Y_integrated_pion_spectra(K_T[iKT], 0.0);
	}
	if (1) return (0);*/
	/*if (USE_INTERPOLATION)
		Compute_direct_resonance_spectra();
	
	for (int iKT = 0; iKT < n_K_T; iKT++)
		Resonance_decay_contributions_to_pion_spectra(K_T[iKT], 0.0, 0.0, 1);*/
	double temp = 0.0;
	for (int iKT = 0; iKT < n_K_T; iKT++)
	{
		//cout << sqrt(K_T[iKT]*K_T[iKT] + m_pion*m_pion) - m_pion << "   " << /*1000.*2.*sqrt(2.*M_PI)**/Y_integrated_direct_resonance_spectra(K_T[iKT], 0) << endl;
		temp = ar.Direct_contributions_to_pion_spectra(K_T[iKT], 0.0, 0.00342436);
	}
	if (1) return (0);
		//Direct_contributions_to_Y_integrated_pion_spectra(K_T[iKT], 0.0);
	//if (USE_INTERPOLATION)
	//	Compute_Y_integrated_direct_resonance_spectra();
	for (int ires = 0; ires < n_resonance; ires++)
	{
	cerr << "Doing resonance #" << ires + 1 << endl;
	for (int iKT = 0; iKT < n_K_T; iKT++)
		temp = ar.Resonance_decay_contributions_to_Y_integrated_pion_spectra(K_T[iKT], ires + 1);
	}

/**************************************************************************/
/***********************END MAIN PART OF THE PROGRAM***********************/
/**************************************************************************/
	return (0);
}

//End of file

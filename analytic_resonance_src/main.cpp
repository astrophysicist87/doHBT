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

#include "analytic_resonance.h"
#include "read_in_resonances.h"
#include "gauss.h"
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

	//set up stuff for integrations
	Create_integrations_points_and_weights();

	//get spectra
	//Compute_direct_pion_spectra_OLD();
	for (int iKT = 0; iKT < n_K_T; iKT++)
		Direct_contributions_to_pion_spectra(K_T[iKT], 0.0, 0.0);
	//Compute_direct_resonance_spectra();

/**************************************************************************/
/***********************END MAIN PART OF THE PROGRAM***********************/
/**************************************************************************/
	return (0);
}

/**************************************************************************/


//End of file

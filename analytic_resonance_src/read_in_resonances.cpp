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

#include "read_in_resonances.h"

int read_in_resonances(resonance_info * resonances)
{
	int number_of_resonances = 0;
	string resonancefilename = "/home/plumberg.1/HBTPlumberg/EOS/temporary_resonance_data.dat";
	ifstream resonanceinput (resonancefilename.c_str());
	resonanceinput >> number_of_resonances;
	resonances->resonance_mass = new double [number_of_resonances];
	resonances->resonance_Gamma = new double [number_of_resonances];
	resonances->resonance_total_br = new double [number_of_resonances];
	resonances->resonance_mu = new double [number_of_resonances];
	resonances->resonance_gspin = new double [number_of_resonances];
	resonances->resonance_sign = new int [number_of_resonances];
	resonances->resonance_decay_masses = new double* [number_of_resonances];
	for (int ir=0; ir<number_of_resonances; ir++)
	{
		resonances->resonance_decay_masses[ir] = new double [2];
		resonances->resonance_decay_masses[ir][0] = 0.0;
		resonances->resonance_decay_masses[ir][1] = 0.0;
		resonances->resonance_mu[ir] = 0.0;
		resonances->resonance_gspin[ir] = 1.0;	//actual g's have been absorbed into definitions of br
		//resonances->resonance_sign[ir] = 1;	//not quite right
	}
	int row_index = 0;
	resonanceinput >> row_index;
	while (!resonanceinput.eof() && row_index != 0)
	{
		//note that we have to convert given table values to GeV
		resonanceinput >> resonances->resonance_mass[row_index-1];
		resonanceinput >> resonances->resonance_decay_masses[row_index-1][0];
		resonanceinput >> resonances->resonance_decay_masses[row_index-1][1];
		resonanceinput >> resonances->resonance_Gamma[row_index-1];
		resonanceinput >> resonances->resonance_total_br[row_index-1];
		resonanceinput >> resonances->resonance_sign[row_index-1];
		//if (DEBUG)
		//	cerr << "Made it through row_index = " << row_index << endl;
		resonanceinput >> row_index;
	}
	resonanceinput.close();
	for (int ir=0; ir<number_of_resonances; ir++)
	{
		resonances->resonance_mass[ir] *= MeVToGeV;
		resonances->resonance_decay_masses[ir][0] *= MeVToGeV;
		resonances->resonance_decay_masses[ir][1] *= MeVToGeV;
		resonances->resonance_Gamma[ir] *= MeVToGeV;
	}
	
	return (number_of_resonances);
}

//End of file
